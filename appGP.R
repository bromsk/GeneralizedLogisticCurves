
library(shiny)
library(ggplot2)
library(mathjaxr)

# Tab that simulates gaussian processes for our time frame
ui <- fluidPage(
    withMathJax(), # Initialize mathJax so the equation renders properly
  
    titlePanel("Gaussian processes (GP)"),
    # helpText("$$\\c(x_i, x_j) = \\tau^2 \\cdot e^{\\left(-\\frac{(x_i - x_j)^2}{2l^2}\\right)$$"),
    helpText("$$c(x_i, x_j) = \\tau^2 \\cdot e^{\\left(-\\frac{(x_i - x_j)^2}{2l^2}\\right)}$$"),
    
    h4("(Squared exponential kernel. Season length of 150 days.)"),
  
    # Sidebar for inputs
    sidebarLayout(
        sidebarPanel(
          numericInput(inputId = "l",
                         label = "l (length scale)",
                        min = 0,
                        max = 500, 
                        step = 10,
                        value = 10),
          numericInput(inputId = "tau",
                         label = "tau (standard deviation of GP)",
                         min = 0,
                         max = 50,
                         step = 1,
                         value = 5)
         ),

        # Show plots
        mainPanel(
          # plot of the variable over time
           plotOutput("GPplot", height = 300, fill = F)
        )
    )
)


server <- function(input, output) {

    output$GPplot <- renderPlot({
      x = seq(1, 150, by = 7)
      d = abs(outer(x, x, "-"))
      Sigma_SE = input$tau^2 * exp(-d^2/(2*input$l^2))
      GP = as.vector(mvtnorm::rmvnorm(1,sigma=Sigma_SE))
      
      # plot(x, GP, type = "l", xlab = "Days after install (DAI)")
      ggplot(data.frame(x, GP), aes(x, GP)) + 
        geom_line() + 
        labs(x = "time (days after install)",
             y = "GP",
             title = "1-D Gaussian process") +
        theme_bw(base_size = 22)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
