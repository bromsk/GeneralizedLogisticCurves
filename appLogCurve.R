
library(shiny)
library(ggplot2)
library(mathjaxr)

# Define UI for inputs and resulting plots
ui <- fluidPage(
    withMathJax(), # Initialize mathJax so the equation renders properly
  
    titlePanel("Generalized logistic curves"),
    helpText("$$\\beta(t) = A + \\frac{(K - A)}{(1 + e^{-B(t-M)})^{1/v}}$$"),
    h4("(For a variable that starts negative and asymptotes at 0.)"),
  
    # Sidebar for inputs
    sidebarLayout(
        sidebarPanel(
          h5("K = 0 (right asymptote)"),
          numericInput(inputId = "A",
                         label = "A (left asymptote)",
                        min = -10,
                        max = 0, 
                        step = 0.5,
                        value = -5),
          numericInput(inputId = "B",
                         label = "B (the growth rate when B>0)",
                         min = 0,
                         max = 2,
                         step = 0.05,
                         value = 0.2),
            numericInput(inputId = "v",
                         label = "v (affects where the inflection point occurs, v > 0)",
                         min = 0,
                         max = 20,
                         step = 0.1,
                         value = 1),
            numericInput(inputId = "M",
                         label = "M (Relates to the time of when the curve begins)",
                         min = 0,
                         max = 150,
                         step = 1,
                         value = 25),
          
          helpText("$$TR = 100 \\cdot (1 - e^\\beta)$$"),
          
            
         ),

        # Show plots
        mainPanel(
          # plot of the variable over time
           plotOutput("betaPlot", height = 300, fill = F),
           
           # plot the resulting TR
           plotOutput("TRplot", height = 300)
        )
    )
)


server <- function(input, output) {
  
  DF <- reactive({
    time = 0:150
    K = 0
    C = 1
    DF = data.frame(time,
                    beta = input$A + ( (K - input$A)  / 
                                         ( (C + exp(-1 * input$B * (time - input$M))) ^ (1 / input$v) ) ))
    DF$TR = 100 * (1 - exp(DF$beta))
    DF
  })

    output$betaPlot <- renderPlot({
        ggplot(DF(), aes(time, beta)) + 
          geom_point() +
          geom_line() + 
          labs(x = "time (days after install)",
              y = expression(beta),
              title = "Regression coefficient") +
          theme_bw(base_size = 22)

    })
    
    output$TRplot <- renderPlot({
      ggplot(DF(), aes(time, TR)) + 
        geom_point() +
        geom_line() + 
        ylim(0, 100) + 
        labs(x = "time (days after install)",
             y = "TR",
             title = "Resulting trap reduction (TR) over time") +
        theme_bw(base_size = 22)
      
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
