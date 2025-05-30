
library(shiny)
library(ggplot2)
library(mathjaxr)
library(truncnorm)
library(magrittr)
library(dplyr)
library(tidyr)

# Define UI for inputs and resulting plots
ui <- navbarPage(
  theme = shinythemes::shinytheme("united"),
    withMathJax(), # Initialize mathJax so the equation renders properly
   
    # app title
    # titlePanel("Generalized logistic curves and Gaussian processes"),
    # h4("Count data as functions of generalized logistic curves and GPs"),
  
    # Navbar tab 1: logistic curves
    tabPanel("Generalized logistic curves",
      sidebarLayout(
        sidebarPanel(
          h5("K = 0 (right asymptote)"),
          numericInput(inputId = "A",
                         label = "A (left asymptote)",
                        min = -10,
                        max = 0, 
                        step = 0.5,
                        value = -3),
          numericInput(inputId = "B",
                         label = "B (the growth rate when B>0)",
                         min = 0,
                         max = 2,
                         step = 0.05,
                         value = 0.1),
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
                         value = 50),
            helpText("$$\\beta(t) = A + \\frac{(K - A)}{(1 + e^{-B(t-M)})^{1/v}}$$"),
            helpText("$$TR = 100 \\cdot (1 - e^\\beta)$$"),
         ),
         mainPanel(
          # plot of beta over time
           plotOutput("betaPlot", height = 300, fill = F),
           
           # plot the resulting TR
           plotOutput("TRplot", height = 300)
         )
      )
    ),
    
    tabPanel("Gaussian processes (GPs)",
      h4("Squared exponential kernel. Season length of 150 days."),
      sidebarLayout(
        sidebarPanel(
          numericInput(inputId = "l",
                              label = "l (length scale)",
                              min = 0,
                              max = 500, 
                              step = 10,
                              value = 30),
           numericInput(inputId = "tau",
                              label = "tau (standard deviation of GP)",
                              min = 0,
                              max = 50,
                              step = 0.5,
                              value = 1),
            helpText("$$c(x_i, x_j) = \\tau^2 \\cdot e^{\\left(-\\frac{(x_i - x_j)^2}{2l^2}\\right)}$$"),
                 
        ),
        mainPanel(
           plotOutput("GPplot", height = 300, fill = F)
        )
      )
    ),
    
    tabPanel("Resulting count data",
      h3("Count data for agriculture treatment plots, over a season."),
      sidebarLayout(
        sidebarPanel(
          numericInput(inputId = "nLocs",
                       label = "Number of locations",
                       min = 0,
                       max = 50, 
                       step = 1,
                       value = 9),
          numericInput(inputId = "b0",
                       label = "b0 (expected log-count for control treatment)",
                       min = -5,
                       max = 6,
                       step = 0.5,
                       value = 1),
          numericInput(inputId = "b0SD",
                       label = "SD of b0 (SD associated with locations)",
                       min = 0,
                       max = 6,
                       step = 0.5,
                       value = 1.5),
          numericInput(inputId = "ASD",
                       label = "SD of A",
                       min = 0,
                       max = 6,
                       step = 0.1,
                       value = 0.5),
          numericInput(inputId = "BSD",
                       label = "SD of B",
                       min = 0,
                       max = 2,
                       step = 0.01,
                       value = 0.02),
          numericInput(inputId = "MSD",
                       label = "SD of M",
                       min = 0,
                       max = 100,
                       step = 10,
                       value = 10),
          numericInput(inputId = "lSD",
                       label = "SD of length GP parameter",
                       min = 0,
                       max = 100,
                       step = 1,
                       value = 10),
          numericInput(inputId = "tauSD",
                       label = "SD of tau GP parameter",
                       min = 0,
                       max = 20,
                       step = 0.1,
                       value = 0.5),
          numericInput(inputId = "negbinom",
                       label = "negative binomial extra variability",
                       min = 0,
                       max = 100,
                       step = 1,
                       value = 30)
        ),
        mainPanel(
          plotOutput("CountPlot", height = 600)
        )
      )
    )
    
)


server <- function(input, output) {
  
  curveDF <- reactive({
    time = 0:150
    K = 0
    C = 1
    curveDF = data.frame(time,
                    beta = input$A + ( (K - input$A)  / 
                                         ( (C + exp(-1 * input$B * (time - input$M))) ^ (1 / input$v) ) ))
    curveDF$TR = 100 * (1 - exp(curveDF$beta))
    curveDF
  })
  
  GPDF <- reactive({
    x = seq(1, 150, by = 7)
    d = abs(outer(x, x, "-"))
    Sigma_SE = input$tau^2 * exp(-d^2/(2*input$l^2))
    GP = as.vector(mvtnorm::rmvnorm(1,sigma=Sigma_SE))
    data.frame(x, GP)
  })

  output$betaPlot <- renderPlot({
        ggplot(curveDF(), aes(time, beta)) + 
          geom_point() +
          geom_line() + 
          labs(x = "time (days after install)",
              y = expression(beta),
              title = "Regression coefficient") +
          theme_bw(base_size = 22)
  })
    
  output$TRplot <- renderPlot({
      ggplot(curveDF(), aes(time, TR)) + 
        geom_point() +
        geom_line() + 
        ylim(0, 100) + 
        labs(x = "time (days after install)",
             y = "TR",
             title = "Resulting trap reduction (TR) over time") +
        theme_bw(base_size = 22)
  })
    
  output$GPplot <- renderPlot({
      ggplot(GPDF(), aes(x, GP)) + 
        geom_line() + 
        labs(x = "time (days after install)",
             y = "GP",
             title = "1-D Gaussian process") +
        theme_bw(base_size = 22)
  })
    
  output$CountPlot <- renderPlot({
    Treatment = c("CGP", "trt A")
    Location = LETTERS[1:input$nLocs]
    Trap = 1:4

    # simulate the GP:
    x = seq(1, 150, by = 7)
    d = abs(outer(x, x, "-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    l = rnorm(n = input$nLocs, mean = input$l, sd = input$lSD) # smaller l is more wiggle
    tau = rnorm(n = input$nLocs, mean = input$tau, sd = input$tauSD)
    for (i in 1:input$nLocs) {
      Sigma_SE = tau[i]^2 * exp(-d^2/(2*l[i]^2)) # squared exponential kernel
      GP = as.vector(mvtnorm::rmvnorm(1, sigma=Sigma_SE))
      if (all(GP<0)) GP = abs(GP)
      if (i == 1) {
        out = data.frame(DAT = x,
                         GP = as.vector(GP),
                         Location = Location[i])
      } else {
        out %<>%
          bind_rows(data.frame(DAT = x,
                               GP = as.vector(GP),
                               Location = Location[i]))
      }
    }
    #

    b0 = rnorm(n = input$nLocs, input$b0, sd = input$b0SD)
    A = rnorm(n = input$nLocs, input$A, sd = input$ASD)

    B = rtruncnorm(n = input$nLocs,a = 0, mean = input$B, sd = input$BSD)
    M = rtruncnorm(n = input$nLocs, a = 0, mean = input$M, sd = input$MSD)

    simdata <- expand_grid(Location, Treatment, Trap, DAT = x) %>%
      mutate(DaysOfCatch = 1,
             numTrt = ifelse(Treatment == "CGP", 0, 1)) %>%
      left_join(data.frame(b0, Location))  %>%
      left_join(data.frame(A, Location)) %>%
      left_join(data.frame(B, Location)) %>%
      left_join(data.frame(M, Location)) %>%
      left_join(out)

    simdata %<>%
      mutate(b1 = A * (1 - 1 / (1 + exp(-1 * B * (DAT - M) ))),
             muTmp = exp(b0 + b1*numTrt + GP),
             nYSB = rnbinom(n=nrow(simdata), mu = muTmp, size = input$negbinom))

    simdata %<>%
      mutate(nYSB = ifelse(is.na(nYSB), 0, nYSB))

    ggplot(simdata, aes(DAT, nYSB, color = Treatment)) +
      geom_point() +
      geom_smooth() + 
      facet_wrap(~Location, scales = "free") + 
      labs(x = "time (days after install)",
           y = "Moth counts per treatment per trap",
           title = "Simulated count data") +
      theme_bw()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
