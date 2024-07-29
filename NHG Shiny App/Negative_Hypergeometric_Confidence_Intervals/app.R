#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Confidence Intervals for the Negative Hypergeometric Distribution"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("parameter_of_interest", 
                  "Parameter of Interest", 
                  choices = list("N: Total number of items" = "N", 
                                 "M: Number of successes in the population" = "M")),
      
      selectInput("procedure", 
                  "Procedure", 
                  choices = c("Normal Approximation / Large Sample", 
                              "Analog to the Clopper-Pearson", 
                              "Modified Sterne", 
                              "Crow & Gardner", 
                              "Bryne and Kabaila", 
                              "Blaker", 
                              "CMC")),
      
      numericInput("confidence_level", 
                   "Confidence Level (%)", 
                   value = 95, 
                   min = 0, 
                   max = 100, 
                   step = 1),
      
      numericInput("m", 
                   "Fixed number of successes to be observed (m)", 
                   value = 1, 
                   min = 1, 
                   step = 1),
      
      conditionalPanel(
        condition = "input.parameter_of_interest == 'N'",
        numericInput("M", 
                     "Number of successes in the population (M)", 
                     value = 0, 
                     min = 0, 
                     step = 1),
        
        numericInput("N", 
                     "Total number of items (N)", 
                     value = 1, 
                     min = 1, 
                     step = 1)
      ),
      
      conditionalPanel(
        condition = "input.parameter_of_interest == 'M'",
        numericInput("N", 
                     "Total number of items (N)", 
                     value = 1, 
                     min = 1, 
                     step = 1),
        
        numericInput("M", 
                     "Number of successes in the population (M)", 
                     value = 0, 
                     min = 0, 
                     step = 1)
      ),
      
      numericInput("x", 
                   "Observed value of x (Number of failures observed before the mth success)", 
                   value = 0, 
                   min = 0, 
                   step = 1),
      
      # Placeholder for error messages
      textOutput("error_message")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observe({
    # Validate confidence level
    if (input$confidence_level < 0 || input$confidence_level > 100) {
      updateNumericInput(session, "confidence_level", value = 95)
      output$error_message <- renderText("Error: Confidence Level must be between 0 and 100.")
    } else {
      output$error_message <- renderText("")
    }
    
    # Validate m
    if (input$m < 1) {
      updateNumericInput(session, "m", value = 1)
      output$error_message <- renderText("Error: m must be an integer greater than or equal to 1.")
    } else {
      output$error_message <- renderText("")
    }
    
    # Validate N
    if (input$N < 1) {
      updateNumericInput(session, "N", value = 1)
      output$error_message <- renderText("Error: N must be an integer greater than or equal to 1.")
    } else {
      output$error_message <- renderText("")
    }
    
    # Validate M
    if (input$M < 0) {
      updateNumericInput(session, "M", value = 0)
      output$error_message <- renderText("Error: M must be an integer greater than or equal to 0.")
    } else {
      output$error_message <- renderText("")
    }
    
    # Validate x
    if (input$x < 0) {
      updateNumericInput(session, "x", value = 0)
      output$error_message <- renderText("Error: x must be an integer greater than or equal to 0.")
    } else {
      output$error_message <- renderText("")
    }
  })
  
  output$distPlot <- renderPlot({
    # Generate bins based on input$bins from ui.R
    x <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # Draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white',
         xlab = 'Waiting time to next eruption (in mins)',
         main = 'Histogram of waiting times')
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

