# TO DO: 
# edit input errors (does is automatically change input, then error message doesnt show)
# maybe make error messages red and bigger?
# make it so when the inputs are empty, it doenst crash, just has an error
# add a submit button


#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(tidyverse)
library(extraDistr)
library(DT)  # Load DT for data tables

source('../../functions.R', encoding = 'UTF-8')

# Define UI for application that calculates and displays confidence intervals
ui <- fluidPage(
  
  # Application title
  titlePanel("Confidence Intervals for the Negative Hypergeometric Distribution"),
  
  # Sidebar with inputs
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "parameter_of_interest", 
                  label = "Parameter of Interest", 
                  choices = list("N: Total number of items" = "N", 
                                 "M: Number of successes in the population" = "M")),
      
      selectInput(inputId = "procedure", 
                  label = "Procedure", 
                  choices = c("Normal Approximation / Large Sample (MLE)", 
                              "Normal Approximation / Large Sample (Unbiased)",
                              "Analog to the Clopper-Pearson", 
                              "Modified Sterne", 
                              "Crow & Gardner", 
                              "Bryne and Kabaila", 
                              "Blaker", 
                              "CMC")),
      
      numericInput(inputId = "confidence_level", 
                   label = "Confidence Level (%)", 
                   value = 95, 
                   min = 0, 
                   max = 100, 
                   step = 1),
      
      numericInput(inputId = "m", 
                   label = "Fixed number of successes to be observed (m)", 
                   value = 1, 
                   min = 1, 
                   step = 1),
      
      conditionalPanel(
        condition = "input.parameter_of_interest == 'N'",
        numericInput(inputId = "M", 
                     label = "Number of successes in the population (M)", 
                     value = 10, 
                     min = 1, 
                     step = 1)
      ),
      
      conditionalPanel(
        condition = "input.parameter_of_interest == 'M'",
        numericInput(inputId = "N", 
                     label = "Total number of items (N)", 
                     value = 20, 
                     min = 1, 
                     step = 1)
      ),
      
      numericInput(inputId = "x", 
                   label = "Observed value of x (Number of failures observed before the mth success)", 
                   value = 1, 
                   min = 0, 
                   step = 1),
      
      # Placeholder for error messages
      textOutput("error_message"),
      
      # Contact information
      HTML('<p style="font-size: x-small; margin-top: 20px;">
             <strong>Contact Information:</strong><br>
             Rachel Roggenkemper<br>
             rroggenk@calpoly.edu
           </p>')
    ),
    
    # Main panel to show the text and data table
    mainPanel(
      verbatimTextOutput("text_output"),  # For the descriptive text
      DTOutput("result_table")  # For the interactive data table
    )
  )
)

# Define server logic to calculate and display the confidence intervals
server <- function(input, output, session) {
  
  observe({
    # Validate confidence level
    if (input$confidence_level < 0 || input$confidence_level > 100) {
      #updateNumericInput(session, "confidence_level", value = 95)
      output$error_message <- renderText("Error: Confidence Level must be between 0 and 100.")
    } 
    else if (input$confidence_level != floor(input$confidence_level)) {
      #updateNumericInput(session, "confidence_level", value = 95)
      output$error_message <- renderText("Error: Confidence Level must be an integer.")
    } 
    else {
      output$error_message <- renderText("")
    }
  })
  
  observeEvent(input$m, {
    # Validate m
    if (input$m < 1) {
      #updateNumericInput(session, "m", value = 1)
      output$error_message <- renderText("Error: m must be an integer greater than or equal to 1.")
    } 
    else if (input$m != floor(input$m)) {
      #updateNumericInput(session, "m", value = 1)
      output$error_message <- renderText("Error: m must be an integer.")
    } 
    else if (input$m > input$M) {
      #updateNumericInput(session, "m", value = input$M)
      output$error_message <- renderText("Error: m must be less than or equal to M.")
    } 
    else {
      output$error_message <- renderText("")
    }
  })
  
  observeEvent(input$N, {
    # Validate N
    if (input$N < 1) {
      #updateNumericInput(session, "N", value = 20)
      output$error_message <- renderText("Error: N must be an integer greater than or equal to 1.")
    } 
    else if (input$N != floor(input$N)) {
      #updateNumericInput(session, "N", value = 20)
      output$error_message <- renderText("Error: N must be an integer.")
    } 
    else if (input$M > input$N) {
      #updateNumericInput(session, "M", value = input$N)
      output$error_message <- renderText("Error: M must be less than or equal to N.")
    } 
    else {
      output$error_message <- renderText("")
    }
  })
  
  observeEvent(input$M, {
    # Validate M
    if (input$M < 0) {
      #updateNumericInput(session, "M", value = 0)
      output$error_message <- renderText("Error: M must be an integer greater than or equal to 0.")
    } 
    else if (input$M != floor(input$M)) {
      #updateNumericInput(session, "M", value = 0)
      output$error_message <- renderText("Error: M must be an integer.")
    } 
    else if (input$M > input$N) {
      #updateNumericInput(session, "M", value = input$N)
      output$error_message <- renderText("Error: M must be less than or equal to N.")
    } 
    else if (input$m > input$M) {
      #updateNumericInput(session, "m", value = input$M)
      output$error_message <- renderText("Error: m must be less than or equal to M.")
    } 
    else {
      output$error_message <- renderText("")
    }
  })
  
  observeEvent(input$x, {
    # Validate x
    if (input$x < 0) {
      #updateNumericInput(session, "x", value = 1)
      output$error_message <- renderText("Error: x must be an integer greater than or equal to 0.")
    } 
    else if (input$x != floor(input$x)) {
      #updateNumericInput(session, "x", value = 1)
      output$error_message <- renderText("Error: x must be an integer.")
    } 
    else if (input$x > (input$N - input$M)) {
      #updateNumericInput(session, "x", value = input$N - input$m)
      output$error_message <- renderText("Error: Observed x must be less than or equal to N - m.")
    } 
    else {
      output$error_message <- renderText("")
    }
  })
  
  output$text_output <- renderText({
    # Handle only if the parameter of interest is "M"
    if (input$parameter_of_interest == "M") {
      # Convert confidence level to a decimal
      conf_level <- input$confidence_level / 100
      
      # Get N, m, and M based on the input
      N <- input$N
      m <- input$m
      M <- input$M
      
      # Text output with parameter details
      paste(
        "Parameter of Interest:", input$parameter_of_interest, "\n",
        "Procedure:", input$procedure, "\n",
        "Confidence Level (%):", input$confidence_level, "\n",
        "m (Number of successes to be observed):", input$m, "\n",
        "N (Total number of items):", N, "\n",
        "Observed x (Number of failures observed):", input$x, "\n"
      )
    } else {
      "Currently only handling cases where the parameter of interest is M."
    }
  })
  
  output$result_table <- renderDT({
    if (input$parameter_of_interest == "M") {
      # Convert confidence level to a decimal
      conf_level <- input$confidence_level / 100
      
      # Get N, m, and M based on the input
      N <- input$N
      m <- input$m
      M <- input$M
      
      # Calculate confidence interval based on the selected procedure
      result_df <- switch(input$procedure,
                          "Normal Approximation / Large Sample (MLE)" = CI_cov_prob_MLE(N = N, m = m, conf_level = conf_level),
                          "Normal Approximation / Large Sample (Unbiased)" = CI_cov_prob_unbiased(N = N, m = m, conf_level = conf_level),
                          "Analog to the Clopper-Pearson" = CI_cov_prob(N = N, m = m, conf_level = conf_level),
                          "Modified Sterne" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "MST"),
                          "Crow & Gardner" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "CG"),
                          "Bryne and Kabaila" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "BK"),
                          "Blaker" = blaker_ci(N = N, m = m, conf_level = conf_level),
                          "CMC" = cmc_ci(N = N, m = m, conf_level = conf_level)
      )
      
      # Extract the row corresponding to the observed x
      result_row <- result_df %>% filter(x == input$x)
      
      if (nrow(result_row) > 0) {
        # Display the result in a data table
        result_row %>%
          datatable(rownames = FALSE, options = list(pageLength = 5, dom = 't', autoWidth = TRUE))
      } else {
        data.frame(Message = "No result found for the given inputs.") %>%
          datatable(rownames = FALSE, options = list(pageLength = 1, dom = 't', autoWidth = TRUE))
      }
    } else {
      data.frame(Message = "Currently only handling cases where the parameter of interest is M.") %>%
        datatable(rownames = FALSE, options = list(pageLength = 1, dom = 't', autoWidth = TRUE))
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)



