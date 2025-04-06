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
library(data.table)
library(DT)  # Load DT for data tables

source('../../functions.R', encoding = 'UTF-8')
source('../../functions_vec.R', encoding = 'UTF-8')

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
                  choices = c("Crow & Gardner", 
                              "CMC")),
      
      # selectInput(inputId = "procedure", 
      #             label = "Procedure", 
      #             choices = c("Normal Approximation / Large Sample (MLE)", 
      #                         "Normal Approximation / Large Sample (Unbiased)",
      #                         "Analog to the Clopper-Pearson", 
      #                         "Modified Sterne", 
      #                         "Crow & Gardner", 
      #                         "Bryne and Kabaila", 
      #                         "Blaker", 
      #                         "CMC")),
      
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
      
      # Submit button
      actionButton(inputId = "submit", label = "Submit"),
      
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
  
  # Reactive values to store validation errors
  validation_errors <- reactiveValues(message = "")
  
  observe({
    # Validate confidence level
    if (input$confidence_level < 0 || input$confidence_level > 100) {
      validation_errors$message <- "Error: Confidence Level must be between 0 and 100."
    } else if (input$confidence_level != floor(input$confidence_level)) {
      validation_errors$message <- "Error: Confidence Level must be an integer."
    } else {
      validation_errors$message <- ""
    }
  })
  
  observeEvent(input$m, {
    # Validate m
    if (input$m < 1) {
      validation_errors$message <- "Error: m must be an integer greater than or equal to 1."
    } else if (input$m != floor(input$m)) {
      validation_errors$message <- "Error: m must be an integer."
    } else if (input$m > input$M) {
      validation_errors$message <- "Error: m must be less than or equal to M."
    } else {
      validation_errors$message <- ""
    }
  })
  
  observeEvent(input$N, {
    # Validate N
    if (input$N < 1) {
      validation_errors$message <- "Error: N must be an integer greater than or equal to 1."
    } else if (input$N != floor(input$N)) {
      validation_errors$message <- "Error: N must be an integer."
    } else if (input$M > input$N) {
      validation_errors$message <- "Error: M must be less than or equal to N."
    } else {
      validation_errors$message <- ""
    }
  })
  
  observeEvent(input$M, {
    # Validate M
    if (input$M < 0) {
      validation_errors$message <- "Error: M must be an integer greater than or equal to 0."
    } else if (input$M != floor(input$M)) {
      validation_errors$message <- "Error: M must be an integer."
    } else if (input$M > input$N) {
      validation_errors$message <- "Error: M must be less than or equal to N."
    } else if (input$m > input$M) {
      validation_errors$message <- "Error: m must be less than or equal to M."
    } else {
      validation_errors$message <- ""
    }
  })
  
  observeEvent(input$x, {
    # Validate x
    if (input$x < 0) {
      validation_errors$message <- "Error: x must be an integer greater than or equal to 0."
    } else if (input$x != floor(input$x)) {
      validation_errors$message <- "Error: x must be an integer."
    } else if (input$x > (input$N - input$M)) {
      validation_errors$message <- "Error: Observed x must be less than or equal to N - m."
    } else {
      validation_errors$message <- ""
    }
  })
  
  # Reactive expression to run calculations after clicking submit button
  calculate_results <- eventReactive(input$submit, {
    if (validation_errors$message != "") {
      return(NULL)
    }
    
    # Parameter of interest is "M"
    if (input$parameter_of_interest == "M") {
      # Convert confidence level to a decimal
      conf_level <- input$confidence_level / 100
      
      # Get N, m, and M based on the input
      N <- input$N
      m <- input$m
      M <- input$M
      
      # Return a list of results and text details
      list(
        text_output = paste(
          "Parameter of Interest:", input$parameter_of_interest, "\n",
          "Procedure:", input$procedure, "\n",
          "Confidence Level (%):", input$confidence_level, "\n",
          "m (Number of successes to be observed):", input$m, "\n",
          "N (Total number of items):", N, "\n",
          "Observed x (Number of failures observed):", input$x, "\n"
        ),
        result_table = switch(input$procedure,
                              # "Normal Approximation / Large Sample (MLE)" = CI_cov_prob_MLE(N = N, m = m, conf_level = conf_level),
                              # "Normal Approximation / Large Sample (Unbiased)" = CI_cov_prob_unbiased(N = N, m = m, conf_level = conf_level),
                              # "Analog to the Clopper-Pearson" = CI_cov_prob(N = N, m = m, conf_level = conf_level),
                              # "Modified Sterne" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "MST"),
                              "Crow & Gardner" = minimal_cardinality_ci_vec(N = N, m = m, conf_level = conf_level, procedure = "CG"),
                              # "Bryne and Kabaila" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "BK"),
                              # "Blaker" = blaker_ci(N = N, m = m, conf_level = conf_level),
                              "CMC" = cmc_ci_vec(N = N, m = m, conf_level = conf_level)
        ) %>% filter(x == input$x)
      )
    } 
    
    # Parameter of interest is "N"
    if (input$parameter_of_interest == "N") {
      # Convert confidence level to a decimal
      conf_level <- input$confidence_level / 100
      
      # Get N, m, and M based on the input
      N <- input$N
      m <- input$m
      M <- input$M
      
      # Return a list of results and text details
      list(
        text_output = paste(
          "Parameter of Interest:", input$parameter_of_interest, "\n",
          "Procedure:", input$procedure, "\n",
          "Confidence Level (%):", input$confidence_level, "\n",
          "m (Number of successes to be observed):", input$m, "\n",
          "M (Total successes in population):", M, "\n",
          "Observed x (Number of failures observed):", input$x, "\n"
        ),
        result_table = switch(input$procedure,
                              # "Analog to the Clopper-Pearson" = CI_Analog_CP_N_Unknown_vec(M = M, m = m, conf_level = conf_level, max_N = max_N),
                              # "Modified Sterne" = minimal_cardinality_ci(M = M, m = m, conf_level = conf_level, max_N = max_N, procedure = "MST"),
                              "Crow & Gardner" = minimal_cardinality_ci_vec(M = M, m = m, conf_level = conf_level, max_N = max_N, procedure = "CG"),
                              # "Bryne and Kabaila" = minimal_cardinality_ci(M = M, m = m, conf_level = conf_level, max_N = max_N, procedure = "BK"),
                              # "Blaker" = blaker_ci(M = M, m = m, conf_level = conf_level, max_N = max_N),
                              "CMC" = cmc_ci_vec(M = M, m = m, conf_level = conf_level, max_N = max_N)
        ) %>% filter(x == input$x)
      )
    } 
    
    
    else {
      list(
        text_output = "ERROR",
        result_table = NULL
      )
    }
  })
  
  # Render the text output
  output$text_output <- renderText({
    req(calculate_results())
    calculate_results()$text_output
  })
  
  # Render the result table
  output$result_table <- renderDT({
    req(calculate_results())
    if (is.null(calculate_results()$result_table) || nrow(calculate_results()$result_table) == 0) {
      data.frame(Message = "No result found for the given inputs.") %>%
        datatable(rownames = FALSE, options = list(pageLength = 1, dom = 't', autoWidth = TRUE))
    } else {
      calculate_results()$result_table %>%
        datatable(rownames = FALSE, options = list(pageLength = 5, dom = 't', autoWidth = TRUE))
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)





