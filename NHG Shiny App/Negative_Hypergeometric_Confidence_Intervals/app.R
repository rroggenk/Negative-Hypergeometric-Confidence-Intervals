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
library(shinybusy)
library(shinyjs)

source('../../functions.R', encoding = 'UTF-8')
source('../../functions_vec.R', encoding = 'UTF-8')

# Define UI for application that calculates and displays confidence intervals
ui <- fluidPage(
  useShinyjs(),
  add_busy_bar(color = "#1c88e5", height = "4px"),
  
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
                              "CMC",
                              "Blaker")),
      
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
      
      sliderInput(inputId="confidence_level",
                  label = "Confidence Level (%):",
                  min = 80, 
                  max = 99, 
                  value = 95, 
                  step = 1),
      
      conditionalPanel(
        condition = "input.parameter_of_interest == 'N'",
        sliderInput(inputId = "M", 
                    label = "Number of successes in the population (M):", 
                    min = 1, 
                    max = 50, 
                    value = 10, 
                    step = 1)
      ),
      
      conditionalPanel(
        condition = "input.parameter_of_interest == 'M'",
        sliderInput(inputId = "N", 
                    label = "Total number of items (N):", 
                    min = 1, 
                    max = 100, 
                    value = 20, 
                    step = 1)
      ),
      
      sliderInput(inputId = "m", 
                  label = "Fixed number of successes to be observed (m):", 
                  min = 1, 
                  max = 50, 
                  value = 5, 
                  step = 1),
      
      sliderInput(inputId = "x", 
                  label = "Observed value of x (Number of failures observed before the mth success):", 
                  min = 0, 
                  max = 50, 
                  value = 1, 
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
  
  # Dynamically update m based on M when parameter of interest is N
  observe({
    if (input$parameter_of_interest == "N") {
      max_m <- input$M
      updateSliderInput(session, "m", min = 1, max = max_m, value = 5)
    }
  })
  
  # Dynamically update m based on N when parameter of interest is M
  observe({
    if (input$parameter_of_interest == "M") {
      max_m <- input$N
      updateSliderInput(session, "m", min = 1, max = max_m, value = 5)
    }
  })
  
  # Dynamically update x based on N when parameter of interest is M
  observe({
    if (input$parameter_of_interest == "M") {
      max_x <- input$N 
      updateSliderInput(session, "x", min = 0, max = max_x, value = 5)
    }
  })
  
  
  # observe({
  #   # Validate confidence level
  #   if (input$confidence_level < 0 || input$confidence_level > 100) {
  #     validation_errors$message <- "Error: Confidence Level must be between 0 and 100."
  #   } else if (input$confidence_level != floor(input$confidence_level)) {
  #     validation_errors$message <- "Error: Confidence Level must be an integer."
  #   } else {
  #     validation_errors$message <- ""
  #   }
  # })
  # 
  # observeEvent(input$m, {
  #   # Validate m
  #   if (input$m < 1) {
  #     validation_errors$message <- "Error: m must be an integer greater than or equal to 1."
  #   } else if (input$m != floor(input$m)) {
  #     validation_errors$message <- "Error: m must be an integer."
  #   } else if (input$m > input$M) {
  #     validation_errors$message <- "Error: m must be less than or equal to M."
  #   } else {
  #     validation_errors$message <- ""
  #   }
  # })
  # 
  # observeEvent(input$N, {
  #   # Validate N
  #   if (input$N < 1) {
  #     validation_errors$message <- "Error: N must be an integer greater than or equal to 1."
  #   } else if (input$N != floor(input$N)) {
  #     validation_errors$message <- "Error: N must be an integer."
  #   } else if (input$M > input$N) {
  #     validation_errors$message <- "Error: M must be less than or equal to N."
  #   } else {
  #     validation_errors$message <- ""
  #   }
  # })
  # 
  # observeEvent(input$M, {
  #   # Validate M
  #   if (input$M < 0) {
  #     validation_errors$message <- "Error: M must be an integer greater than or equal to 0."
  #   } else if (input$M != floor(input$M)) {
  #     validation_errors$message <- "Error: M must be an integer."
  #   } else if (input$M > input$N) {
  #     validation_errors$message <- "Error: M must be less than or equal to N."
  #   } else if (input$m > input$M) {
  #     validation_errors$message <- "Error: m must be less than or equal to M."
  #   } else {
  #     validation_errors$message <- ""
  #   }
  # })
  # 
  # observeEvent(input$x, {
  #   # Validate x
  #   if (input$x < 0) {
  #     validation_errors$message <- "Error: x must be an integer greater than or equal to 0."
  #   } else if (input$x != floor(input$x)) {
  #     validation_errors$message <- "Error: x must be an integer."
  #   } else if (input$x > (input$N - input$M)) {
  #     validation_errors$message <- "Error: Observed x must be less than or equal to N - m."
  #   } else {
  #     validation_errors$message <- ""
  #   }
  # })
  
  # Reactive expression to run calculations after clicking submit button
  calculate_results <- eventReactive(input$submit, {
    if (validation_errors$message != "") {
      return(NULL)
    }
    
    disable("submit")  # Disable the Submit button
    on.exit(enable("submit"))  # Ensure re-enabling happens even on error
    
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
          "N (Total number of items):", N, "\n",
          "m (Number of successes to be observed):", input$m, "\n",
          "Observed x (Number of failures observed):", input$x, "\n"
        ),
        result_table = switch(input$procedure,
                              # "Normal Approximation / Large Sample (MLE)" = CI_cov_prob_MLE(N = N, m = m, conf_level = conf_level),
                              # "Normal Approximation / Large Sample (Unbiased)" = CI_cov_prob_unbiased(N = N, m = m, conf_level = conf_level),
                              # "Analog to the Clopper-Pearson" = CI_cov_prob(N = N, m = m, conf_level = conf_level),
                              # "Modified Sterne" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "MST"),
                              "Crow & Gardner" = minimal_cardinality_ci_vec(N = N, m = m, conf_level = conf_level, procedure = "CG"),
                              # "Bryne and Kabaila" = minimal_cardinality_ci(N = N, m = m, conf_level = conf_level, procedure = "BK"),
                              "Blaker" = blaker_ci(N = N, m = m, conf_level = conf_level),
                              "CMC" = cmc_ci_vec(N = N, m = m, conf_level = conf_level)
        ) %>% filter(x == input$x)
      )
    } 
    
    # Parameter of interest is "N"
    else if (input$parameter_of_interest == "N") {
      conf_level <- input$confidence_level / 100
      m <- input$m
      x <- input$x
      M <- input$M
      
      # Translate NHG to HG
      n_HG <- m + x
      M_HG <- M
      x_HG <- m
      
      # Get HG CI for N
      hg_ci_df <- CI.N(n = n_HG, M = M_HG, level = conf_level)
      hg_upper <- hg_ci_df[x_HG + 1, "upper"]
      
      # Handle Inf case
      if (is.infinite(hg_upper)) {
        max_N <- 10000
      } else {
        # Try different multipliers
        multiplier_list <- c(1.5, 3, 10)
        max_N <- NULL
        nhg_ci <- NULL
        for (mult in multiplier_list) {
          candidate_max_N <- ceiling(mult * hg_upper)
          nhg_ci <- switch(input$procedure,
                           "Crow & Gardner" = minimal_cardinality_ci_N_unkown_vec(M = M, m = m, conf_level = conf_level, max_N = candidate_max_N, procedure = "CG"),
                           "Blaker" = blaker_ci_N_unkown_vec(M = M, m = m, conf_level = conf_level, max_N = candidate_max_N),
                           "CMC" = cmc_ci_N_unkown_vec(M = M, m = m, conf_level = conf_level, max_N = candidate_max_N)
          )
          if (nrow(nhg_ci %>% filter(x == input$x)) > 0) {
            max_N <- candidate_max_N
            break
          }
        }
        # If all failed, fall back to 10000
        if (is.null(max_N)) {
          max_N <- 10000
        }
      }
      
      # Final computation with chosen max_N
      nhg_ci_final <- switch(input$procedure,
                             "Crow & Gardner" = minimal_cardinality_ci_N_unkown_vec(M = M, m = m, conf_level = conf_level, max_N = max_N, procedure = "CG"),
                             "Blaker" = blaker_ci_N_unkown_vec(M = M, m = m, conf_level = conf_level, max_N = max_N),
                             "CMC" = cmc_ci_N_unkown_vec(M = M, m = m, conf_level = conf_level, max_N = max_N)
      ) %>% filter(x == input$x)
      
      list(
        text_output = paste(
          "Parameter of Interest:", input$parameter_of_interest, "\n",
          "Procedure:", input$procedure, "\n",
          "Confidence Level (%):", input$confidence_level, "\n",
          "M (Total successes in population):", M, "\n",
          "m (Number of successes to be observed):", m, "\n",
          "Observed x (Number of failures observed):", x
        ),
        result_table = nhg_ci_final
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





