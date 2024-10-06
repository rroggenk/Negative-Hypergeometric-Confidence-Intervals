###################################
#---------------------------------#
# Load Packages                   #
#---------------------------------#
###################################
library(shiny) 
library(DT) #fancy data table
library(rsconnect) #to deploy (update) app to shinyapp.io


###################################
#---------------------------------#
# Source R Scripts                #
#---------------------------------#
###################################

# Preliminary functions                                       
source('PreliminaryFunctions.R', encoding = 'UTF-8')

# Minimal Cardinality Procedures (CG, OC, BK)
source('MinCard.R', encoding = 'UTF-8')

# Minimal Cardinality Procedures (CG, OC, BK)
source('CMC.R', encoding = 'UTF-8')

###################################
#---------------------------------#
# Shiny App                       #
#---------------------------------#
###################################

ui <- fluidPage(
  
  # Application title
  titlePanel("Confidence Intervals for the Count Parameters in Binomial and Negative Binomial Experiments"),
  
  sidebarLayout(

    sidebarPanel(
      
      #h5 5th level of a hiearchy of decreasing subheadings;i.e.subsubsubsubheading
      #br() creates a line break

      div(
        selectInput(inputId="parameter", 
                    label="parameter of interest:", 
                    choices = list("binomial n, x=# of successes observed", "negative binomial r, x=# of failures observed", "negative binomial r, N=# of trials observed")),
        
        selectInput(inputId="procedure", 
                    label="procedure:", 
                    choices = list("CMC", "CG"),
                    width='40%'),
        
        #Remove minor tick marks
        tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
        
        #remove bold from labels
        tags$head(tags$style(HTML("label {font-weight:normal;}"))),
        
        sliderInput(inputId="level",
                    label = "confidence level %:",
                    min = 80, max = 99, value = 95, step=1),
        
        sliderInput(inputId="prob_success",
                    label = "probability of success (p):",
                    min = 0.01, max = 0.99, value = 0.80, step=0.01),

      
        # sliderInput(inputId="max_rv",
        #             label = "largest value of the RV desired:",
        #             min = 0, max = 200, value = 10, step=1),br(),
        
      
        numericInput(inputId="rv_value", 
                     label=NULL, 
                     value = 10, step=1, width='40%'),br(),
      
      style="font-size:90%"),

      #div(submitButton("Submit"),align="right"),br(), 
      
      div("*We take a negative binomial variable to be either the number of failures x 
          in successive independent Bernoulli trials until the rth success, or the 
          number of trials N = r + x until the rth success.", align="left", style = "font-size: 8pt"),br(),
      
      div(
      "Confidence intervals for sample size parameter n of the hypergeometric distribution
      can be obtained from,",
      a(href="https://mfschilling.shinyapps.io/hgci-success-count/",target="_blank","https://mfschilling.shinyapps.io/hgci-success-count/"),
      ", by swapping the roles of the sample size (n) and success count (M).",
      align="left", style = "font-size: 8pt")
      
      # ,div("maintained by",
      #     a(href="https://statistics.calpoly.edu/bret-holladay",target="_blank",
      #       "Bret Holladay"),align="center", style = "font-size: 8pt"), br(),
      # 
      # div(strong("Authors:"),align="center", style = "font-size: 8pt"),
      # div(a(href="https://statistics.calpoly.edu/bret-holladay",target="_blank",
      #       "Bret Holladay"), ", Cal Poly, San Luis Obispo",align="center", style = "font-size: 8pt"),
      # div(a(href="http://www.csun.edu/~hcmth031/",target="_blank",
      #       "Mark Schilling"),", California State University, Northridge",align="center", style = "font-size: 8pt")
    ),

    mainPanel(

      verbatimTextOutput("info"),
      #br(),
      verbatimTextOutput("results"),
      #tableOutput("results"),
      #br(),
    )
  )
)

server <- function(input, output, session) {

  observe({
    if(input$parameter=="negative binomial r, N=# of trials observed"){
      RV='N'
    }else{
      RV='X'
    }
    updateSliderInput(session, inputId="rv_value", label=paste("observed value of", RV,':'))
  })
  

  output$info <- renderText({ 
    paste(
      paste("parameter: ", input$parameter),
      paste("procedure: ", input$procedure),
      paste("level: ", input$level/100),
      paste("p: ", input$prob_success),
      sep="\n"
    )
  })
  
  #output$tbl <- DT::renderDataTable({
  #output$results<- renderTable(digits=0,{
  output$results<- renderPrint({
    parameter<-input$parameter
    procedure <- input$procedure
    conf.level <- input$level/100
    #max_rv <- input$max_rv
    max_rv <- input$rv_value
    p <- input$prob_success
    
    if(parameter=='binomial n, x=# of successes observed'){
      if(procedure=='CMC'){
        CI=CMC.binomn(p=p, x.max=max_rv, conf.level=conf.level)
        x=(CI$CMC)$x; lower=(CI$CMC)$LL; upper=(CI$CMC)$UL
        df=data.frame(x=x,lower=lower,upper=upper)
      }else if(procedure=='CG'){
        CI=MC.binomn(p=p,x.max=max_rv, conf.level=conf.level)
        x=(CI$CG)$x; lower=(CI$CG)$LL; upper=(CI$CG)$UL
        df=data.frame(x=x,lower=lower,upper=upper)
      }
    }else if(parameter=='negative binomial r, x=# of failures observed'){
      if(procedure=='CMC'){
        CI=CMC.negbinomr(p=p,x.max=max_rv, conf.level=conf.level)
        x=(CI$CMC)$x; lower=(CI$CMC)$LL; upper=(CI$CMC)$UL
        df=data.frame(x=x,lower=lower,upper=upper)  
      }else if(procedure=='CG'){
        CI=MC.negbinomr(p=p,x.max=max_rv, conf.level=conf.level)
        x=(CI$CG)$x; lower=(CI$CG)$LL; upper=(CI$CG)$UL
        df=data.frame(x=x,lower=lower,upper=upper)
      }
    }else if(parameter=='negative binomial r, N=# of trials observed'){
     if(procedure=='CMC'){
        CI=CMC.negbinomrN(p=p,y.max=max_rv, conf.level=conf.level)
        N=(CI$CMC)$N; lower=(CI$CMC)$LL; upper=(CI$CMC)$UL
        df=data.frame(N=N,lower=lower,upper=upper)  
      }else if(procedure=='CG'){
        CI=MC.negbinomrN(p=p,y.max=max_rv, conf.level=conf.level)
        N=(CI$CG)$N; lower=(CI$CG)$LL; upper=(CI$CG)$UL
        df=data.frame(N=N,lower=lower,upper=upper)
      }
    }
      #df[nrow(df),]
      df=df[nrow(df),]
      print(df, row.names = F)
    }) 
    #caption=paste0("parameter: ", input$parameter,"; procedure: ", input$procedure, "; level: ", input$level/100, ";   p: ", input$prob_success)
    
    # DT::datatable(df,rownames = FALSE,
    #               extensions = "Buttons",filter = "none",
    #               caption=caption,
    #               options = list(searching=FALSE,
    #                              lengthChange = TRUE,pageLength=10,lengthMenu = c(10, 25, 50, 75, 100, 200),
    #                              dom = "lrtB", 
    #                              buttons = c("copy", "csv", "excel", "pdf", "print")
    #               )
    # )
  # })
}

shinyApp(ui = ui, server = server)