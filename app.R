#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

n <- 10
levels <- c("Susceptible", "Infected", "Recovered", "Vaccinated")
cell_cols <- c("Susceptible" = "white", "Infected" = "red", "Recovered" = "gray", "Vaccinated" = "gray3")

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme="styles.css",
  titlePanel("Epidemic explorer"),
  sidebarLayout(
    sidebarPanel(
      actionButton("reSet","Reset model"),
      h3("Set up index cases before running model"),
      fluidRow(
        column(8,sliderInput("indexCaseN","number of index cases",min=1,max=10,step=1,value=1)),
        column(4,actionButton("indexCaseSet","Set index cases")),
      ),
      fluidRow(
        column(8,sliderInput("vaccinateN","number vaccinated",min=0,max=90,step=5,value=0)),
        column(4,actionButton("vaccinateSet","Vaccinate")),
      ),
      h3("The workflow for simulating one turn of the model"),
      fluidRow(
        column(8,sliderInput("beta","transmission rate (pathogen load)",min=0.5,max=8,step=0.5,value=2)),
        column(4,selectInput("betaDist", "Distribution for pathogen load", choices = c("delta", "Poisson")))
      ),
      fluidRow(
        column(4,actionButton("attemptsNSet","Set pathogen load")),
        column(8,numericInput("attemptsN", "Next transmission attempts",0,min=0,max=100,step=1))
      ),
      fluidRow(
        column(4,actionButton("transmitSet","Transmission attempts")),
        column(4,actionButton("removeSet","Remove infecteds")),
        column(4,actionButton("infectSet","New infections")),
      ),
      
      hr()
    ),
    mainPanel(
      fluidRow(
        column(8,plotOutput("gridPlot", height="480px")),
        column(4,
               h3("Model output values"),
               div(class = "info_row",
                   span(class = "info_label", "Susceptible"),
                   textOutput("S", inline=T)
               ),
               div(class = "info_row",
                   span(class = "info_label", "Infected"),
                   textOutput("I", inline=T)
               ),
               div(class = "info_row",
                   span(class = "info_label", "Removed"),
                   textOutput("R", inline=T)
               ),
               div(class = "info_row",
                   span(class = "info_label", "Vaccinated"),
                   textOutput("V", inline=T)
               ),
               div(class = "info_row",
                   span(class = "info_label", "Final epidemic size"),
                   textOutput("F", inline=T)
               ),
               div(class = "info_row",
                   span(class = "info_label", "Time"),
                   textOutput("Time", inline=T)
               ),
               div(class = "info_row",
                   span(class = "info_label", "Individuals in receipt of pathogen"),
                   textOutput("pX", inline=T)
               ),               
        )
      ),
      img(src='ioa_logo.png',style="width: 256px; align: left; margin-right: 2em"),
      "Darren Green (2026)",
      img(src='parasite_2.png',style="width: 64px; align: right; margin-left: 2em")
    )
  )
)

server <- function(input, output, session) {
  D <- reactiveValues()
  D$M <- matrix(factor(sample(levels, n^2, replace = TRUE), levels = levels), nrow = n)
  D$X <- matrix(sample(c(F, T), n^2, replace = TRUE), nrow = n)
  D$S<-0
  D$I<-0
  D$R<-0
  D$V<-0
  D$time <- 0
  output$S <- renderText({sum(D$M==levels[1])})
  output$I <- renderText({sum(D$M==levels[2])})
  output$R <- renderText({sum(D$M==levels[3])})
  output$V <- renderText({sum(D$M==levels[4])})
  output$F <- renderText({sum(D$M==levels[2]) + sum(D$M==levels[3])})
  output$pX <- renderText({sum(D$X)})
  output$Time <- renderText({D$time})
  output$gridPlot <- renderPlot({
    plot(NULL, xlim = c(0, n), ylim = c(0, n), asp = 1, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    for (row in 1:n) {
      for (col in 1:n) {
        val <- D$M[row, col]
        fill <- cell_cols[val]
        rect(col - 1, n - row, col, n - row + 1, col = fill, border = "black", lwd = 1)
        val <- D$X[row, col]
        if (val) points(col - 1 + 0.5, n - row + 0.5, col = "darkred", pch = 15)
      }
    }
  })
  observeEvent(input$reSet, {
    D$X[] <- F
    D$M[] <- levels[1]
    D$time <- 0
  })
  observeEvent(input$indexCaseSet, {
    if (input$indexCaseN==0) return
    ic <- sample(1:(n^2), input$indexCaseN)
    ix <- (ic-1) %/% n
    iy <- (ic-1) %% n
    D$M[cbind(ix+1,iy+1)] <- levels[2]
  })
  observeEvent(input$vaccinateSet, {
    if (input$vaccinateN==0) return
    ic <- sample(1:(n^2), input$vaccinateN)
    ix <- (ic-1) %/% n
    iy <- (ic-1) %% n
    D$M[cbind(ix+1,iy+1)] <- levels[4]
  })
  observeEvent(input$attemptsNSet, {
    xu <- sum(D$M==levels[2]) * input$beta
    x <- ifelse(input$betaDist == "delta", xu, rpois(1, xu))
    updateNumericInput(session, "attemptsN", value = x)
  })
  observeEvent(input$transmitSet, {
    if (input$attemptsN==0) return
    ic <- sample(1:(n^2), input$attemptsN, replace=T)
    ix <- (ic-1) %/% n
    iy <- (ic-1) %% n
    D$X[cbind(ix+1,iy+1)] <- T
  })
  observeEvent(input$removeSet, {
    D$M[D$M==levels[2]] <- levels[3]
  })
  observeEvent(input$infectSet, {
    D$M[D$X & (D$M==levels[1])] <- levels[2]
    D$X[] <- F
    D$time <- D$time + 1
    if (sum(D$M==levels[2])==0) showModal(modalDialog(title="No new infections - simulation is finished"))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
