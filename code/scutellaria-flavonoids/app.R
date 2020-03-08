# Interactive visualization of flavonoid data collected for various species of Scutellaria

# Load necessary libraries, functions, and data ----
library(shiny)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(cowplot)

source("helpers.R")

load("data/flavonoidConcs.rda") # name of variable is "allData"

# Define UI ----
ui <- fluidPage(
  titlePanel("Scutellaria flavonoids"),
  sidebarLayout(
    sidebarPanel(
      helpText("Interactive visualization of medicinally relevant flavonoid concentrations for
               several species of Scutellaria"),
      selectInput("organ", "Select organ",
                  choices=list("Leaves", "Shoots", "Roots"), selected=1),
      checkboxGroupInput("metabolites", "Select metabolites",
                         choices=list("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin",
                                      "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside",
                                      "scutellarein", "scutellarin", "baicalin", "baicalein",
                                      "wogonin", "wogonoside"))
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)


# Define server logic ----
server <- function(input, output){
  output$plot <- renderPlot({
    if(length(input$metabolites) > 0){
      filteredData <- allData %>%
        filter(grepl(input$organ, organ)) %>%
        filter(grepl(paste(paste("^", input$metabolites, "$", sep=""), collapse="|"), metabolite))
      return(createPlot(filteredData, createLegend(allData, metaboliteColors)))
    }else{
      return("")
    }
  })
}


# Run the app
shinyApp(ui=ui, server=server)