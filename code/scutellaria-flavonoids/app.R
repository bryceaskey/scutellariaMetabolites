# Interactive visualization of flavonoid data collected for various species of Scutellaria
# Before running, working directory must be set to parent directory of app.R

# Load necessary libraries, functions, and data ----
packageList <- c("shiny", "shinyjs", "plyr", "dplyr", "tibble", "ggplot2", "ggrepel", "cowplot")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0){
  install.packages(newPackages)
}

library(shiny)
library(shinyjs)
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
  useShinyjs(),
  titlePanel("Scutellaria flavonoids"),
  fluidRow(
    column(3,
      wellPanel(
        helpText("Interactive visualization of medicinally relevant flavonoid concentrations for
                 several species of Scutellaria"),
        selectInput("organ", "Select organ",
                    choices=list("Leaves", "Shoots", "Roots"), selected=1),
        checkboxGroupInput("metabolites", "Select metabolites",
                           choices=list("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin",
                                        "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside",
                                        "scutellarein", "scutellarin", "baicalin", "baicalein",
                                        "wogonin", "wogonoside")),
        actionButton("resetButton", "Reset selections")
      )
    ),
    
    column(9,
      fluidRow(
        column(10, plotOutput("plot", click="plotClick")),
        column(2, plotOutput("legend"))
      ),
      br(),
      h4(textOutput("clickInfo")),
      tableOutput("selectedData")
    )
  )
)


# Define server logic ----
server <- function(input, output){
  # Return reactive plot
  output$plot <- renderPlot({
    if(length(input$metabolites) > 0){
      filteredData <- filterData(allData, input$organ, input$metabolites)
      return(createPlot(filteredData, metaboliteColors))
    }else{
      # Should return an empty plot with a gray background - y-axis scaled from 0 to 1
      return(createEmptyPlot(allData, input$organ))
    }
  })
    
  # Return legend
  output$legend <- renderPlot({
    return(createLegend(allData, metaboliteColors))
  })
  
  # Functionality for reset button - uncheck metabolite boxes
  observeEvent(input$resetButton, {
    reset("metabolites")
  })
  
  # Return species selected by click on plot
  output$clickInfo <- renderText({
    if(is.null(input$plotClick$x)){
      paste0("Selected species: ")
    }else if(abs(input$plotClick$x - round(input$plotClick$x)) > 0.275){
      paste0("Selected species: ")
    }else{
      filteredData <- filterData(allData, input$organ, input$metabolites)
      paste0("Selected species: ", levels(filteredData$species)[round(input$plotClick$x)])
    }
  })
  
  # Return data about selected species in table
  output$selectedData <- renderTable({expr=
    if(is.null(input$plotClick$x)){
      emptyDF <- data.frame(a = character(), b = character(), c = character(), d = character())
      colnames(emptyDF) <- c("Flavonoid #", "Flavonoid name", "Concentration (ppm)", "Standard error")
      return(emptyDF)
    }else if(abs(input$plotClick$x - round(input$plotClick$x)) > 0.275){
      emptyDF <- data.frame(a = character(), b = character(), c = character(), d = character())
      colnames(emptyDF) <- c("Flavonoid #", "Flavonoid name", "Concentration (ppm)", "Standard error")
      return(emptyDF)
    }else{
      selectedSpecies <- levels(allData$species)[round(input$plotClick$x)]
      return(refineData(allData, input$organ, input$metabolites, selectedSpecies))
    }
  width=600
  })
}


# Run the app
shinyApp(ui=ui, server=server)