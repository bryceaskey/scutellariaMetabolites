library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel("censusVis"),
  sidebarLayout(
    sidebarPanel(
      helpText("Create demographic maps with information form the 2010 US Census."),
      selectInput("var", h3("Choose a variable to display"),
                  choices=list("Percent White", "Percent Black", "Percent Hispanic", "Percent Asian"),
                  selected=1),
      sliderInput("range", h3("Range of interest"),
                  min=0, max=100, value=c(0, 100))
    ),
    mainPanel(
      textOutput("selected_var"),
      textOutput("selected_range")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  output$selected_var <- renderText({
    paste("You have selected", input$var)
  })
  output$selected_range <- renderText({
    paste("You have chosen a range that goes from", input$range[1], "to", input$range[2])
  })
}

# Run the app
shinyApp(ui=ui, server=server)