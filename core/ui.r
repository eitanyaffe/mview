library(shiny)
library(DT)
library(ggplot2)
library(shinyjqui)
library(plotly)

ui <- fluidPage(
  titlePanel("mview"),
  tags$head(
    keyboard_initialize() # Add the keyboard event listener
  ),
  tags$style(HTML("
    .ui-resizable-s {
      height: 8px;
      background: #cccccc;
      cursor: ns-resize;
    }
    .state-info-box {
      border: 1px solid #ddd;
      padding: 8px;
      background-color: #f8f8f8;
      border-radius: 4px;
      margin-bottom: 10px;
      font-family: monospace;
    }
  ")),
  fluidRow(
    column(
      width = 2,
      h4("State"),
      uiOutput("state_info"),
      shiny::selectInput("states_module-assembly_select", "Select Assembly:",
        choices = get_assemblies(),
        selected = get_assemblies()[1],
        multiple = FALSE,
        width = "100%"
      ),
      actionButton("clearLog", "Clear Log"),
      verbatimTextOutput("log"),
      h4("Last Key Press"),
      verbatimTextOutput("last_key_output"),
      actionButton("helpBtn", "Help")
    ),
    column(
      width = 10,
      fluidRow(
        column(width = 6, uiOutput("viewSelect")),
        column(width = 6, div(class = "state-info-box", verbatimTextOutput("current_state_display")))
      ),
      uiOutput("profilePlots"),
      hr(),
      tabsetPanel(
        id = "mainTabs",
        tabPanel(
          "States",
          states_ui("states_module")
        ),
        tabPanel(
          "Contigs",
          actionButton("addContigsBtn", "Add selected contigs"),
          DTOutput("contigTable")
        ),
        tabPanel(
          "Genomes",
          actionButton("addGenomesBtn", "Add contigs of selected genomes"),
          DTOutput("genomeTable")
        ),
        tabPanel(
          "Contig Map",
          DTOutput("mapTable")
        ),
        tabPanel(
          "Selected Contigs",
          actionButton("removeContigsBtn", "Remove selected contigs"),
          DTOutput("selectedTable")
        ),
        tabPanel(
          "Genes",
          div(
            style = "margin-bottom: 10px;",
            actionButton("showGeneDetailsBtn", "Show Details"),
            actionButton("zoomToGeneBtn", "Zoom to Gene")
          ),
          DTOutput("genesTable")
        ),
        tabPanel(
          "Options",
          numericInput("log.length", "Max log messages:",
            value = 10,
            min = 1,
            max = 1000
          )
        ),
        tabPanel(
          "Parameters",
          uiOutput("parameters_ui")
        )
      )
    )
  )
)
