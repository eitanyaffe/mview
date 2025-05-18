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
  ")),
  fluidRow(
    column(
      width = 2,
      h4("State"),
      uiOutput("state_info"),
      actionButton("clearLog", "Clear Log"),
      verbatimTextOutput("log"),
      h4("Last Key Press"), # Label for the new output
      verbatimTextOutput("last_key_output"), # New output for last key press
      actionButton("helpBtn", "Help")
    ),
    column(
      width = 10,
      uiOutput("viewSelect"),
      uiOutput("profilePlots"), # Plot area with zoom brush
      hr(),
      tabsetPanel(
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
        ),
        tabPanel(
          "States",
          states_ui("states_module")
        )
      )
    )
  )
)
