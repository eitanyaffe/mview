library(shiny)
library(DT)
library(ggplot2)
library(shinyjqui)
library(plotly)
library(shinyjs)

ui <- fluidPage(
  titlePanel("mview"),
  useShinyjs(), # Enable shinyjs functionality
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
    .parameter-panel {
      border: 1px solid #ddd;
      border-radius: 4px;
      background-color: #f9f9f9;
      margin-left: 10px;
      transition: all 0.3s ease;
      overflow: hidden;
    }
    .parameter-panel-header {
      background-color: #e9e9e9;
      border-bottom: 1px solid #ddd;
      padding: 5px 10px;
      cursor: pointer;
    }
    .parameter-panel-header:hover {
      background-color: #e0e0e0;
    }
    .parameter-toggle-btn {
      transition: transform 0.3s ease;
    }
    .parameter-panel.collapsed .parameter-toggle-btn {
      transform: rotate(180deg);
    }
    .parameter-panel.collapsed #parameter-panel-content {
      display: none;
    }
    #profilePlots {
      min-height: 3in;
    }
    .parameter-panel-content {
      padding: 10px;
      max-height: 6in;
      overflow-y: auto;
      overflow-x: hidden;
    }
    #parameterTabs .tab-pane {
      max-height: 5.5in;
      overflow-y: auto;
      overflow-x: hidden;
    }
    .parameter-group-content {
      max-height: 5.5in;
      overflow-y: auto;
      overflow-x: hidden;
      padding-left: 10px;
    }
    .parameter-panel.collapsed {
      width: 50px !important;
      min-width: 50px;
    }
    .parameter-panel.collapsed .parameter-panel-header h5 {
      display: none !important;
      visibility: hidden !important;
      width: 0 !important;
      margin: 0 !important;
      padding: 0 !important;
    }
    #parameter-panel-column.collapsed {
      width: 50px !important;
      max-width: 50px;
      flex: 0 0 50px;
    }
    #profile-plots-column.expanded {
      width: calc(100% - 60px) !important;
      flex: 1 1 auto;
    }
    .expanded #profilePlots {
      width: 100% !important;
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
      uiOutput("viewSelect"),
      h4("Basic Info"),
      verbatimTextOutput("basic_info"),
      h4("Last Key Press"),
      verbatimTextOutput("last_key_output"),
      actionButton("helpBtn", "Help")
    ),
    column(
      width = 10,

      # Profile plots area with collapsible parameter panel
      fluidRow(
        column(
          width = 8,
          id = "profile-plots-column",
          uiOutput("profilePlots")
        ),
        column(
          width = 4,
          id = "parameter-panel-column",
          div(
            id = "parameter-panel",
            class = "parameter-panel",
            div(
              class = "parameter-panel-header",
              actionButton("toggleParameterPanel", "", 
                icon = icon("chevron-left"),
                class = "btn-sm parameter-toggle-btn",
                style = "float: right; margin: 5px;"
              ),
              h5("Parameters", style = "margin: 5px 10px; display: inline-block;")
            ),
            div(
              id = "parameter-panel-content",
              uiOutput("parameter_tabs_ui")
            )
          )
        )
      ),
      hr(),
      uiOutput("mainTabsPanel")
    )
  )
)
