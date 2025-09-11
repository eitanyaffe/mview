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
    keyboard_initialize(), # Add the keyboard event listener
    tags$script(HTML("
      $(document).ready(function() {
        var plotContainer = null;
        var hoverTimeout = null;
        var HOVER_DELAY = 30;
        
        function findPlotContainer() {
          return $('#combined_plot').find('.plotly').first();
        }
        
        function getPlotCoords(mouseX, mouseY) {
          if (!plotContainer || !plotContainer.length) return null;
          
          var plotOffset = plotContainer.offset();
          var plotWidth = plotContainer.width();
          var plotHeight = plotContainer.height();
          var relX = mouseX - plotOffset.left;
          var relY = mouseY - plotOffset.top;
          
          if (relX < 0 || relX > plotWidth || relY < 0 || relY > plotHeight) {
            return null;
          }
          
          return {
            x: mouseX,
            y: mouseY,
            plotX: relX,
            plotY: relY,
            plotWidth: plotWidth,
            plotHeight: plotHeight
          };
        }
        
        function updateMouseCoords(e) {
          if (!plotContainer) plotContainer = findPlotContainer();
          
          var coords = getPlotCoords(e.pageX, e.pageY);
          Shiny.setInputValue('mouse_coords', coords, {priority: 'event'});
        }
        
        $(document).on('mousemove', function(e) {
          if (hoverTimeout) clearTimeout(hoverTimeout);
          hoverTimeout = setTimeout(function() { updateMouseCoords(e); }, HOVER_DELAY);
        });
        
        $(document).on('mouseleave', function() {
          if (hoverTimeout) clearTimeout(hoverTimeout);
        });
        
        $(document).on('shiny:value', function(event) {
          if (event.target.id === 'combined_plot') {
            plotContainer = null;
            if (hoverTimeout) clearTimeout(hoverTimeout);
            setTimeout(function() { plotContainer = findPlotContainer(); }, 100);
          }
        });
      });
      
      // handle focus input message from regions module
      Shiny.addCustomMessageHandler('focusInput', function(inputId) {
        document.getElementById(inputId).focus();
      });
      
      // handle updating readonly inputs
      Shiny.addCustomMessageHandler('updateReadonlyInput', function(data) {
        document.getElementById(data.id).value = data.value;
      });
      
      // handle triggering save dialog from keyboard shortcuts
      Shiny.addCustomMessageHandler('triggerStateSave', function(stateNumber) {
        document.getElementById('save_state_' + stateNumber).click();
      });
    "))
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
      max-height: 8in;
      overflow-y: auto;
      overflow-x: hidden;
    }
    #parameterTabs .tab-pane {
      max-height: 7.5in;
      overflow-y: auto;
      overflow-x: hidden;
    }
    .parameter-group-content {
      max-height: 7.5in;
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
      flex: 0 0 50px !important;
    }
    #profile-plots-column.expanded {
      width: calc(100% - 60px) !important;
      flex: 1 1 auto !important;
    }
    .expanded #profilePlots {
      width: 100% !important;
    }
  ")),
  fluidRow(
    column(
      width = 2,
      uiOutput("state_info"),
      verbatimTextOutput("project_info"),
      shiny::selectInput("regions_module-assembly_select", "Assembly:",
        choices = get_assemblies(),
        selected = get_assemblies()[1],
        multiple = FALSE,
        width = "100%"
      ),
      uiOutput("viewSelect"),
      h5("Info"),
      verbatimTextOutput("basic_info"),
      h5("Last Key Press"),
      verbatimTextOutput("last_key_output"),
      uiOutput("refreshBtnUI"),
      actionButton("plotViewBtn", "Export", icon = icon("file-pdf")),
      actionButton("plotRegionsBtn", "Export All", icon = icon("file-pdf")),
      actionButton("helpBtn", "Help"),
      actionButton("aboutBtn", "About"),
      br(), br(),
      uiOutput("state_buttons_ui")
    ),
    column(
      width = 10,

      # Profile plots area with collapsible parameter panel
      div(
        style = "display: flex; width: 100%;",
        div(
          id = "profile-plots-column",
          style = "flex: 1 1 auto; width: 79.17%;",
          uiOutput("profilePlots")
        ),
        div(
          id = "parameter-panel-column",
          style = "flex: 0 0 20.83%; width: 20.83%;",
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
