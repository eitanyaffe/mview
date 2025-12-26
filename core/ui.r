library(shiny)
library(DT)
library(ggplot2)
library(shinyjqui)
library(plotly)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(), # Enable shinyjs functionality
  div(style = "height: 10px;"), # Small gap at top
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
    ")),
    
    # ResizeObserver for dynamic container height
    tags$script(HTML("
      (function() {
        function debounce(fn, ms){
          let t; 
          return function(arg){
            clearTimeout(t);
            t = setTimeout(function(){ fn(arg); }, ms);
          };
        }
        function sendHeight(h){
          if (window.Shiny && Shiny.setInputValue) {
            Shiny.setInputValue('container_height', Math.round(h), {priority:'event'});
          }
        }
        const sendHeightDebounced = debounce(sendHeight, 200); // debounced updates
        
        function initResizeObserver(){
          const el = document.getElementById('resizable_profile_container');
          if(!el) { 
            setTimeout(initResizeObserver, 50); 
            return; 
          }
          
          // Send initial height after layout
          requestAnimationFrame(function(){
            sendHeight(el.getBoundingClientRect().height);
          });
          
          const ro = new ResizeObserver(entries => {
            for (const entry of entries) {
              const h = (entry.contentBoxSize && entry.contentBoxSize[0])
                        ? entry.contentBoxSize[0].blockSize
                        : entry.contentRect.height;
              // Only fire when user stops dragging (debounced)
              sendHeightDebounced(h);
            }
          });
          ro.observe(el);
          // Keep a reference in case hot-reload occurs
          window.__profile_container_ro__ = ro;
        }
        
        if (document.readyState === 'loading') {
          document.addEventListener('DOMContentLoaded', initResizeObserver);
        } else {
          initResizeObserver();
        }
      })();
    "))
  ),
  tags$style(HTML("
    .ui-resizable-s {
      height: 10px !important;
      background: #e9ecef;
      border-top: 1px solid #d0d7de;
      cursor: ns-resize;
    }
    /* Make the resizable container fill properly */
    #resizable_profile_container {
      border: 1px solid #ddd;
      border-radius: 6px;
      display: flex;
      width: 100%;
    }
    /* Make the inner two-pane fill the resizable container */
    #profile_two_pane {
      height: 100%;
      display: flex;
      gap: 12px;
      width: 100%;
    }
    /* Left plot column flexes; right sidebar fixed width */
    #profile-plots-column {
      flex: 1 1 auto;
      min-width: 400px;
      height: 100%;
      display: flex;
      flex-direction: column;
    }
    #parameter-panel-column {
      flex: 0 0 280px;
      max-width: 360px;
      min-width: 220px;
      height: 100%;
      overflow-y: auto;
      border-left: 1px solid #eee;
      padding-left: 12px;
    }
    /* Make plotly output stretch vertically */
    #plot_wrap { height: 100%; padding: 0; }
    #combined_plot { height: 100% !important; }
    /* Minimum height for resizable profile container */
    #resizable_profile_container { min-height: 600px !important; }
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
      flex: 0 0 50px !important;
      max-width: 50px;
      min-width: 50px;
    }
    #profile-plots-column.expanded {
      flex: 1 1 auto !important;
    }
    .navigation-panel {
      background-color: #f5f5f5;
      border: 1px solid #ddd;
      border-radius: 4px;
      padding: 5px 10px;
      margin-bottom: 10px;
      display: flex;
      align-items: center;
      gap: 15px;
    }
    .nav-button-group {
      display: flex;
      gap: 3px;
    }
    .nav-button-group .btn-sm {
      width: 30px;
      height: 30px;
      padding: 0;
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .nav-button-group .context-zoom-btn {
      width: 40px;
    }
  ")),
  fluidRow(
    column(
      width = 2,
      uiOutput("state_info"),
      h4("mview 1.0", style = "margin-bottom: 5px; margin-top: 10px; color: #333;"),
      verbatimTextOutput("project_info"),
      shiny::selectInput("regions_module-assembly_select", "Assembly:",
        choices = get_assemblies(),
        selected = {
          assemblies <- get_assemblies()
          default_assembly <- if (length(assemblies) > 0) assemblies[1] else ""
          cached_assembly <- cache_get_if_exists("assembly.selected", default_assembly)
          if (cached_assembly %in% assemblies) cached_assembly else default_assembly
        },
        multiple = FALSE,
        width = "100%"
      ),
      uiOutput("viewSelect"),
      h5("Info"),
      verbatimTextOutput("basic_info"),
      actionButton("plotViewBtn", "Export", icon = icon("file-pdf")),
      actionButton("plotRegionsBtn", "Export All", icon = icon("file-pdf")),
      actionButton("helpBtn", "Help"),
      actionButton("aboutBtn", "About"),
      br(), br(),
      uiOutput("state_buttons_ui")
    ),
    column(
      width = 10,


      # Profile plots area with resizable container and parameter panel
      uiOutput("resizableContainer"),
      hr(),
      uiOutput("mainTabsPanel")
    )
  )
)
