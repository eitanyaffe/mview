# selector tab
# minimal entry point that sources components and defines UI

# source component files
source("core/organizer_state.r", local = TRUE)
source("core/organizer_server.r", local = TRUE)
source("core/graph_state.r", local = TRUE)
source("core/graph_utils.r", local = TRUE)
source("core/graph_server.r", local = TRUE)

# Organizer tab UI (merged with graph)
organizer_tab_ui <- function() {
  tabPanel(
    "Organizer (O)",
    # load sortable.js scripts
    singleton(tags$script(src = "https://unpkg.com/sortablejs@1.15.0/Sortable.min.js")),
    tags$style(HTML("
      #selectorSegmentList { padding: 0; margin: 0; max-height: 800px; overflow-y: auto; width: 100%; }
      #selectorSegmentList li {
        list-style: none;
        margin: 3px 0;
        padding: 6px 10px;
        border-radius: 4px;
        cursor: grab;
        user-select: none;
      }
      #selectorSegmentList li.selected {
        outline: 2px solid #333;
        opacity: 0.75;
      }
      .selector-buttons {
        display: flex;
        gap: 5px;
        margin-top: 8px;
      }
      .selector-buttons .btn {
        padding: 4px 10px;
        font-size: 12px;
      }
      .selector-buttons .btn.faded {
        opacity: 0.4;
        pointer-events: none;
      }
      .organizer-container {
        display: flex;
        width: 100%;
        gap: 0;
      }
      .organizer-left {
        width: 350px;
        min-width: 150px;
        max-width: 500px;
        padding-right: 10px;
      }
      .organizer-splitter {
        width: 8px;
        background: #e9ecef;
        cursor: ew-resize;
        flex-shrink: 0;
        border-left: 1px solid #d0d7de;
        border-right: 1px solid #d0d7de;
      }
      .organizer-splitter:hover {
        background: #dee2e6;
      }
      .organizer-right {
        flex: 1;
        padding-left: 10px;
        min-width: 400px;
      }
    ")),
    tags$script(HTML("
      $(document).ready(function() {
        var splitter = document.querySelector('.organizer-splitter');
        var left = document.querySelector('.organizer-left');
        if (!splitter || !left) return;
        
        var dragging = false;
        splitter.addEventListener('mousedown', function(e) {
          dragging = true;
          document.body.style.cursor = 'ew-resize';
          document.body.style.userSelect = 'none';
        });
        document.addEventListener('mousemove', function(e) {
          if (!dragging) return;
          var container = document.querySelector('.organizer-container');
          var rect = container.getBoundingClientRect();
          var newWidth = e.clientX - rect.left;
          newWidth = Math.max(150, Math.min(500, newWidth));
          left.style.width = newWidth + 'px';
        });
        document.addEventListener('mouseup', function() {
          dragging = false;
          document.body.style.cursor = '';
          document.body.style.userSelect = '';
        });
      });
    ")),
    div(class = "organizer-container",
      # left column: segment organizer
      div(class = "organizer-left",
        h4("Segments in view"),
        selectInput("segmentColorScheme", "Color by:",
                   choices = c("order"),
                   selected = "order",
                   width = "100%"),
        checkboxInput("segmentGrayscale", "Grayscale", value = FALSE),
        p("Click to select, Alt/Cmd+click for multi-select, drag to reorder", style = "color: #666; font-size: 11px; margin-bottom: 8px;"),
        uiOutput("selectorSegmentListUI"),
        uiOutput("selectorActionButtons"),
        uiOutput("selectorApplyResetButtons")
      ),
      # splitter
      div(class = "organizer-splitter"),
      # right column: graph with controls
      div(class = "organizer-right",
        # graph parameters box
        wellPanel(style = "padding: 10px;",
          fluidRow(
            column(3,
              numericInput("selectorMinSupport", "Min Support:", 
                         value = cache_get_if_exists("selector.min_support", 0),
                         min = 0, step = 1)
            ),
            column(3,
              numericInput("selectorMinPercent", "Min Percent:", 
                         value = cache_get_if_exists("selector.min_percent", 0),
                         min = 0, max = 100, step = 0.1)
            ),
            column(3,
              numericInput("selectorNeighborDepth", "Neighbor depth:", 
                         value = 1, min = 0, max = 3, step = 1)
            )
          )
        ),
        # edge controls
        wellPanel(style = "padding: 10px; margin-top: 10px;",
          fluidRow(
            column(3,
              selectInput("graphEdgeMetric", "Edge metric:",
                         choices = c("none" = "none", "support" = "support", 
                                   "percent" = "percent", "change" = "change"),
                         selected = cache_get_if_exists("graph.edge_metric", "none"),
                         width = "100%")
            ),
            column(3,
              conditionalPanel(
                condition = "input.graphEdgeMetric == 'change'",
                selectInput("graphEdgeLib1", "Library 1:",
                           choices = list(),
                           selected = NULL,
                           width = "100%")
              )
            ),
            column(3,
              conditionalPanel(
                condition = "input.graphEdgeMetric == 'change'",
                selectInput("graphEdgeLib2", "Library 2:",
                           choices = list(),
                           selected = NULL,
                           width = "100%")
              )
            ),
            column(3,
              checkboxInput("graphEdgeLabels", "Edge labels", 
                          value = cache_get_if_exists("graph.edge_labels", FALSE))
            )
          )
        ),
        # controls split: left half (update graph + selection), right half (hover info)
        div(style = "display: flex; margin-top: 10px; gap: 15px;",
          div(style = "flex: 2;",
            wellPanel(style = "padding: 10px; height: 100%; margin: 0;",
              tags$strong("Update graph"),
              div(style = "margin-top: 8px; display: flex; align-items: center; gap: 10px;",
                actionButton("selectorUpdateToView", "Current view", class = "btn btn-sm btn-primary"),
                actionButton("selectorUpdateToZoom", "Current zoom", class = "btn btn-sm btn-default"),
                radioButtons("selectorAutoUpdate", "Automatic:",
                            choices = c("Off" = "off", "Track View" = "view", "Track Zoom" = "zoom"),
                            selected = "off", inline = TRUE)
              ),
              div(style = "margin-top: 10px;",
                tags$strong("Selected: "),
                textOutput("selectorSelectedSegmentText", inline = TRUE)
              ),
              div(style = "margin-top: 10px;",
                actionButton("selectorGotoBtn", "Goto", class = "btn btn-sm btn-primary"),
                actionButton("selectorAddBtn", "Add", class = "btn btn-sm btn-default"),
                actionButton("selectorRemoveBtn", "Remove", class = "btn btn-sm btn-default"),
                actionButton("selectorSelectAllBtn", "Select All", class = "btn btn-sm btn-default"),
                actionButton("selectorClearSelectionBtn", "Clear", class = "btn btn-sm btn-default")
              )
            )
          ),
          div(style = "flex: 1.1;",
            wellPanel(style = "padding: 10px; height: 100%; margin: 0;",
              uiOutput("selectorHoverInfo")
            )
          )
        ),
        # graph controls
        div(style = "display: flex; align-items: center; margin-top: 10px; margin-bottom: 5px; gap: 10px;",
          actionButton("selectorGraphReset", "Fit", class = "btn btn-sm btn-default"),
          actionButton("selectorGraphZoomIn", "+", class = "btn btn-sm btn-default"),
          actionButton("selectorGraphZoomOut", "−", class = "btn btn-sm btn-default"),
          tags$span("Label:", style = "margin-left: 10px;"),
          selectInput("graphNodeLabel", "", choices = c("index" = "index", "id" = "id"), selected = "index", width = "80px"),
          tags$span("Font:", style = "margin-left: 10px;"),
          selectInput("graphFontSize", "", choices = c("S" = "28", "M" = "38", "L" = "48", "XL" = "56"), selected = "38", width = "70px"),
          checkboxInput("graphDirectedEdges", "Directed", value = FALSE)
        ),
        # graph at bottom
        div(
          style = "border: 1px solid #ccc; border-radius: 4px; padding: 5px;",
          visNetwork::visNetworkOutput("selectorGraph", height = "600px")
        )
      )
    )
  )
}
