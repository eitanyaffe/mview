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
    ")),
    fluidRow(
      # left column: segment organizer
            column(3,
        h4("Segments in view"),
        p("Click to select, Alt/Cmd+click for multi-select, drag to reorder", style = "color: #666; font-size: 11px; margin-bottom: 8px;"),
        uiOutput("selectorSegmentListUI"),
        uiOutput("selectorOrganizerButtons")
      ),
      # right column: graph with controls
      column(9,
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
              selectInput("selectorColorBy", "Color by:",
                         choices = c("bin" = "bin"),
                         selected = "bin")
            ),
            column(3,
              numericInput("selectorNeighborDepth", "Neighbor depth:", 
                         value = 0, min = 0, max = 3, step = 1)
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
                actionButton("selectorClearSelectionBtn", "Clear Selection", class = "btn btn-sm btn-default")
              )
            )
          ),
          div(style = "flex: 1;",
            wellPanel(style = "padding: 10px; height: 100%; margin: 0;",
              uiOutput("selectorHoverInfo")
            )
          )
        ),
        # graph zoom controls
        div(style = "display: flex; align-items: center; margin-top: 10px; margin-bottom: 5px; gap: 5px;",
          actionButton("selectorGraphReset", "Reset", class = "btn btn-sm btn-default"),
          actionButton("selectorGraphZoomIn", "+", class = "btn btn-sm btn-default"),
          actionButton("selectorGraphZoomOut", "−", class = "btn btn-sm btn-default")
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
