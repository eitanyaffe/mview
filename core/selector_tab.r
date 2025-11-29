# selector tab implementation
# Core tab for selecting and ordering segments via interactive graph

# source component files
source("core/selector_tab_state.r", local = TRUE)
source("core/selector_tab_utils.r", local = TRUE)
source("core/selector_tab_server.r", local = TRUE)

# UI function (will be called from server_tabs.r)
selector_tab_ui <- function() {
  tabPanel(
    "Selector",
    fluidRow(
      column(12,
        wellPanel(
          h5("Controls"),
          fluidRow(
            column(3,
              textInput("selectorBinInput", "Bin:", 
                       value = "",
                       placeholder = "Enter bin ID (e.g., b1)")
            ),
            column(3,
              actionButton("selectorLoadBinBtn", "Load Bin", 
                         class = "btn btn-primary")
            ),
            column(3,
              actionButton("selectorAddNeighborsBtn", "Add Neighbors", 
                         class = "btn btn-default")
            ),
            column(3,
              actionButton("selectorResetBtn", "Reset", 
                         class = "btn btn-default"),
              actionButton("selectorUpdateBtn", "Update", 
                         class = "btn btn-primary")
            )
          ),
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
            column(6,
              h6("Sequence Actions:"),
              actionButton("selectorAddAfterBtn", "Add After", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorAddBeforeBtn", "Add Before", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorAddFirstBtn", "Add First", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorAddLastBtn", "Add Last", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorRemoveBtn", "Remove", 
                         class = "btn btn-sm btn-default")
            )
          )
        )
      )
    ),
    fluidRow(
      column(8,
        h4("Graph (click nodes/edges to select)"),
        visNetwork::visNetworkOutput("selectorGraph", height = "500px")
      ),
      column(4,
        h4("Sequence (drag to reorder, click ✕ to remove)"),
        uiOutput("selectorSequenceUI")
      )
    ),
    fluidRow(
      column(12,
        h4("Edge Data"),
        DTOutput("selectorEdgeTable")
      )
    )
  )
}
