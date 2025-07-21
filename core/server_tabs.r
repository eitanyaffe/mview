# core/server_tabs.r
# Server-side rendering for dynamic tabs

# render the main tabs panel with both fixed and dynamic tabs
output$mainTabsPanel <- renderUI({
  # fixed tabs that are always present
  fixed_tabs <- list(
    tabPanel(
      "States",
      states_ui("states_module")
    ),
    tabPanel(
      "Contigs",
      actionButton("addContigsBtn", "Add selected contigs"),
      actionButton("setContigsBtn", "Set selected contigs"),
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
    )
  )
  
  # get dynamic tabs from registration
  dynamic_tab_panels <- get_tab_panels()
  
  # combine fixed and dynamic tabs
  all_tabs <- c(fixed_tabs, dynamic_tab_panels)
  
  # create the tabsetPanel with all tabs
  do.call(tabsetPanel, c(list(id = "mainTabs"), all_tabs))
})
