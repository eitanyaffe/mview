# core/server_tabs.r
# Server-side rendering for dynamic tabs

# render the main tabs panel with both fixed and dynamic tabs
output$mainTabsPanel <- renderUI({
  # fixed tabs that are always present
  fixed_tabs <- list(
    tabPanel(
      "Regions",
      regions_ui("regions_module")
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
      h4("Plot Height Settings"),
      numericInput(
        inputId = "min_plot_height",
        label = "Minimum total plot height (px)",
        value = 200,
        min = 100,
        max = 1000,
        step = 50
      ),
      numericInput(
        inputId = "max_plot_height", 
        label = "Maximum total plot height (px)",
        value = 1200,
        min = 500,
        max = 3000,
        step = 100
      ),
      p("These settings control the minimum and maximum allowed height for the combined profile plots.")
    )
  )
  
  # get dynamic tabs from registration
  dynamic_tab_panels <- get_tab_panels()
  
  # combine fixed and dynamic tabs
  all_tabs <- c(fixed_tabs, dynamic_tab_panels)
  
  # create the tabsetPanel with all tabs
  do.call(tabsetPanel, c(list(id = "mainTabs"), all_tabs))
})
