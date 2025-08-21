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
      "Legends",
      uiOutput("legendsContent")
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
      p("These settings control the minimum and maximum allowed height for the combined profile plots."),
      br(),
      h4("Legend Settings"),
      numericInput(
        inputId = "legend_scale",
        label = "Legend scale factor",
        value = 1.0,
        min = 0.5,
        max = 3.0,
        step = 0.1
      ),
      p("Scale factor for all legend sizes. Use values < 1 to make legends smaller, > 1 to make them larger."),
      br(),
      h4("Contig Length Filtering"),
      numericInput(
        inputId = "min_contig_length",
        label = "Minimum contig length (kb)",
        value = cache_get_if_exists("input_min_contig_length", 0),
        min = 0,
        step = 1
      ),
      numericInput(
        inputId = "max_contig_length",
        label = "Maximum contig length (kb)",
        value = cache_get_if_exists("input_max_contig_length", 100000),
        min = 0,
        step = 1000
      ),
      p("Filter contigs displayed in the contig table by length. Adjust values to focus on contigs within specific size ranges.")
    )
  )
  
  # get dynamic tabs from registration
  dynamic_tab_panels <- get_tab_panels()
  
  # combine fixed and dynamic tabs
  all_tabs <- c(fixed_tabs, dynamic_tab_panels)
  
  # create the tabsetPanel with all tabs
  do.call(tabsetPanel, c(list(id = "mainTabs"), all_tabs))
})

# observers to cache contig length filter values
observeEvent(input$min_contig_length, {
  cache_set("input_min_contig_length", input$min_contig_length)
})

observeEvent(input$max_contig_length, {
  cache_set("input_max_contig_length", input$max_contig_length)
})


