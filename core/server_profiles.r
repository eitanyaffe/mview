# ---- Dynamic Profile Rendering ----

# Store info_df for each profile for efficient hover
profile_info_dfs <- reactiveValues()

# UI side: dynamically generate plot outputs for all registered profiles
output$profilePlots <- renderUI({
  # Make this reactive to view changes
  req(input$view_id)
  req(state$contigs)
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }
  weights <- sapply(profiles, function(p) p$height)
  total <- sum(weights)
  brush_opts <- brushOpts(id = "zoomBrush", direction = "x")

  plot_list <- lapply(names(profiles), function(id) {
    prof <- profiles[[id]]
    plotOutput(
      outputId = paste0("plot_", id),
      height = paste0(300 * prof$height / total, "px"),
      hover = hoverOpts(id = paste0("hover_", id), delay = 50),
      brush = brush_opts,
      dblclick = "plot_dblclick"
    )
  })
  tagList(plot_list)
})

# Server side: render plots for all registered profiles and cache info_df for hover
observe({
  cxt <- build_context(
    state_contigs = state$contigs,
    contig_table = contigs,
    zoom = state$zoom
  )
  req(cxt)
  plots <- plot_profiles(cxt)
  for (id in names(plots)) {
    local({
      myid <- id
      myplot <- plots[[id]]$plot
      myinfo <- plots[[id]]$info_df
      output[[paste0("plot_", myid)]] <- renderPlot({
        myplot
      })
      profile_info_dfs[[myid]] <- myinfo
    })
  }
})
