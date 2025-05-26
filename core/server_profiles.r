# ---- Dynamic Profile Rendering ----

library(plotly)

# UI side: dynamically generate plot outputs for all registered profiles
output$profilePlots <- renderUI({
  # Make this reactive to view changes
  req(input$view_id)
  req(state$contigs)
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }
  plotlyOutput(
    outputId = "combined_plot",
    height = "800px" # Adjust height as needed, or make dynamic
  )
})

# Server side: render plots for all registered profiles and cache info_df for hover
observe({
  cxt <- build_context(
    state_contigs = state$contigs,
    contig_table = get_contigs(state$assembly),
    zoom = state$zoom,
    assembly = state$assembly
  )
  req(input$view_id)
  req(cxt)

  state$plotly_registered <- FALSE

  # Get combined plotly object
  combined_plotly_obj <- plot_profiles(cxt)

  # Register the relayout event for the plot so zoom events can be captured with the correct source ID
  if (!is.null(combined_plotly_obj)) {
    combined_plotly_obj <- plotly::event_register(combined_plotly_obj, "plotly_relayout")
    # Register click for alignment details (source is implicitly the plot's outputId "combined_plot")
    # combined_plotly_obj <- plotly::event_register(combined_plotly_obj, "plotly_selected")

    # add a reactive flag to indicate that the plotly object has been registered
    state$plotly_registered <- TRUE
  }

  # Render the combined plot
  output$combined_plot <- renderPlotly({
    cat(sprintf("rendering combined plot\n"))
    combined_plotly_obj
  })
})

# Handle plotly zoom/relayout events to update state$zoom
observeEvent(eventExpr = plotly::event_data("plotly_relayout"), {
  req(state$plotly_registered)
  event_data <- suppressWarnings(plotly::event_data("plotly_relayout"))
  req(event_data)
  if (!is.null(event_data[["xaxis.range[0]"]]) && !is.null(event_data[["xaxis.range[1]"]])) {
    # Zoom event with defined range
    new_xlim <- c(event_data[["xaxis.range[0]"]], event_data[["xaxis.range[1]"]])
    # Check if xlim has actually changed
    if (!isTRUE(all.equal(state$current_xlim, new_xlim))) {
      state$current_xlim <- new_xlim
    }
  }
})
