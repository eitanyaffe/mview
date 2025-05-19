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
    height = "500px" # Adjust height as needed, or make dynamic
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
    # add a reactive flag to indicate that the plotly object has been registered
    state$plotly_registered <- TRUE
  }

  # Render the combined plot
  output$combined_plot <- renderPlotly({
    combined_plotly_obj
  })
})

# Handle plotly zoom/relayout events to update state$zoom
observeEvent(eventExpr = plotly::event_data("plotly_relayout"), {
  return()

  req(state$plotly_registered)
  event_data <- suppressWarnings(plotly::event_data("plotly_relayout"))
  req(event_data)
  cat(sprintf("relayout event_data: %s\n", paste(names(event_data), collapse = ", ")))
  print(event_data)
  if (!is.null(event_data[["xaxis.range[0]"]]) && !is.null(event_data[["xaxis.range[1]"]])) {
    # Zoom event with defined range
    new_zoom <- c(event_data[["xaxis.range[0]"]], event_data[["xaxis.range[1]"]])
    # Check if zoom has actually changed to avoid infinite loops
    if (!isTRUE(all.equal(state$zoom, new_zoom))) {
      state$zoom <- new_zoom
      cat(sprintf("zoom set by plotly: %.1f - %.1f\n", new_zoom[1], new_zoom[2]))
    }
  }
})
