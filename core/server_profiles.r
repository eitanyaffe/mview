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
  
  # Calculate total height from profile heights (including dynamic parameter values)
  total_height <- sum(sapply(profiles, function(p) {
    height <- p$height  # start with static height
    
    # update height from parameters if it exists
    if ("height" %in% names(p$params)) {
      height_param <- p$params[["height"]]
      if (!is.null(height_param) && !is.null(height_param$group)) {
        height <- get_param(height_param$group, "height")()
      }
    }
    
    return(height)
  }))
  if (total_height <= 0) {
    total_height <- 200  # Default fallback
  }
  
  # Apply min/max constraints from options
  min_height <- if (!is.null(input$min_plot_height)) input$min_plot_height else 200
  max_height <- if (!is.null(input$max_plot_height)) input$max_plot_height else 1200
  total_height <- max(min_height, min(max_height, total_height))
  cat(sprintf("total_height: %dpx\n", total_height))
  plotlyOutput(
    outputId = "combined_plot",
    height = paste0(total_height, "px")
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

  # Get combined plotly object and height info
  plot_result <- plot_profiles(cxt)
  combined_plotly_obj <- if (!is.null(plot_result)) plot_result$plot else NULL
  
  # Store legends in state for the legend tab
  state$current_legends <- if (!is.null(plot_result)) plot_result$legends else list()

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

# Convert pixel coordinates to biological coordinates
convert_mouse_coords <- function(mouse_data, state) {
  if (is.null(mouse_data$plotX) || is.null(mouse_data$plotWidth)) {
    return(list(x = NULL, y = NULL))
  }
  
  # get current x-axis range
  xlim <- get_current_xlim(state)
  
  # convert pixel position to biological coordinates
  rel_x <- mouse_data$plotX / mouse_data$plotWidth
  plot_x <- xlim[1] + rel_x * (xlim[2] - xlim[1])
  
  list(x = plot_x, y = mouse_data$plotY)
}

# Get current x-axis limits from state
get_current_xlim <- function(state) {
  if (!is.null(state$current_xlim)) {
    return(state$current_xlim)
  }
  if (!is.null(state$zoom)) {
    return(state$zoom)
  }
  
  # build context to get full range
  cxt <- build_context(
    state_contigs = state$contigs,
    contig_table = get_contigs(state$assembly),
    zoom = NULL,
    assembly = state$assembly
  )
  
  if (!is.null(cxt) && !is.null(cxt$mapper)) {
    return(cxt$mapper$xlim)
  }
  
  c(0, 1000)  # fallback
}

# Handle JavaScript mouse coordinates from UI
observeEvent(input$mouse_coords, {
  mouse_data <- input$mouse_coords
  req(mouse_data)
  
  plot_coords <- convert_mouse_coords(mouse_data, state)
  
  state$mouse_coords <- list(
    screen_x = mouse_data$x,
    screen_y = mouse_data$y,
    plot_x = plot_coords$x,
    plot_y = plot_coords$y,
    timestamp = Sys.time()
  )
})
