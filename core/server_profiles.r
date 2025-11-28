# ---- Dynamic Profile Rendering ----

library(plotly)

# Reactive value to store cached plotly objects
profile_plotly_objects <- reactiveVal(NULL)

# UI side: generate resizable container with cached height
output$resizableContainer <- renderUI({
  height <- cache_get_if_exists("profile_panel_height", 700)  # Default 600px
  
  jqui_resizable(
    div(
      id = "resizable_profile_container",
      style = paste0("height: ", height, "px; padding: 8px;"),
      div(
        id = "profile_two_pane",
        # Left = plot area
        div(
          id = "profile-plots-column",
          uiOutput("navigationPanelUI"),
          div(
            id = "plot_wrap",
            style = "flex: 1 1 auto; height: 100%; padding: 0;",
            plotlyOutput("combined_plot", height = "100%")
          )
        ),
        # Right = parameters
        div(
          id = "parameter-panel-column",
          div(
            id = "parameter-panel",
            class = "parameter-panel",
            div(
              class = "parameter-panel-header",
              actionButton("toggleParameterPanel", "", 
                icon = icon("chevron-left"),
                class = "btn-sm parameter-toggle-btn",
                style = "float: right; margin: 5px;"
              ),
              h5("Parameters", style = "margin: 5px 10px; display: inline-block;")
            ),
            div(
              id = "parameter-panel-content",
              uiOutput("parameter_tabs_ui")
            )
          )
        )
      )
    ),
    options = list(handles = "s", minHeight = 200, maxHeight = 1200)
  )
})

# UI side: dynamically generate plot outputs for all registered profiles
output$profilePlots <- renderUI({
  # Make this reactive to view changes
  req(input$view_id)
  req(nrow(get_state_segments()) > 0)
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }
  
  # Calculate total height from profile heights (including dynamic parameter values)
  total_height <- sum(sapply(profiles, function(p) {
    # apply default height fields if missing
    if (is.null(p$is_fixed)) p$is_fixed <- FALSE
    if (is.null(p$height)) p$height <- 150
    
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
  
  # Apply max constraint (1500px max, no min constraint for initial height)
  total_height <- min(total_height, 1500)
  div(
    id = "plot_wrap",
    style = "flex: 1 1 auto; height: 100%; padding: 0;",
    plotlyOutput("combined_plot", height = "100%")
  )
})

# reactive trigger for manual refresh
refresh_trigger <- reactiveVal(0)

# reactive value to track plot sync status (TRUE = up to date, FALSE = needs refresh)
plot_updated <- reactiveVal(TRUE)

# global function to invalidate plots (called from tabs when changes occur)
invalidate_plot <- function() {
  plot_updated(FALSE)
}

# Observer to update cached ggplot objects when data/profiles change
observe({
  refresh_trigger()  # establish reactive dependency on manual refresh trigger
  
  # Update context state
  cxt_set_assembly(state$assembly)
  cxt_set_view(get_state_segments())
  
  if (!is.null(state$zoom)) {
    cxt_set_zoom(state$zoom)
  }
  
  req(input$view_id)

  # Update the cached plotly objects
  update_gg_objects(profile_plotly_objects)
  
  # mark plots as up to date after successful plotting
  plot_updated(TRUE)
})

# Observer to render plots when cached objects OR container height changes
observe({
  plotly_objects_result <- profile_plotly_objects()  # Reactive dependency on cached objects
  container_height <- input$container_height  # Reactive dependency on resize (from ResizeObserver)
  
  # Update cached height when user resizes (but not on initial load)
  if (!is.null(container_height) && !is.na(container_height) && is.numeric(container_height)) {
    cache_set("profile_panel_height", as.integer(container_height))
  }
  
  # Convert to integer for use in calculations
  container_height <- if (!is.null(container_height) && !is.na(container_height)) as.integer(container_height) else NULL
  

  req(plotly_objects_result)  # Don't proceed if objects not ready
  
  # Context already set by previous observer
  req(cxt_get_assembly())

  state$plotly_registered <- FALSE

  # Get combined plotly object using cached objects
  plot_result <- plot_profiles_cached(plotly_objects_result, container_height)
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
    # Make it reactive to container height changes
    container_height_reactive <- input$container_height
    
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
  
  # get full range from context
  xlim <- cxt_get_xlim()
  if (!is.null(xlim) && length(xlim) == 2) {
    return(xlim)
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

# refresh button functionality moved to server_navigation.r

