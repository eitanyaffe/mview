# core/server_navigation.r
# Server-side logic for navigation panel

# navigation panel UI
output$navigationPanelUI <- renderUI({
  div(
    class = "navigation-panel",
    
    # refresh button (first on left)
    div(
      class = "nav-button-group",
      uiOutput("navRefreshBtnUI")
    ),
    
    # navigation buttons group
    div(
      class = "nav-button-group",
      actionButton("navZoomToWindowBtn", "", 
        icon = icon("crosshairs"),
        class = "btn-sm",
        title = "zoom to window (Alt+Z)"
      ),
      actionButton("navMoveLeftBtn", "", 
        icon = icon("arrow-left"),
        class = "btn-sm",
        title = "move left (Alt+Left)"
      ),
      actionButton("navMoveRightBtn", "", 
        icon = icon("arrow-right"),
        class = "btn-sm",
        title = "move right (Alt+Right)"
      ),
      actionButton("navZoomInBtn", "", 
        icon = icon("search-plus"),
        class = "btn-sm",
        title = "zoom in (Alt+Plus)"
      ),
      actionButton("navZoomOutBtn", "", 
        icon = icon("search-minus"),
        class = "btn-sm",
        title = "zoom out (Alt+Minus)"
      ),
      actionButton("navShowAllBtn", "", 
        icon = icon("expand"),
        class = "btn-sm",
        title = "show all (Alt+0)"
      )
    ),
    
    # context zoom buttons group
    div(
      class = "nav-button-group",
      actionButton("navZoom1kbBtn", "1kb", 
        class = "btn-sm context-zoom-btn",
        title = "zoom to 1kb window (Alt+3)"
      ),
      actionButton("navZoom10kbBtn", "10kb", 
        class = "btn-sm context-zoom-btn",
        title = "zoom to 10kb window (Alt+4)"
      ),
      actionButton("navZoom100kbBtn", "100kb", 
        class = "btn-sm context-zoom-btn",
        title = "zoom to 100kb window (Alt+5)"
      ),
      actionButton("navZoom1mbBtn", "1mb", 
        class = "btn-sm context-zoom-btn",
        title = "zoom to 1mb window (Alt+6)"
      )
    ),
    
    # undo button
    div(
      class = "nav-button-group",
      actionButton("navUndoBtn", "", 
        icon = icon("undo"),
        class = "btn-sm",
        title = "undo last action (Alt+Backspace)"
      )
    )
  )
})

# refresh button UI with conditional styling (moved from server_profiles.r)
output$navRefreshBtnUI <- renderUI({
  # observe plot_updated to determine button style
  is_updated <- plot_updated()
  btn_class <- if (is_updated) "btn-default" else "btn-warning"
  
  actionButton("navRefreshBtn", "", 
    icon = icon("refresh"), 
    class = paste("btn-sm", btn_class),
    title = "refresh plots"
  )
})

# helper function to get middle coordinate when zoom is null
get_middle_coordinate <- function(main_state_rv) {
  # get full range from context
  xlim <- cxt_get_xlim()
  
  if (!is.null(xlim) && length(xlim) == 2) {
    return((xlim[1] + xlim[2]) / 2)
  }
  
  # fallback
  return(50000)  # arbitrary center
}

# helper function to zoom to specific window size around center
zoom_to_window_size <- function(window_size, regions_module_output, main_state_rv) {
  # determine current center
  current_center <- if (!is.null(main_state_rv$zoom)) {
    # if zoomed, use center of current zoom
    (main_state_rv$zoom[1] + main_state_rv$zoom[2]) / 2
  } else {
    # if not zoomed, get middle coordinate of entire range
    get_middle_coordinate(main_state_rv)
  }
  
  if (is.null(current_center)) {
    cat("cannot determine center coordinate for zoom\n")
    return()
  }
  
  # save current region before changing zoom
  regions_module_output$push_undo_state()
  
  # set new zoom range centered on current center
  half_window <- window_size / 2
  new_start <- current_center - half_window
  new_end <- current_center + half_window
  main_state_rv$zoom <- c(new_start, new_end)
  
  cat(sprintf("zoomed to %s window around center %.0f\n", 
              format(window_size, big.mark = ","), current_center))
}

# navigation button handlers

# refresh button
observeEvent(input$navRefreshBtn, {
  current_val <- refresh_trigger()
  refresh_trigger(current_val + 1)
  cat("plot refresh triggered from navigation panel\n")
})

# zoom to window button
observeEvent(input$navZoomToWindowBtn, {
  if (!is.null(state$current_xlim)) {
    regions_module_output$push_undo_state()
    state$zoom <- state$current_xlim
    cat("zoomed to selected window\n")
  }
})

# move left button
observeEvent(input$navMoveLeftBtn, {
  if (!is.null(state$zoom)) {
    regions_module_output$push_undo_state()
    zoom_range <- state$zoom[2] - state$zoom[1]
    # move 3/4 window to the left
    shift <- zoom_range * 0.75
    state$zoom <- c(state$zoom[1] - shift, state$zoom[2] - shift)
    cat("moved zoom area left\n")
  }
})

# move right button
observeEvent(input$navMoveRightBtn, {
  if (!is.null(state$zoom)) {
    regions_module_output$push_undo_state()
    zoom_range <- state$zoom[2] - state$zoom[1]
    # move 3/4 window to the right
    shift <- zoom_range * 0.75
    state$zoom <- c(state$zoom[1] + shift, state$zoom[2] + shift)
    cat("moved zoom area right\n")
  }
})

# zoom in button
observeEvent(input$navZoomInBtn, {
  if (!is.null(state$zoom)) {
    regions_module_output$push_undo_state()
    zoom_range <- state$zoom[2] - state$zoom[1]
    # halve the range while keeping the same center
    state$zoom <- c(
      state$zoom[1] + zoom_range * 0.25,
      state$zoom[2] - zoom_range * 0.25
    )
    cat("zoomed in\n")
  }
})

# zoom out button
observeEvent(input$navZoomOutBtn, {
  if (!is.null(state$zoom)) {
    regions_module_output$push_undo_state()
    zoom_range <- state$zoom[2] - state$zoom[1]
    # double the range
    state$zoom <- c(
      state$zoom[1] - zoom_range * 0.5,
      state$zoom[2] + zoom_range * 0.5
    )
    cat("zoomed out\n")
  }
})

# show all button
observeEvent(input$navShowAllBtn, {
  regions_module_output$push_undo_state()
  state$zoom <- NULL
  cat("reset zoom to show all\n")
})

# context zoom buttons
observeEvent(input$navZoom1kbBtn, {
  zoom_to_window_size(1000, regions_module_output, state)
})

observeEvent(input$navZoom10kbBtn, {
  zoom_to_window_size(10000, regions_module_output, state)
})

observeEvent(input$navZoom100kbBtn, {
  zoom_to_window_size(100000, regions_module_output, state)
})

observeEvent(input$navZoom1mbBtn, {
  zoom_to_window_size(1000000, regions_module_output, state)
})

# undo button
observeEvent(input$navUndoBtn, {
  regions_module_output$undo_last_action()
})
