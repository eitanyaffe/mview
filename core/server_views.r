output$viewSelect <- renderUI({
  ids <- get_view_ids()
  selected_id <- isolate(input$view_id)
  
  # use cached recent view if no current selection
  if (is.null(selected_id)) {
    cached_view <- cache_get_if_exists("recent_view", NULL)
    if (!is.null(cached_view) && cached_view %in% ids) {
      selected_id <- cached_view
    } else {
      selected_id <- ids[1]
    }
  }
  
  selectInput(
    inputId = "view_id",
    label = "View:",
    choices = ids,
    selected = selected_id
  )
})

# set the view when changed
observeEvent(input$view_id,
  {
    req(input$view_id)
    # don't trigger set_view if we're updating programmatically from state load
    if (exists("updating_view_programmatically") && updating_view_programmatically()) {
      cat("skipping set_view - programmatic update\n")
      return()
    }
    
    # cache the selected view
    cache_set("recent_view", input$view_id)
    
    set_view(input$view_id)
  },
  ignoreInit = TRUE
)
