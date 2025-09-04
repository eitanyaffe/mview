output$viewSelect <- renderUI({
  ids <- get_view_ids()
  selected_id <- isolate(input$view_id)
  if (is.null(selected_id)) {
    selected_id <- ids[1]
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
    set_view(input$view_id)
  },
  ignoreInit = TRUE
)
