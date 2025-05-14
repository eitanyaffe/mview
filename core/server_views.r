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
    id <- input$view_id
    set_view(id)
  },
  ignoreInit = TRUE
)
