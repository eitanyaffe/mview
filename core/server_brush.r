
# Zoom logic
observeEvent(input$plot_dblclick, {
  b <- input$zoomBrush
  click_x <- input$plot_dblclick$x
  if (!is.null(b) && click_x >= b$xmin && click_x <= b$xmax) {
    state$zoom <- c(b$xmin, b$xmax)
    addLog(sprintf("Zoom set: %.1f - %.1f", b$xmin, b$xmax))
  } else {
    state$zoom <- NULL
    addLog("Zoom reset")
  }
})

# Hover info
output$profileHoverInfo <- renderText({
  hover_data <- reactiveValuesToList(input)
  for (id in names(profile_info_dfs)) {
    hname <- paste0("hover_", id)
    if (!is.null(hover_data[[hname]])) {
      h <- hover_data[[hname]]
      gcoord <- as.numeric(h$x)
      info_df <- profile_info_dfs[[id]]
      if (!is.null(info_df) && nrow(info_df) > 0) {
        idx <- which.min(abs(info_df$gcoord - gcoord))
        if (length(idx) == 1 && is.finite(info_df$gcoord[idx])) {
          return(info_df$description[idx])
        }
      }
    }
  }
  NULL
})
