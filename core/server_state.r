# ---- State and Logging Management ----

state <- reactiveValues(
  contigs = character(),
  zoom = NULL
)

log_messages <- reactiveVal(character())

addLog <- function(msg) {
  n <- input$log.length
  if (is.null(n) || n <= 0) n <- 20
  current <- log_messages()
  log_messages(c(msg, head(current, n - 1)))
}

output$log <- renderText({
  paste(rev(log_messages()), collapse = "\n")
})

observeEvent(input$clearLog, {
  log_messages(character())
})

output$contig_count <- renderUI({
  zoom_text <- if (is.null(state$zoom)) {
    "Zoom: full range"
  } else {
    sprintf("Zoom: %d – %d", state$zoom[1], state$zoom[2])
  }

  contig_text <- sprintf("Contigs: %d", length(state$contigs))
  list(
    tags$div(style = "margin-bottom: 4px;", contig_text),
    tags$div(style = "white-space: pre-wrap;", paste("  ", paste(state$contigs, collapse = ", "))),
    tags$div(style = "margin-top: 6px;", zoom_text)
  )
})
