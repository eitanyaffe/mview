# ---- State and Logging Management ----

state <- reactiveValues(
  contigs = character(),
  zoom = NULL,
  current_xlim = NULL,
  assembly = get_assemblies()[1]
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

  assembly_text <- if (is.null(state$assembly) || state$assembly == "") {
    "Assembly: none selected"
  } else {
    sprintf("Assembly: %s", state$assembly)
  }

  contig_text <- sprintf("Contigs: %d", length(state$contigs))
  list(
    tags$div(style = "margin-bottom: 4px;", assembly_text),
    tags$div(style = "margin-bottom: 4px;", contig_text),
    tags$div(style = "white-space: pre-wrap;", paste("  ", paste(state$contigs, collapse = ", "))),
    tags$div(style = "margin-top: 6px;", zoom_text)
  )
})

# Add a new output renderer for the current state display box
output$current_state_display <- renderText({
  # Capture output from print_state in a string
  output_text <- capture.output({
    print_state(state, show_title = FALSE)
  })

  # Return the captured output as a single string
  paste(output_text, collapse = "\n")
})
