# ---- State and Logging Management ----

state <- reactiveValues(
  contigs = character(),
  zoom = NULL,
  current_xlim = NULL,
  assembly = get_assemblies()[1],
  clicked_read_alignments = NULL
)

# Basic info output with same font as top state display
output$basic_info <- renderText({
  assembly_text <- if (is.null(state$assembly) || state$assembly == "") {
    "Assembly: none selected"
  } else {
    sprintf("Assembly: %s", state$assembly)
  }

  contig_len = length(state$contigs)
  contig_text <- sprintf("Contigs: n=%d", contig_len)
  
  contig_list <- if (contig_len > 0 && contig_len <= 10) {
    paste(paste(" ", state$contigs, collapse = "\n"))
  }

  # Combine all text with newlines
  paste(c(assembly_text, contig_text, contig_list), collapse = "\n")
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
