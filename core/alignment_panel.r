# core/alignment_panel.r

# Handle plotly click events from the combined plot to display alignment details
observeEvent(eventExpr = plotly::event_data("plotly_click"), {
  cat("plotly_click event detected\n")
  event_data <- plotly::event_data("plotly_click")
  req(event_data)
  req(event_data$key)

  read_ids <- event_data$key[[1]]
  aln <- cache_get("aln_obj")
  cat(sprintf("read_ids: %s\n", paste(read_ids, collapse = ", ")))

  df <- aln_alignments_from_read_ids(aln, read_ids)
  if (!is.null(df) && nrow(df) > 0) {
    state$clicked_read_alignments <- df
  } else {
    state$clicked_read_alignments <- NULL
  }
})

# Render the DataTable for alignment details
output$alignmentDetailsTable <- renderDT({
  align_data <- state$clicked_read_alignments
  if (!is.null(align_data) && nrow(align_data) > 0) {
    # Add a unique row identifier if it doesn't exist, for button ID generation
    # Ensure row_id is unique even if align_data is filtered/subset later
    # However, for this specific case, simple sequence should be fine as it's regenerated on click
    align_data$row_id <- seq_len(nrow(align_data))

    # Add the "Go to" button column
    align_data$`Go to` <- sapply(align_data$row_id, function(id) {
      # The inputId for setInputValue must be unique and not clash with other inputs
      as.character(actionButton(paste0("goto_aln_btn_", id), "Go",
        onclick = sprintf("Shiny.setInputValue(\"goto_alignment_location_trigger\", %d, {priority: \"event\"})", id)
      ))
    })

    datatable(
      align_data[, setdiff(colnames(align_data), "row_id")], # Exclude row_id from display
      escape = FALSE, # Important to render HTML for the button
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        lengthMenu = c(5, 10, 25, 50),
        dom = "lftip"
      ),
      selection = "none"
    )
  } else {
    # Return an empty datatable or a message if no data
    datatable(data.frame(Message = "Click on an alignment in the plot to see details here."), rownames = FALSE)
  }
})

# Observer for the "Go to" button clicks
observeEvent(input$goto_alignment_location_trigger, {
  req(input$goto_alignment_location_trigger) # Reacts to the button click
  row_idx <- as.integer(input$goto_alignment_location_trigger)

  # Ensure that clicked_read_alignments is not NULL and row_idx is valid
  req(state$clicked_read_alignments, row_idx <= nrow(state$clicked_read_alignments), row_idx > 0)

  selected_alignment <- state$clicked_read_alignments[row_idx, ]

  # Ensure the necessary columns exist in the selected alignment data
  req(
    selected_alignment,
    "contig_id" %in% colnames(selected_alignment),
    "contig_start" %in% colnames(selected_alignment),
    "contig_end" %in% colnames(selected_alignment)
  )

  states_module_output$push_state()

  state$contigs <- selected_alignment$contig_id
  gstart <- selected_alignment$contig_start
  gend <- selected_alignment$contig_end
  dd <- (gend - gstart) / 2
  state$zoom <- c(gstart - dd, gend + dd)
})
