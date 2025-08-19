# ---- State and Logging Management ----

# Get mouse coordinate info for display
get_mouse_info <- function(state) {
  if (is.null(state$mouse_coords)) {
    return(NULL)
  }
  
  plot_x <- state$mouse_coords$plot_x
  plot_y <- state$mouse_coords$plot_y
  
  if (is.null(plot_x) || is.null(plot_y)) {
    return(NULL)
  }
  
  # build context for coordinate mapping
  cxt <- build_context(
    state_contigs = state$contigs,
    contig_table = get_contigs(state$assembly),
    zoom = state$zoom,
    assembly = state$assembly
  )
  
  if (is.null(cxt) || is.null(cxt$mapper)) {
    return(sprintf("plot(%.1f,%.1f)", plot_x, plot_y))
  }
  
  # convert to contig coordinates
  tryCatch({
    local_info <- cxt$mapper$g2l(plot_x)
    sprintf("Coord: %d", round(local_info$coord[1]))
  }, error = function(e) {
    sprintf("out of range, y:%.1f", plot_y)
  })
}

state <- reactiveValues(
  contigs = character(),
  zoom = NULL,
  current_xlim = NULL,
  assembly = get_assemblies()[1],
  clicked_read_alignments = NULL,
  mouse_coords = NULL
)

# Basic info output with same font as top region display
output$basic_info <- renderText({
  assembly_text <- if (is.null(state$assembly) || state$assembly == "") {
    "Assembly: none selected"
  } else {
    sprintf("Assembly: %s", state$assembly)
  }

  # zoom region information
  zoom_info <- if (is.null(state$zoom)) {
    "Zoom: full range"
  } else {
    window_size <- round(state$zoom[2]) - round(state$zoom[1])
    window_size_text <- format_bp(window_size)
    cxt <- build_context(
      state_contigs = state$contigs,
      contig_table = get_contigs(state$assembly),
      zoom = state$zoom,
      assembly = state$assembly)
    
    if (!is.null(cxt) && !is.null(cxt$mapper)) {
      # find which contigs are covered by zoom using mapper$cdf
      cdf <- cxt$mapper$cdf
      zoom_contigs <- cdf[cdf$start < state$zoom[2] & cdf$end > state$zoom[1], ]
      
      if (nrow(zoom_contigs) > 1) {
        zoom_text <- "Zoom: multiple contigs"
      } else if (nrow(zoom_contigs) == 1) {
        # single contig - convert to local coordinates
        contig_name <- zoom_contigs$contig[1]
        contig_start <- zoom_contigs$start[1]
        local_start <- round(state$zoom[1] - contig_start) + 1
        local_end <- round(state$zoom[2] - contig_start)
        zoom_text <- sprintf("%s: %d-%d", contig_name, local_start, local_end)
      } else {
        zoom_text <- sprintf("global: %d – %d", round(state$zoom[1]), round(state$zoom[2]))
      }
    } else {
      zoom_text <- sprintf("global: %d – %d", round(state$zoom[1]), round(state$zoom[2]))
    }
    
    c(zoom_text, sprintf("Window: %s", window_size_text))
  }

  # hover coord functionality disabled - coordinates are a bit off and need debugging
  # mouse_info <- get_mouse_info(state)
  
  # combine all text with newlines
  info_parts <- c(assembly_text, zoom_info)
  # if (!is.null(mouse_info)) {
  #   info_parts <- c(info_parts, mouse_info)
  # }
  paste(info_parts, collapse = "\n")
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

# Add a new output renderer for the current region display box
output$current_region_display <- renderText({
  # Capture output from print_state in a string
  output_text <- capture.output({
  print_region(state, show_title = FALSE)
  })

  # Return the captured output as a single string
  paste(output_text, collapse = "\n")
})