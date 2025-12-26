# ---- State and Logging Management ----

# Get state segments data.frame
get_state_segments <- function() {
  return(state$segments)
}

# Get unique contigs from state segments
get_state_contigs <- function() {
  if (is.null(state$segments) || nrow(state$segments) == 0) {
    return(character())
  }
  return(unique(state$segments$contig))
}

# Get intervals from state segments, optionally merging adjacent segments
get_state_intervals <- function(merge_adjacent_segments = FALSE) {
  if (is.null(state$segments) || nrow(state$segments) == 0) {
    return(data.frame(contig = character(), start = integer(), end = integer(), stringsAsFactors = FALSE))
  }
  
  if (!merge_adjacent_segments) {
    return(data.frame(
      contig = state$segments$contig,
      start = state$segments$start,
      end = state$segments$end,
      stringsAsFactors = FALSE
    ))
  }
  
  # merge adjacent segments from same contig
  segments <- state$segments
  segments <- segments[order(segments$contig, segments$start), ]
  
  merged_intervals <- list()
  current_contig <- NULL
  current_start <- NULL
  current_end <- NULL
  
  for (i in seq_len(nrow(segments))) {
    seg <- segments[i, ]
    
    if (is.null(current_contig) || seg$contig != current_contig || seg$start != current_end + 1) {
      # start new interval
      if (!is.null(current_contig)) {
        merged_intervals[[length(merged_intervals) + 1]] <- list(
          contig = current_contig,
          start = current_start,
          end = current_end
        )
      }
      current_contig <- seg$contig
      current_start <- seg$start
      current_end <- seg$end
    } else {
      # extend current interval
      current_end <- seg$end
    }
  }
  
  # add last interval
  if (!is.null(current_contig)) {
    merged_intervals[[length(merged_intervals) + 1]] <- list(
      contig = current_contig,
      start = current_start,
      end = current_end
    )
  }
  
  # convert to data.frame
  if (length(merged_intervals) == 0) {
    return(data.frame(contig = character(), start = integer(), end = integer(), stringsAsFactors = FALSE))
  }
  
  do.call(rbind, lapply(merged_intervals, function(x) {
    data.frame(contig = x$contig, start = x$start, end = x$end, stringsAsFactors = FALSE)
  }))
}

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
  
  # convert to contig coordinates using context service
  tryCatch({
    local_info <- cxt_global2contig(plot_x)
    sprintf("Coord: %d", round(local_info$coord[1]))
  }, error = function(e) {
    sprintf("out of range, y:%.1f", plot_y)
  })
}

state <- reactiveValues(
  segments = data.frame(
    segment = character(),
    contig = character(),
    start = integer(),
    end = integer(),
    stringsAsFactors = FALSE
  ),
  zoom = NULL,
  current_xlim = NULL,
  assembly = get_assemblies()[1],
  clicked_read_alignments = NULL,
  mouse_coords = NULL
)

# reactive variable for last selected genome (gid)
# updated when user clicks goto on a genome
last_selected_genome <- reactiveVal(NULL)

# Project info output
output$project_info <- renderText({
  sprintf("Project: %s", get_project_id())
})

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
    
    # find which contigs are covered by zoom using cdf
    cdf <- cxt_get_entire_view()
    if (!is.null(cdf) && nrow(cdf) > 0) {
      zoom_contigs <- cdf[cdf$vstart < state$zoom[2] & cdf$vend > state$zoom[1], ]
      
      if (nrow(zoom_contigs) > 1) {
        zoom_text <- "Zoom: multiple contigs"
      } else if (nrow(zoom_contigs) == 1) {
        # single contig - convert to local coordinates
        contig_name <- zoom_contigs$contig[1]
        vstart <- zoom_contigs$vstart[1]
        local_start_offset <- zoom_contigs$start[1]
        local_start <- round(state$zoom[1] - vstart) + local_start_offset
        local_end <- round(state$zoom[2] - vstart) + local_start_offset
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

  segment_text <- sprintf("Segments: %d", nrow(get_state_segments()))
  list(
    tags$div(style = "margin-bottom: 4px;", assembly_text),
    tags$div(style = "margin-bottom: 4px;", segment_text),
    tags$div(style = "white-space: pre-wrap;", paste("  ", paste(get_state_segments()$segment, collapse = ", "))),
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

# update assembly from dropdown
observeEvent(input$`regions_module-assembly_select`, {
  if (!identical(state$assembly, input$`regions_module-assembly_select`)) {
    state$assembly <- input$`regions_module-assembly_select`
    cache_set("assembly.selected", state$assembly)
  }
})

# sync dropdown when assembly changes
observe({
  current_input <- input$`regions_module-assembly_select`
  selected <- state$assembly
  if (!is.null(selected) && !identical(current_input, selected)) {
    updateSelectInput(session, "regions_module-assembly_select", selected = selected)
  }
})

# persist assembly changes that originate outside the dropdown
observe({
  current_assembly <- state$assembly
  if (!is.null(current_assembly) && current_assembly != "") {
    cache_set("assembly.selected", current_assembly)
  }
})