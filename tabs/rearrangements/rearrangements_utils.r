# tabs/rearrangements/rearrangements_utils.r
# standalone rearrangements processing utilities

# query rearrangements for a specific context (extracted from rearrangements_tab.r)
query_rearrangements_for_context <- function(assembly, contigs, zoom, tab_config) {
  # early return if no context
  if (is.null(assembly) || length(contigs) == 0) {
    return(NULL)
  }
  
  # extract required parameters from tab config
  get_aln_f <- tab_config$get_aln_f
  library_ids <- tab_config$library_ids
  max_margin <- tab_config$max_margin %||% 10
  min_element_length <- tab_config$min_element_length %||% 50
  min_anchor_length <- tab_config$min_anchor_length %||% 200
  max_anchor_mutations_percent <- tab_config$max_anchor_mutations_percent %||% 1.0
  max_element_mutation_percent <- tab_config$max_element_mutation_percent %||% 1.0
  
  # get alignment stores from all configured library_ids
  stores <- list()
  
  for (library_id in library_ids) {
    tryCatch({
      aln <- get_aln_f(assembly, library_id)
      if (!is.null(aln) && inherits(aln, "externalptr")) {
        stores[[library_id]] <- aln
      } else {
        warning(sprintf("get_aln_f returned invalid alignment for library_id '%s'", library_id))
      }
    }, error = function(e) {
      warning(sprintf("error getting alignment for library_id '%s': %s", library_id, e$message))
    })
  }
  
  if (length(stores) == 0) {
    return(NULL)
  }
  
  # get segments for selected contigs and build context for intervals
  all_segments <- get_segments(assembly)
  if (is.null(all_segments)) {
    return(NULL)
  }
  # get intervals from global context
  intervals <- cxt_get_zoom_view()
  if (is.null(intervals) || nrow(intervals) == 0) {
    return(NULL)
  }
  
  # get reference contigs for seam resolution
  reference_contigs <- get_fasta(assembly)
  if (is.null(reference_contigs)) {
    reference_contigs <- list()
  }
  
  # call aln_rearrange
  rearrange_results <- aln_rearrange(
    store_list = stores,
    intervals_df = intervals,
    max_margin = as.integer(max_margin),
    min_element_length = as.integer(min_element_length),
    min_anchor_length = as.integer(min_anchor_length),
    max_anchor_mutations_percent = as.numeric(max_anchor_mutations_percent),
    max_element_mutation_percent = as.numeric(max_element_mutation_percent),
    should_verify = FALSE,
    reference_contigs = reference_contigs,
    resolve_seams = "reference_only"
  )

  # add frequency column to events dataframe
  if (!is.null(rearrange_results$events) && nrow(rearrange_results$events) > 0) {
    events_df <- rearrange_results$events
    events_df$frequency <- ifelse(events_df$total_coverage > 0, 
                                 events_df$total_support / events_df$total_coverage, 0)
    rearrange_results$events <- events_df
  }
  
  return(rearrange_results)
}

# filter rearrangements by span (frequency range)
filter_rearrangements_by_span <- function(rearrange_data, min_span) {
  if (is.null(rearrange_data) || is.null(rearrange_data$events)) {
    return(rearrange_data)
  }
  
  events_df <- rearrange_data$events
  keep_events <- events_df$frequency >= min_span
  
  # filter all components
  filtered_data <- list(
    events = events_df[keep_events, , drop = FALSE],
    read_events = if (!is.null(rearrange_data$read_events)) {
      rearrange_data$read_events[rearrange_data$read_events$event_id %in% events_df$event_id[keep_events], , drop = FALSE]
    } else NULL,
    support = if (!is.null(rearrange_data$support)) rearrange_data$support[keep_events, , drop = FALSE] else NULL,
    coverage = if (!is.null(rearrange_data$coverage)) rearrange_data$coverage[keep_events, , drop = FALSE] else NULL,
    library_ids = rearrange_data$library_ids
  )
  
  return(filtered_data)
}

# filter rearrangements by total support
filter_rearrangements_by_support <- function(rearrange_data, min_support) {
  if (is.null(rearrange_data) || is.null(rearrange_data$events)) {
    return(rearrange_data)
  }
  
  events_df <- rearrange_data$events
  keep_events <- events_df$total_support >= min_support
  
  # filter all components
  filtered_data <- list(
    events = events_df[keep_events, , drop = FALSE],
    read_events = if (!is.null(rearrange_data$read_events)) {
      rearrange_data$read_events[rearrange_data$read_events$event_id %in% events_df$event_id[keep_events], , drop = FALSE]
    } else NULL,
    support = if (!is.null(rearrange_data$support)) rearrange_data$support[keep_events, , drop = FALSE] else NULL,
    coverage = if (!is.null(rearrange_data$coverage)) rearrange_data$coverage[keep_events, , drop = FALSE] else NULL,
    library_ids = rearrange_data$library_ids
  )
  
  return(filtered_data)
}

# add color information to rearrangements dataframe
add_rearrangement_colors <- function(events_df) {
  if (is.null(events_df) || nrow(events_df) == 0) {
    return(events_df)
  }
  
  # define colors for rearrangement types (matching profile colors)
  type_colors <- c(
    "large_insert" = "red",
    "large_delete" = "blue",  
    "large_invert" = "orange",
    "unknown" = "black"
  )
  
  # add color column
  events_df$color <- type_colors[events_df$type]
  events_df$color[is.na(events_df$color)] <- "black"  # default to black for unknown types
  
  return(events_df)
}

# filter rearrangements by region (contigs and zoom coordinates)
filter_rearrangements_by_region <- function(rearrange_data, contigs, zoom, assembly) {
  # early return if no data
  if (is.null(rearrange_data) || is.null(rearrange_data$events)) {
    return(NULL)
  }
  
  # early return if no filtering criteria
  if ((is.null(contigs) || length(contigs) == 0) && is.null(zoom)) {
    return(rearrange_data)
  }
  
  events_df <- rearrange_data$events
  
  # filter by contigs if specified
  if (!is.null(contigs) && length(contigs) > 0) {
    # check which column name is used for contig
    contig_col <- if ("contig" %in% colnames(events_df)) "contig" else "contig_id"
    
    keep_contigs <- events_df[[contig_col]] %in% contigs
    events_df <- events_df[keep_contigs, ]
    
    if (nrow(events_df) == 0) {
      return(NULL)
    }
    
    # filter support and coverage matrices
    if (!is.null(rearrange_data$support)) {
      rearrange_data$support <- rearrange_data$support[keep_contigs, , drop = FALSE]
    }
    if (!is.null(rearrange_data$coverage)) {
      rearrange_data$coverage <- rearrange_data$coverage[keep_contigs, , drop = FALSE]
    }
    
    # filter read_events if present
    if (!is.null(rearrange_data$read_events)) {
      kept_event_ids <- events_df$event_id
      rearrange_data$read_events <- rearrange_data$read_events[rearrange_data$read_events$event_id %in% kept_event_ids, ]
    }
  }
  
  # filter by zoom coordinates if specified
  if (!is.null(zoom) && length(zoom) == 2) {
    # need to convert local coordinates to global coordinates
    contigs_table <- get_contigs(assembly)
    if (!is.null(contigs_table)) {
      contig_col <- if ("contig" %in% colnames(events_df)) "contig" else "contig_id"
      
      # first check if coordinates are in plotted view before converting
      in_view <- cxt_coords_in_view(events_df[[contig_col]], events_df$out_clip, limit_to_zoom = FALSE)
      events_df <- events_df[in_view, ]
      
      if (nrow(events_df) == 0) {
        return(NULL)
      }
      
      # filter support and coverage matrices to match in_view filter
      if (!is.null(rearrange_data$support)) {
        rearrange_data$support <- rearrange_data$support[in_view, , drop = FALSE]
      }
      if (!is.null(rearrange_data$coverage)) {
        rearrange_data$coverage <- rearrange_data$coverage[in_view, , drop = FALSE]
      }
      
      # filter read_events if present
      if (!is.null(rearrange_data$read_events)) {
        kept_event_ids <- events_df$event_id
        rearrange_data$read_events <- rearrange_data$read_events[rearrange_data$read_events$event_id %in% kept_event_ids, ]
      }
      
      # convert event coordinates to global using context services
      events_df$gcoord <- cxt_contig2global(events_df[[contig_col]], events_df$out_clip)
      
      # filter by zoom range
      keep_zoom <- events_df$gcoord >= zoom[1] & events_df$gcoord <= zoom[2]
      events_df <- events_df[keep_zoom, ]
      
      # remove temporary gcoord column
      events_df$gcoord <- NULL
      
      if (nrow(events_df) == 0) {
        return(NULL)
      }
      
      # filter support and coverage matrices by zoom
      if (!is.null(rearrange_data$support)) {
        rearrange_data$support <- rearrange_data$support[keep_zoom, , drop = FALSE]
      }
      if (!is.null(rearrange_data$coverage)) {
        rearrange_data$coverage <- rearrange_data$coverage[keep_zoom, , drop = FALSE]
      }
      
      # filter read_events if present
      if (!is.null(rearrange_data$read_events)) {
        kept_event_ids <- events_df$event_id
        rearrange_data$read_events <- rearrange_data$read_events[rearrange_data$read_events$event_id %in% kept_event_ids, ]
      }
    }
  }
  
  # update events in the data structure
  rearrange_data$events <- events_df
  
  return(rearrange_data)
}

# save rearrangements table to tab-delimited file
save_rearrangements_table <- function(events_df, output_dir, region_name, use_simple_name = FALSE) {
  # create filename based on export type
  if (use_simple_name) {
    table_filename <- "rearrangements.txt"
  } else {
    table_filename <- paste0(region_name, "_rearrangements_table.txt")
  }
  filename <- file.path(output_dir, table_filename)
  
  if (is.null(events_df) || nrow(events_df) == 0) {
    # create empty file with headers
    empty_df <- data.frame(
      event_id = character(0),
      contig = character(0),
      out_clip = character(0),
      type = character(0),
      desc = character(0),
      frequency = numeric(0),
      stringsAsFactors = FALSE
    )
    write.table(empty_df, filename, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("saved empty rearrangements table to: %s\n", filename))
    return(filename)
  }
  
  # prepare table for export (remove internal columns like color)
  export_df <- events_df
  
  # remove color column if it exists
  if ("color" %in% names(export_df)) {
    export_df$color <- NULL
  }
  
  # save to file
  write.table(export_df, filename, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(sprintf("saved rearrangements table (%d events) to: %s\n", nrow(events_df), filename))
  return(filename)
}

# load rearrangements from files for static mode
load_rearrangements_from_files <- function(assembly, contigs, zoom, tab_config) {
  # early return if no assembly
  if (is.null(assembly)) {
    return(NULL)
  }
  
  # extract required parameters from tab config
  library_ids <- tab_config$library_ids
  get_rearrange_events_f <- tab_config$get_rearrange_events_f
  get_rearrange_support_f <- tab_config$get_rearrange_support_f
  get_rearrange_coverage_f <- tab_config$get_rearrange_coverage_f
  
  tryCatch({
    # load events table
    events_table <- get_rearrange_events_f(assembly)
    if (is.null(events_table) || nrow(events_table) == 0) {
      return(NULL)
    }
    
    # validate events table columns
    expected_events_cols <- c("event_id", "type", "contig_id", "out_clip", "in_clip")
    missing_events_cols <- expected_events_cols[!expected_events_cols %in% colnames(events_table)]
    if (length(missing_events_cols) > 0) {
      warning(sprintf("rearrangements events table missing expected columns: %s", 
                     paste(missing_events_cols, collapse = ", ")))
    }
    
    # check for summary statistics columns (may need to be calculated)
    summary_cols <- c("total_support", "total_coverage", "library_count", "frequency")
    missing_summary_cols <- summary_cols[!summary_cols %in% colnames(events_table)]
    if (length(missing_summary_cols) > 0) {
      warning(sprintf("rearrangements events table missing summary columns (will be calculated): %s", 
                     paste(missing_summary_cols, collapse = ", ")))
    }
    
    # load support and coverage matrices
    support_matrix <- get_rearrange_support_f(assembly)
    coverage_matrix <- get_rearrange_coverage_f(assembly)
    
    if (is.null(support_matrix) || is.null(coverage_matrix)) {
      return(NULL)
    }
    
    # validate matrix structures
    if (!"event_id" %in% colnames(support_matrix)) {
      warning("rearrangements support matrix missing 'event_id' column")
      return(NULL)
    }
    if (!"event_id" %in% colnames(coverage_matrix)) {
      warning("rearrangements coverage matrix missing 'event_id' column")
      return(NULL)
    }
    
    # check matrix dimensions match
    if (nrow(support_matrix) != nrow(coverage_matrix)) {
      warning(sprintf("rearrangements matrix dimension mismatch: support has %d rows, coverage has %d rows",
                     nrow(support_matrix), nrow(coverage_matrix)))
    }
    if (ncol(support_matrix) != ncol(coverage_matrix)) {
      warning(sprintf("rearrangements matrix dimension mismatch: support has %d columns, coverage has %d columns",
                     ncol(support_matrix), ncol(coverage_matrix)))
    }
    
    # filter events by selected contigs only (ignore zoom)
    # require contigs to be selected - if no contigs, show no events
    if (is.null(contigs) || length(contigs) == 0) {
      return(NULL)
    }
    
    # filter by selected contigs (assuming events have contig_id column)
    events_table <- events_table[events_table$contig_id %in% contigs, ]
    
    if (nrow(events_table) == 0) {
      return(NULL)
    }
    
    # filter support and coverage matrices to match filtered events
    event_ids <- events_table$event_id
    support_filtered <- support_matrix[support_matrix$event_id %in% event_ids, ]
    coverage_filtered <- coverage_matrix[coverage_matrix$event_id %in% event_ids, ]
    
    # warn if some events are missing from matrices
    missing_in_support <- event_ids[!event_ids %in% support_matrix$event_id]
    missing_in_coverage <- event_ids[!event_ids %in% coverage_matrix$event_id]
    
    if (length(missing_in_support) > 0) {
      warning(sprintf("rearrangements: %d event(s) from events table not found in support matrix", 
                     length(missing_in_support)))
    }
    if (length(missing_in_coverage) > 0) {
      warning(sprintf("rearrangements: %d event(s) from events table not found in coverage matrix", 
                     length(missing_in_coverage)))
    }
    
    # convert matrices to proper format (event_id as rownames, library columns only)
    rownames(support_filtered) <- support_filtered$event_id
    rownames(coverage_filtered) <- coverage_filtered$event_id
    
    # rename support/coverage matrix columns to match configured library_ids
    # skip first column (event_id) and map remaining columns to library_ids
    support_cols <- colnames(support_filtered)
    actual_lib_cols <- support_cols[support_cols != "event_id"]
    
    if (length(actual_lib_cols) == 0) {
      warning("no library columns found in support matrix (only event_id column)")
      return(NULL)
    }
    
    # map actual columns to configured library_ids (in order)
    num_libs <- min(length(actual_lib_cols), length(library_ids))
    if (num_libs < length(library_ids)) {
      warning(sprintf("fewer library columns (%d) than configured library_ids (%d)", 
                     length(actual_lib_cols), length(library_ids)))
    }
    
    # rename columns in both matrices
    old_names <- actual_lib_cols[1:num_libs]
    new_names <- library_ids[1:num_libs]
    
    for (i in 1:num_libs) {
      colnames(support_filtered)[colnames(support_filtered) == old_names[i]] <- new_names[i]
      colnames(coverage_filtered)[colnames(coverage_filtered) == old_names[i]] <- new_names[i]
    }
    
    available_libs <- new_names
    
    support_matrix_final <- as.matrix(support_filtered[, available_libs, drop = FALSE])
    coverage_matrix_final <- as.matrix(coverage_filtered[, available_libs, drop = FALSE])
    
    # reorder events table to match matrix row order
    events_table <- events_table[match(rownames(support_matrix_final), events_table$event_id), ]
    
    # calculate missing summary statistics from matrices
    if (!"total_support" %in% colnames(events_table)) {
      events_table$total_support <- rowSums(support_matrix_final, na.rm = TRUE)
    }
    if (!"total_coverage" %in% colnames(events_table)) {
      events_table$total_coverage <- rowSums(coverage_matrix_final, na.rm = TRUE)
    }
    if (!"library_count" %in% colnames(events_table)) {
      events_table$library_count <- rowSums(support_matrix_final > 0, na.rm = TRUE)
    }
    if (!"frequency" %in% colnames(events_table)) {
      events_table$frequency <- ifelse(events_table$total_coverage > 0, 
                                      events_table$total_support / events_table$total_coverage, 0)
    }
    
    # ensure column names match what profile expects
    if ("contig_id" %in% colnames(events_table) && !"contig" %in% colnames(events_table)) {
      events_table$contig <- events_table$contig_id
    }
    
    # return data in same format as query_rearrangements_for_context
    result <- list(
      events = events_table,
      support = support_matrix_final,
      coverage = coverage_matrix_final,
      library_ids = available_libs
    )
    
    # filter by zoom coordinates if specified
    if (!is.null(zoom) && length(zoom) == 2) {
      result <- filter_rearrangements_by_region(result, contigs, zoom, assembly)
    }
    
    return(result)
    
  }, error = function(e) {
    warning(sprintf("error loading rearrangements from files: %s", e$message))
    return(NULL)
  })
}
