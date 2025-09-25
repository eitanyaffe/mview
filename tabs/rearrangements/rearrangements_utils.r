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
  max_anchor_mutations_percent <- tab_config$max_anchor_mutations_percent %||% 0.001
  max_element_mutation_percent <- tab_config$max_element_mutation_percent %||% 0.01
  
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
  
  # get contigs data and build context for intervals
  contigs_table <- get_contigs(assembly)
  if (is.null(contigs_table)) {
    return(NULL)
  }
  
  # build context for intervals
  cxt <- build_context(contigs, contigs_table, zoom, assembly)
  if (is.null(cxt) || is.null(cxt$intervals) || nrow(cxt$intervals) == 0) {
    return(NULL)
  }
  
  # call aln_rearrange
  rearrange_results <- aln_rearrange(
    store_list = stores,
    intervals_df = cxt$intervals,
    max_margin = as.integer(max_margin),
    min_element_length = as.integer(min_element_length),
    min_anchor_length = as.integer(min_anchor_length),
    max_anchor_mutations_percent = as.numeric(max_anchor_mutations_percent),
    max_element_mutation_percent = as.numeric(max_element_mutation_percent),
    should_verify = FALSE,
    reference_contigs = list(),
    reference_reads = list()
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
