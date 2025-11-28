# Context API - Service-oriented design
# Profiles and tabs call context services; they never create contexts.

# Private storage for contexts - PERSIST across reloads
if (!exists(".context_env")) {
  .context_env <- new.env(parent = emptyenv())
  .context_env$current_assembly <- NULL
  .context_env$contexts <- list()  # keyed by assembly ID
  .context_env$current_context <- NULL
  .context_env$current_segments <- NULL
  cat("[context.r] Created fresh context environment\n")
} else {
  cat("[context.r] Reusing existing context environment\n")
}

# Merge C++ result back to input df using input_index field
merge_context_result <- function(df_input, df_result) {
  if (is.null(df_input) || nrow(df_input) == 0) {
    return(df_result)
  }
  
  if (is.null(df_result) || nrow(df_result) == 0) {
    return(df_input)
  }
  
  if (!"input_index" %in% names(df_result)) {
    stop("merge_context_result: df_result missing required 'input_index' column")
  }
  
  input_ix <- df_result$input_index
  
  if (any(is.na(input_ix))) {
    stop("merge_context_result: input_index contains NA values")
  }
  
  if (any(input_ix < 1) || any(input_ix > nrow(df_input))) {
    stop(sprintf("merge_context_result: input_index out of range [1, %d]", nrow(df_input)))
  }
  
  # columns from result that should be added to merged output
  result_cols <- setdiff(names(df_result), c(names(df_input), "input_index"))
  
  # index directly into input (input_index is 1-based)
  merged <- df_input[input_ix, , drop = FALSE]
  
  for (col in result_cols) {
    merged[[col]] <- df_result[[col]]
  }
  
  return(merged)
}

################################################################################
# SETTERS
################################################################################

# Set current assembly and initialize context with contig/segment tables
cxt_set_assembly <- function(assembly_id) {
  if (is.null(assembly_id) || assembly_id == "") {
    .context_env$current_assembly <- NULL
    .context_env$current_context <- NULL
    .context_env$current_segments <- NULL
    cat("[cxt_set_assembly] Assembly ID is NULL or empty, context reset.\n")
    return(invisible(FALSE))
  }
  
  
  # Check if context already exists for this assembly
  if (!is.null(.context_env$contexts[[assembly_id]])) {
    cat(sprintf("[cxt_set_assembly] Reusing existing context for assembly '%s'\n", assembly_id))
    .context_env$current_assembly <- assembly_id
    .context_env$current_context <- .context_env$contexts[[assembly_id]]
  } else {
    cat(sprintf("[cxt_set_assembly] Creating new context for assembly '%s'\n", assembly_id))
    # Create new context pointer (once per assembly)
    .context_env$contexts[[assembly_id]] <- context_create()  # No parameters
    .context_env$current_assembly <- assembly_id
    .context_env$current_context <- .context_env$contexts[[assembly_id]]
  }
  
  # Always update tables (in case data changed)
  contig_table <- get_contigs(assembly_id)
  full_segments <- get_segments(assembly_id)
  
  # make sure contig table only contains contigs that are in the segments table
  ix = match(unique(full_segments$contig), contig_table$contig)
  if (any(is.na(ix))) {
    stop(sprintf("Cannot get contigs for assembly '%s'", assembly_id))
  }
  contig_table = contig_table[ix, ]

  if (is.null(contig_table) || nrow(contig_table) == 0) {
    stop(sprintf("Cannot get contigs for assembly '%s'", assembly_id))
  }
  if (is.null(full_segments) || nrow(full_segments) == 0) {
    stop(sprintf("Cannot get segments for assembly '%s'", assembly_id))
  }
  
  # Ensure full_segments has 'length' column
  if (!"length" %in% names(full_segments)) {
    full_segments$length <- full_segments$end - full_segments$start + 1
  }
  
  # Validate contig and segment table organization
  segment_contigs <- unique(full_segments$contig)
  contig_table_contigs <- contig_table$contig
  
  # Check that unique contigs in segments match contigs in contig table
  contigs_only_in_contig_table <- setdiff(contig_table_contigs, segment_contigs)
  contigs_only_in_segments_table <- setdiff(segment_contigs, contig_table_contigs)
  if (length(segment_contigs) != length(contig_table_contigs) || !all(segment_contigs == contig_table_contigs)) {
    stop(
      sprintf(
        paste0(
          "Assembly error: contig sets in contig table and segment table do not match.\n",
          " - Number of contigs in contig table: %d\n",
          " - Number of unique contigs in segments table: %d\n",
          " - Contigs in contig table only: %s\n",
          " - Contigs in segments table only: %s\n",
          "This usually indicates a problem in your assembly data.\n",
          "Please check that the 'contig' column in both the contig table and segment table have exactly matching entries."
        ),
        length(contig_table_contigs),
        length(segment_contigs),
        paste(contigs_only_in_contig_table, collapse = ", "),
        paste(contigs_only_in_segments_table, collapse = ", ")
      )
    )
  }
  
  # Check that total segment lengths match total contig lengths
  total_segment_length <- sum(full_segments$length)
  total_contig_length <- sum(contig_table$length)
  if (total_segment_length != total_contig_length) {
    stop("error: total segment lengths do not match total contig lengths")
  }
  
  cat(sprintf("[cxt_set_assembly] Setting tables: %d contigs, %d segments\n", 
              nrow(contig_table), nrow(full_segments)))
  
  # Initialize/update the context with current tables
  context_set_tables(.context_env$current_context, full_segments, contig_table)
  
  # Reset segments when assembly changes
  .context_env$current_segments <- NULL
  
  # Reset zoom to full range of the new assembly
  cxt_set_zoom(NULL)
  
  invisible(TRUE)
}

# Set selected segments for the current view
cxt_set_view <- function(selected_segments) {
  if (is.null(.context_env$current_assembly)) {
    warning("cxt_set_view: assembly not set, call cxt_set_assembly first")
    return(invisible(FALSE))
  }
  
  if (is.null(selected_segments) || nrow(selected_segments) == 0) {
    cat("[cxt_set_view] empty segments, clearing view\n")
    .context_env$current_segments <- data.frame(
      segment = character(),
      contig = character(),
      start = integer(),
      end = integer(),
      stringsAsFactors = FALSE
    )
    context_update_selected_segments(.context_env$current_context, 
                                    .context_env$current_segments, 
                                    numeric())
    return(invisible(TRUE))
  }
  
  # Validate columns
  if (!all(c("segment", "contig", "start", "end") %in% names(selected_segments))) {
    warning("cxt_set_view: selected_segments must have columns: segment, contig, start, end")
    return(invisible(FALSE))
  }
  
  # Filter to valid contigs
  contig_table <- get_contigs(.context_env$current_assembly)
  plotted_segments <- selected_segments[selected_segments$contig %in% contig_table$contig, ]
  
  cat(sprintf("[cxt_set_view] Setting view with %d segments (from %d input segments)\n", 
              nrow(plotted_segments), nrow(selected_segments)))
  
  # Store for cxt_set_zoom
  .context_env$current_segments <- plotted_segments
  
  # Update C++ context (no zoom yet)
  context_update_selected_segments(.context_env$current_context, plotted_segments, numeric())
  
  invisible(TRUE)
}

# Set zoom level (xlim) for the current view
cxt_set_zoom <- function(xlim) {
  if (is.null(.context_env$current_assembly)) {
    warning("cxt_set_zoom: assembly not set, call cxt_set_assembly first")
    return(invisible(FALSE))
  }
  
  if (is.null(.context_env$current_segments)) {
    warning("cxt_set_zoom: view not set, call cxt_set_view first")
    return(invisible(FALSE))
  }
  
  zoom_vec <- if (!is.null(xlim) && length(xlim) == 2) {
    as.numeric(xlim)
  } else {
    numeric()
  }
  
  if (length(zoom_vec) == 2) {
    cat(sprintf("[cxt_set_zoom] Setting zoom to [%.0f, %.0f]\n", zoom_vec[1], zoom_vec[2]))
  } else {
    cat("[cxt_set_zoom] Resetting zoom to full range\n")
  }
  
  # Update context with same segments but new zoom
  context_update_selected_segments(.context_env$current_context, 
                                  .context_env$current_segments, 
                                  zoom_vec)
  
  invisible(TRUE)
}

################################################################################
# GETTERS
################################################################################

# Get current assembly ID
cxt_get_assembly <- function() {
  .context_env$current_assembly
}

# Get current xlim (zoom range in vcoord)
cxt_get_xlim <- function() {
  if (is.null(.context_env$current_context)) {
    warning("cxt_get_xlim: no current context")
    return(c(0, 0))
  }
  
  context_get_xlim(.context_env$current_context)
}

# Get entire view intervals (all plotted segments, no xlim clipping)
cxt_get_entire_view <- function(merge_adjacent = TRUE) {
  if (is.null(.context_env$current_context)) {
    warning("cxt_get_entire_view: no current context")
    return(data.frame(vstart = numeric(), vend = numeric(), contig = character(), 
                      start = integer(), end = integer(), segment_ids = character(), 
                      n_segments = integer()))
  }
  
  context_get_view_intervals(.context_env$current_context, FALSE, merge_adjacent)
}

# Get zoom view intervals (clipped to xlim)
cxt_get_zoom_view <- function(merge_adjacent = TRUE) {
  if (is.null(.context_env$current_context)) {
    warning("cxt_get_zoom_view: no current context")
    return(data.frame(vstart = numeric(), vend = numeric(), contig = character(), 
                      start = integer(), end = integer(), segment_ids = character(), 
                      n_segments = integer()))
  }
  
  context_get_view_intervals(.context_env$current_context, TRUE, merge_adjacent)
}

# Get raw plotted segments from C++ (low-level)
cxt_get_plotted_segments <- function() {
  if (is.null(.context_env$current_context)) {
    warning("cxt_get_plotted_segments: no current context")
    return(data.frame())
  }
  
  context_get_plotted_segments(.context_env$current_context)
}

# Get unique contigs in current view
cxt_get_contigs <- function() {
  if (is.null(.context_env$current_segments)) {
    warning("cxt_get_contigs: no current segments")
    return(character())
  }
  
  unique(.context_env$current_segments$contig)
}

# Get segment IDs in current view
cxt_get_segments <- function() {
  if (is.null(.context_env$current_segments)) {
    warning("cxt_get_segments: no current segments")
    return(character())
  }
  
  .context_env$current_segments$segment
}

################################################################################
# COORD TRANSFORMERS
################################################################################

# Convert contig coords to global vcoords
cxt_contig2global <- function(contigs, coords) {
  if (is.null(.context_env$current_context)) {
    warning("cxt_contig2global: no current context")
    return(rep(NA_real_, length(contigs)))
  }
  
  df <- data.frame(contig = contigs, coord = coords, stringsAsFactors = FALSE)
  df_minimal <- df[c("contig", "coord")]
  result_minimal <- context_contig2view_point(.context_env$current_context, df_minimal)
  result <- merge_context_result(df, result_minimal)
  result$vcoord
}

# Convert global vcoords to contig coords
cxt_global2contig <- function(gcoords) {
  if (is.null(.context_env$current_context)) {
    warning("cxt_global2contig: no current context")
    return(data.frame(contig = character(), coord = integer(), gcoord = numeric()))
  }
  
  df <- data.frame(vcoord = gcoords)
  df_minimal <- df[c("vcoord")]
  result_minimal <- context_view2contig_point(.context_env$current_context, df_minimal)
  result <- merge_context_result(df, result_minimal)
  data.frame(
    contig = result$contig,
    coord = result$coord,
    gcoord = gcoords,
    stringsAsFactors = FALSE
  )
}

################################################################################
# INTERVAL FUNCTIONS
################################################################################

# Convert contig intervals to vcoord intervals
cxt_contig2view_interval <- function(df, crossing_intervals = "drop") {
  if (is.null(.context_env$current_context)) {
    warning("cxt_contig2view_interval: no current context")
    return(data.frame())
  }
  
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }
  
  df_minimal <- df[c("contig", "start", "end")]
  result_minimal <- context_contig2view_interval(.context_env$current_context, df_minimal, crossing_intervals)
  merge_context_result(df, result_minimal)
}

# Filter point coords to current xlim and add vcoord
cxt_filter_coords <- function(df) {
  if (is.null(.context_env$current_context)) {
    warning("cxt_filter_coords: no current context")
    return(NULL)
  }
  
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  
  if (!all(c("contig", "coord") %in% names(df))) {
    warning("cxt_filter_coords: df must have columns: contig, coord")
    return(NULL)
  }
  
  df_minimal <- df[c("contig", "coord")]
  xlim <- cxt_get_xlim()
  xlim_vec <- if (length(xlim) == 2) as.numeric(xlim) else numeric()
  result_minimal <- context_filter_coords(.context_env$current_context, df_minimal, xlim_vec)
  
  if (is.null(result_minimal) || nrow(result_minimal) == 0) {
    return(NULL)
  }
  
  result <- merge_context_result(df, result_minimal)
  
  # Add gcoord for backward compatibility
  if ("vcoord" %in% names(result)) {
    result$gcoord <- result$vcoord
  }
  
  result
}

# Filter interval segments to current xlim and add vstart/vend
cxt_filter_segments <- function(df) {
  if (is.null(.context_env$current_context)) {
    warning("cxt_filter_segments: no current context")
    return(NULL)
  }
  
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  
  if (!all(c("contig", "start", "end") %in% names(df))) {
    warning("cxt_filter_segments: df must have columns: contig, start, end")
    return(NULL)
  }
  
  df_minimal <- df[c("contig", "start", "end")]
  xlim <- cxt_get_xlim()
  result_minimal <- context_filter_segments(.context_env$current_context, df_minimal, xlim)
  
  if (is.null(result_minimal) || nrow(result_minimal) == 0) {
    return(NULL)
  }
  
  result <- merge_context_result(df, result_minimal)
  
  # Add gstart/gend for backward compatibility
  if ("vstart" %in% names(result)) {
    result$gstart <- result$vstart
  }
  if ("vend" %in% names(result)) {
    result$gend <- result$vend
  }
  
  result
}

################################################################################
# VIEW INTERVALS
################################################################################

