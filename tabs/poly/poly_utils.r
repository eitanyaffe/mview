# tabs/poly/poly_utils.r
# utility functions for poly tab data loading and filtering

# load unified data from getter functions
load_unified_data <- function(assembly, tab_config) {
  if (is.null(assembly)) {
    return(NULL)
  }
  
  # extract required parameters from tab config
  get_unify_table_f <- tab_config$get_unify_table_f
  get_unify_support_f <- tab_config$get_unify_support_f
  get_unify_coverage_f <- tab_config$get_unify_coverage_f
  library_ids <- tab_config$library_ids
  
  tryCatch({
    # load unified table
    unify_table <- get_unify_table_f(assembly)
    if (is.null(unify_table) || nrow(unify_table) == 0) {
      return(NULL)
    }
    
    # validate required columns
    required_cols <- c("uid", "aid", "type", "host_bin", "contig", 
    "coord1", "coord2", "is_sweeping", "desc")
    missing_cols <- required_cols[!required_cols %in% colnames(unify_table)]
    if (length(missing_cols) > 0) {
      warning(sprintf("unified table missing expected columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # load support and coverage matrices
    support_matrix <- get_unify_support_f(assembly)
    coverage_matrix <- get_unify_coverage_f(assembly)
    
    if (is.null(support_matrix) || is.null(coverage_matrix)) {
      return(NULL)
    }
    
    # validate matrix structures
    if (!"uid" %in% colnames(support_matrix)) {
      warning("unified support matrix missing 'uid' column")
      return(NULL)
    }
    if (!"uid" %in% colnames(coverage_matrix)) {
      warning("unified coverage matrix missing 'uid' column")
      return(NULL)
    }
    
    # match support and coverage to unified table by uid
    support_ix <- match(unify_table$uid, support_matrix$uid)
    coverage_ix <- match(unify_table$uid, coverage_matrix$uid)
    
    if (any(is.na(support_ix))) {
      missing <- unify_table$uid[is.na(support_ix)]
      warning(sprintf("unified table: %d uid(s) not found in support matrix", length(missing)))
    }
    if (any(is.na(coverage_ix))) {
      missing <- unify_table$uid[is.na(coverage_ix)]
      warning(sprintf("unified table: %d uid(s) not found in coverage matrix", length(missing)))
    }
    
    # get library columns from support matrix
    lib_cols <- setdiff(colnames(support_matrix), "uid")
    
    # align library columns between support and coverage
    common_libs <- intersect(lib_cols, setdiff(colnames(coverage_matrix), "uid"))
    if (length(common_libs) == 0) {
      warning("no common library columns found between support and coverage matrices")
      return(NULL)
    }
    
    # build support and coverage matrices aligned to unified table
    support_mat <- as.matrix(support_matrix[support_ix, common_libs, drop = FALSE])
    coverage_mat <- as.matrix(coverage_matrix[coverage_ix, common_libs, drop = FALSE])
    
    # set rownames to uid for easy matching
    rownames(support_mat) <- unify_table$uid
    rownames(coverage_mat) <- unify_table$uid
    
    # ensure column names match library_ids if provided
    if (!is.null(library_ids) && length(library_ids) > 0) {
      # map actual columns to configured library_ids
      num_libs <- min(length(common_libs), length(library_ids))
      if (num_libs < length(library_ids)) {
        warning(sprintf("fewer library columns (%d) than configured library_ids (%d)", 
                       length(common_libs), length(library_ids)))
      }
      
      # rename columns to match library_ids
      old_names <- common_libs[1:num_libs]
      new_names <- library_ids[1:num_libs]
      
      colnames(support_mat)[1:num_libs] <- new_names
      colnames(coverage_mat)[1:num_libs] <- new_names
      
      # keep only matched columns
      if (num_libs < length(common_libs)) {
        support_mat <- support_mat[, 1:num_libs, drop = FALSE]
        coverage_mat <- coverage_mat[, 1:num_libs, drop = FALSE]
      }
      
      available_libs <- new_names
    } else {
      available_libs <- common_libs
    }
    
    # return data structure
    result <- list(
      table = unify_table,
      support = support_mat,
      coverage = coverage_mat,
      library_ids = available_libs
    )
    
    return(result)
    
  }, error = function(e) {
    warning(sprintf("error loading unified data: %s", e$message))
    return(NULL)
  })
}

# filter unified data by mode (host, view, zoom, freeze)
filter_by_mode <- function(unified_data, mode, 
  host_bin = NULL, contigs = NULL, 
  zoom = NULL, 
  assembly = NULL, frozen_set = NULL) {
  if (is.null(unified_data) || is.null(unified_data$table)) {
    return(NULL)
  }
  
  # freeze mode: return frozen set if available
  if (mode == "freeze" && !is.null(frozen_set)) {
    return(frozen_set)
  }
  
  table_df <- unified_data$table
  
  # filter by mode
  if (mode == "host") {
    if (is.null(host_bin) || host_bin == "") {
      return(NULL)
    }
    keep_rows <- table_df$host_bin == host_bin
  } else if (mode == "view") {
    # view mode: filter by contigs from state$segments (all segments, not just plotted)
    if (is.null(contigs) || length(contigs) == 0) {
      return(NULL)
    }
    keep_rows <- table_df$contig %in% contigs
  } else if (mode == "zoom") {
    # zoom mode: use cxt_filter_intervals which filters to plotted segments AND by xlim
    # update context xlim to current zoom
    if (!is.null(zoom) && length(zoom) == 2) {
      cxt_set_zoom(zoom)
    }
    
    # prepare interval data for cxt_filter_intervals (needs: contig, start, end)
    interval_df <- data.frame(
      contig = table_df$contig,
      start = table_df$coord1,
      end = table_df$coord2,
      stringsAsFactors = FALSE
    )
    
    filtered_intervals <- cxt_filter_intervals(interval_df, merge_adjacent = FALSE)
    
    if (is.null(filtered_intervals) || nrow(filtered_intervals) == 0) {
      return(NULL)
    }
    
    # get the row indices that passed the filter by matching contig/coord pairs
    # create unique keys for matching
    table_keys <- paste(table_df$contig, table_df$coord1, table_df$coord2, sep = "|")
    filtered_keys <- paste(filtered_intervals$contig, filtered_intervals$start, filtered_intervals$end, sep = "|")
    keep_rows <- table_keys %in% filtered_keys
  } else {
    # freeze mode without frozen set, or unknown mode: return all
    keep_rows <- rep(TRUE, nrow(table_df))
  }
  
  if (sum(keep_rows) == 0) {
    return(NULL)
  }
  
  # filter all components
  filtered_data <- list(
    table = table_df[keep_rows, , drop = FALSE],
    support = unified_data$support[keep_rows, , drop = FALSE],
    coverage = unified_data$coverage[keep_rows, , drop = FALSE],
    library_ids = unified_data$library_ids
  )
  
  return(filtered_data)
}

# filter unified data by span (frequency range across libraries)
filter_by_span <- function(unified_data, min_span) {
  if (is.null(unified_data) || is.null(unified_data$table) || 
      is.null(unified_data$support) || is.null(unified_data$coverage)) {
    return(unified_data)
  }
  
  # validate min_span - must be a valid numeric value
  min_span_numeric <- suppressWarnings(as.numeric(min_span))
  if (is.na(min_span_numeric) || !is.finite(min_span_numeric)) {
    warning(sprintf("invalid min_span value '%s', using 0", as.character(min_span)))
    min_span <- 0
  } else {
    min_span <- min_span_numeric
  }
  
  # calculate frequency for each variant across libraries
  # support and coverage are matrices with rownames (uid) and numeric library columns
  support_data <- unified_data$support
  coverage_data <- unified_data$coverage
  
  # ensure we have matrices (they should already be matrices from load_unified_data)
  if (is.data.frame(support_data)) {
    # exclude uid column if present
    if ("uid" %in% names(support_data)) {
      support_matrix <- as.matrix(support_data[, names(support_data) != "uid", drop = FALSE])
    } else {
      support_matrix <- as.matrix(support_data)
    }
  } else {
    # already a matrix, just ensure it's numeric
    support_matrix <- support_data
  }
  
  if (is.data.frame(coverage_data)) {
    # exclude uid column if present
    if ("uid" %in% names(coverage_data)) {
      coverage_matrix <- as.matrix(coverage_data[, names(coverage_data) != "uid", drop = FALSE])
    } else {
      coverage_matrix <- as.matrix(coverage_data)
    }
  } else {
    # already a matrix, just ensure it's numeric
    coverage_matrix <- coverage_data
  }
  
  # ensure matrices are numeric and have same dimensions
  if (nrow(support_matrix) != nrow(coverage_matrix) || ncol(support_matrix) != ncol(coverage_matrix)) {
    warning("support and coverage matrices have different dimensions")
    return(unified_data)
  }
  
  # convert to numeric (handles any character/factor issues)
  # preserve dimensions and convert values
  support_numeric <- tryCatch({
    as.numeric(support_matrix)
  }, warning = function(w) {
    warning(sprintf("warning converting support to numeric: %s", w$message))
    as.numeric(support_matrix)
  }, error = function(e) {
    stop(sprintf("error converting support to numeric: %s", e$message))
  })
  
  coverage_numeric <- tryCatch({
    as.numeric(coverage_matrix)
  }, warning = function(w) {
    warning(sprintf("warning converting coverage to numeric: %s", w$message))
    as.numeric(coverage_matrix)
  }, error = function(e) {
    stop(sprintf("error converting coverage to numeric: %s", e$message))
  })
  
  support_matrix <- matrix(support_numeric, nrow = nrow(support_matrix), ncol = ncol(support_matrix))
  coverage_matrix <- matrix(coverage_numeric, nrow = nrow(coverage_matrix), ncol = ncol(coverage_matrix))
  
  # avoid division by zero - use matrix operations
  freq_matrix <- support_matrix / coverage_matrix
  freq_matrix[coverage_matrix <= 0 | is.na(coverage_matrix)] <- 0
  
  # calculate span (max - min frequency) for each variant
  variant_spans <- apply(freq_matrix, 1, function(row) {
    valid_freqs <- row[!is.na(row) & is.finite(row)]
    if (length(valid_freqs) == 0) return(0)
    max(valid_freqs) - min(valid_freqs)
  })
  
  # filter variants meeting span threshold
  keep_variants <- variant_spans >= min_span
  keep_variants[is.na(keep_variants)] <- FALSE
  
  if (sum(keep_variants, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  # filter all components
  filtered_data <- list(
    table = unified_data$table[keep_variants, , drop = FALSE],
    support = unified_data$support[keep_variants, , drop = FALSE],
    coverage = unified_data$coverage[keep_variants, , drop = FALSE],
    library_ids = unified_data$library_ids
  )
  
  return(filtered_data)
}

# calculate percent (support / coverage * 100) for plotting
calculate_percent <- function(unified_data) {
  if (is.null(unified_data) || is.null(unified_data$support) || is.null(unified_data$coverage)) {
    return(NULL)
  }
  
  support_matrix <- unified_data$support
  coverage_matrix <- unified_data$coverage
  
  # ensure matrices are numeric
  if (is.data.frame(support_matrix)) {
    support_matrix <- as.matrix(support_matrix)
  }
  if (is.data.frame(coverage_matrix)) {
    coverage_matrix <- as.matrix(coverage_matrix)
  }
  
  support_matrix <- matrix(as.numeric(support_matrix), nrow = nrow(support_matrix), ncol = ncol(support_matrix))
  coverage_matrix <- matrix(as.numeric(coverage_matrix), nrow = nrow(coverage_matrix), ncol = ncol(coverage_matrix))
  
  # avoid division by zero - use matrix operations
  percent_matrix <- (support_matrix / coverage_matrix) * 100
  percent_matrix[coverage_matrix <= 0 | is.na(coverage_matrix)] <- 0
  
  return(percent_matrix)
}

# split unified data by type into locals, rearrangements, elements
split_by_type <- function(unified_data) {
  if (is.null(unified_data) || is.null(unified_data$table)) {
    return(list(
      locals = NULL,
      rearrangements = NULL,
      elements = NULL
    ))
  }
  
  table_df <- unified_data$table
  
  # split by type
  locals_ix <- table_df$type == "local"
  rearrange_ix <- table_df$type == "rearrange"
  elements_ix <- table_df$type == "element"
  
  result <- list()
  
  # locals
  if (any(locals_ix)) {
    result$locals <- list(
      table = table_df[locals_ix, , drop = FALSE],
      support = unified_data$support[locals_ix, , drop = FALSE],
      coverage = unified_data$coverage[locals_ix, , drop = FALSE],
      library_ids = unified_data$library_ids
    )
  } else {
    result$locals <- NULL
  }
  
  # rearrangements
  if (any(rearrange_ix)) {
    result$rearrangements <- list(
      table = table_df[rearrange_ix, , drop = FALSE],
      support = unified_data$support[rearrange_ix, , drop = FALSE],
      coverage = unified_data$coverage[rearrange_ix, , drop = FALSE],
      library_ids = unified_data$library_ids
    )
  } else {
    result$rearrangements <- NULL
  }
  
  # elements
  if (any(elements_ix)) {
    result$elements <- list(
      table = table_df[elements_ix, , drop = FALSE],
      support = unified_data$support[elements_ix, , drop = FALSE],
      coverage = unified_data$coverage[elements_ix, , drop = FALSE],
      library_ids = unified_data$library_ids
    )
  } else {
    result$elements <- NULL
  }
  
  return(result)
}

