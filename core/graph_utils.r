# graph utilities
# data loading and graph building functions (pure functions, no Shiny dependencies)

#########################################################################
# color utilities
#########################################################################

# returns black or white text color based on background luminance
get_text_color_for_background <- function(bg_colors) {
  result <- rep("#000000", length(bg_colors))
  for (i in seq_along(bg_colors)) {
    col <- bg_colors[i]
    if (is.na(col) || col == "") {
      next
    }
    rgb_vals <- tryCatch(col2rgb(col), error = function(e) NULL)
    if (is.null(rgb_vals)) {
      next
    }
    luminance <- (0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]) / 255
    result[i] <- if (luminance < 0.5) "#FFFFFF" else "#000000"
  }
  result
}

#########################################################################
# data loading (uses registered functions from data.r)
#########################################################################

# registered function: register_seg_bins_f
load_bin_segment_table <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  get_seg_bins(assembly)
}

# registered function: register_seg_adj_f
load_seg_adj_count <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  adj <- get_seg_adj(assembly)
  if (is.null(adj)) return(NULL)
  adj$count
}

# registered function: register_seg_adj_f
load_seg_adj_total <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  adj <- get_seg_adj(assembly)
  if (is.null(adj)) return(NULL)
  adj$total
}

# registered function: register_seg_adj_f
load_seg_adj_associated <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  adj <- get_seg_adj(assembly)
  if (is.null(adj)) return(NULL)
  adj$associated
}

#########################################################################
# neighbor lookup
#########################################################################

add_neighbors <- function(segment_id) {
  if (is.null(segment_id) || segment_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  count_mat <- load_seg_adj_count(state$assembly)
  if (is.null(count_mat)) {
    return(NULL)
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(NULL)
  }
  
  # find neighbors where segment is source
  neighbors_src <- count_mat[count_mat$seg_src == segment_id, ]
  neighbors_from_src <- unique(neighbors_src$seg_tgt)
  
  # find neighbors where segment is target
  neighbors_tgt <- count_mat[count_mat$seg_tgt == segment_id, ]
  neighbors_from_tgt <- unique(neighbors_tgt$seg_src)
  
  # combine all neighbors
  all_neighbors <- unique(c(neighbors_from_src, neighbors_from_tgt))
  
  return(all_neighbors)
}

#########################################################################
# edge metrics aggregation
#########################################################################

aggregate_edge_metrics <- function(edge_rows, count_mat, total_mat, associated_mat) {
  if (is.null(edge_rows) || nrow(edge_rows) == 0) {
    return(NULL)
  }
  
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  lib_cols <- names(count_mat)[!names(count_mat) %in% metadata_cols]
  key_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt")
  
  # row sums per matrix limited to library columns
  count_sums <- data.frame(count_mat[, key_cols], 
                           support = rowSums(count_mat[, lib_cols, drop = FALSE], na.rm = TRUE))
  total_sums <- data.frame(total_mat[, key_cols], 
                           total = rowSums(total_mat[, lib_cols, drop = FALSE], na.rm = TRUE))
  associated_sums <- data.frame(associated_mat[, key_cols], 
                                total_associated = rowSums(associated_mat[, lib_cols, drop = FALSE], na.rm = TRUE))
  
  # merge with edge_rows
  result <- edge_rows[, key_cols]
  result <- merge(result, count_sums, by = key_cols, all.x = TRUE)
  result <- merge(result, total_sums, by = key_cols, all.x = TRUE)
  result <- merge(result, associated_sums, by = key_cols, all.x = TRUE)
  
  result$support[is.na(result$support)] <- 0
  result$total[is.na(result$total)] <- 0
  result$total_associated[is.na(result$total_associated)] <- 0
  result$percent <- ifelse(result$total_associated > 0, 
                           (result$support / result$total_associated) * 100, 0)
  
  return(result)
}

#########################################################################
# graph building
#########################################################################

build_graph_nodes <- function(segment_ids) {
  if (length(segment_ids) == 0) {
    return(data.frame(id = character(), label = character(), stringsAsFactors = FALSE))
  }
  
  nodes <- data.frame(
    id = segment_ids,
    label = segment_ids,
    title = segment_ids,
    stringsAsFactors = FALSE
  )
  
  return(nodes)
}

build_graph_edges <- function(segment_ids, count_mat, total_mat, associated_mat, bin_segment_table, min_support_val, min_percent_val) {
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  if (length(segment_ids) == 0) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # filter to edges where both src and tgt are in graph_segments
  edge_rows <- count_mat[count_mat$seg_src %in% segment_ids & 
                        count_mat$seg_tgt %in% segment_ids, ]
  
  if (nrow(edge_rows) == 0) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # aggregate metrics
  edge_metrics <- aggregate_edge_metrics(edge_rows, count_mat, total_mat, associated_mat)
  
  if (is.null(edge_metrics)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # add bin information
  if (!is.null(bin_segment_table) && "segment" %in% names(bin_segment_table) && "bin" %in% names(bin_segment_table)) {
    seg_to_bin <- setNames(bin_segment_table$bin, bin_segment_table$segment)
    edge_metrics$bin_src <- seg_to_bin[edge_metrics$seg_src]
    edge_metrics$bin_tgt <- seg_to_bin[edge_metrics$seg_tgt]
    edge_metrics$bin_src[is.na(edge_metrics$bin_src)] <- "N/A"
    edge_metrics$bin_tgt[is.na(edge_metrics$bin_tgt)] <- "N/A"
  } else {
    edge_metrics$bin_src <- "N/A"
    edge_metrics$bin_tgt <- "N/A"
  }
  
  # apply filters
  edge_metrics <- edge_metrics[edge_metrics$support >= min_support_val & 
                              edge_metrics$percent >= min_percent_val, ]
  
  if (nrow(edge_metrics) == 0) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # build edges data frame (include sides in id to ensure uniqueness)
  edges <- data.frame(
    id = paste0(edge_metrics$seg_src, "_", edge_metrics$side_src, "_", edge_metrics$seg_tgt, "_", edge_metrics$side_tgt),
    from = edge_metrics$seg_src,
    to = edge_metrics$seg_tgt,
    stringsAsFactors = FALSE
  )
  
  # simple tooltip - just segment IDs
  edges$title <- paste0(edge_metrics$seg_src, " - ", edge_metrics$seg_tgt)
  
  # store metrics in edge data for table
  edges$..support.. <- edge_metrics$support
  edges$..total.. <- edge_metrics$total
  edges$..total_associated.. <- edge_metrics$total_associated
  edges$..percent.. <- edge_metrics$percent
  edges$..side_src.. <- edge_metrics$side_src
  edges$..side_tgt.. <- edge_metrics$side_tgt
  edges$..bin_src.. <- edge_metrics$bin_src
  edges$..bin_tgt.. <- edge_metrics$bin_tgt
  
  return(edges)
}

#########################################################################
# table building
#########################################################################

build_edge_table <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0) {
    return(data.frame(
      seg_src = character(),
      seg_tgt = character(),
      bin_src = character(),
      bin_tgt = character(),
      side_src = character(),
      side_tgt = character(),
      support = numeric(),
      total = numeric(),
      total_associated = numeric(),
      percent = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  table_data <- data.frame(
    seg_src = edges$from,
    seg_tgt = edges$to,
    bin_src = edges$..bin_src..,
    bin_tgt = edges$..bin_tgt..,
    side_src = edges$..side_src..,
    side_tgt = edges$..side_tgt..,
    support = edges$..support..,
    total = edges$..total..,
    total_associated = edges$..total_associated..,
    percent = edges$..percent..,
    stringsAsFactors = FALSE
  )
  
  return(table_data)
}

