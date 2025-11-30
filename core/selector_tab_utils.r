# selector tab utilities
# data loading and graph building functions (pure functions, no Shiny dependencies)

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
  
  # get library columns (exclude metadata columns)
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  
  # find library columns in count_mat
  lib_cols <- names(count_mat)[!names(count_mat) %in% metadata_cols]
  
  # aggregate metrics for each edge
  result <- data.frame(
    seg_src = edge_rows$seg_src,
    seg_tgt = edge_rows$seg_tgt,
    side_src = edge_rows$side_src,
    side_tgt = edge_rows$side_tgt,
    support = numeric(nrow(edge_rows)),
    total = numeric(nrow(edge_rows)),
    total_associated = numeric(nrow(edge_rows)),
    percent = numeric(nrow(edge_rows)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(edge_rows))) {
    seg_src <- edge_rows$seg_src[i]
    seg_tgt <- edge_rows$seg_tgt[i]
    side_src <- edge_rows$side_src[i]
    side_tgt <- edge_rows$side_tgt[i]
    
    # find matching rows in all three matrices
    count_row <- count_mat[count_mat$seg_src == seg_src & 
                          count_mat$seg_tgt == seg_tgt &
                          count_mat$side_src == side_src &
                          count_mat$side_tgt == side_tgt, ]
    total_row <- total_mat[total_mat$seg_src == seg_src & 
                          total_mat$seg_tgt == seg_tgt &
                          total_mat$side_src == side_src &
                          total_mat$side_tgt == side_tgt, ]
    associated_row <- associated_mat[associated_mat$seg_src == seg_src & 
                                    associated_mat$seg_tgt == seg_tgt &
                                    associated_mat$side_src == side_src &
                                    associated_mat$side_tgt == side_tgt, ]
    
    # aggregate across libraries
    support_sum <- 0
    total_sum <- 0
    associated_sum <- 0
    
    if (nrow(count_row) > 0) {
      for (col in lib_cols) {
        if (col %in% names(count_row)) {
          support_sum <- support_sum + as.numeric(count_row[[col]][1])
        }
      }
    }
    
    if (nrow(total_row) > 0) {
      for (col in lib_cols) {
        if (col %in% names(total_row)) {
          total_sum <- total_sum + as.numeric(total_row[[col]][1])
        }
      }
    }
    
    if (nrow(associated_row) > 0) {
      for (col in lib_cols) {
        if (col %in% names(associated_row)) {
          associated_sum <- associated_sum + as.numeric(associated_row[[col]][1])
        }
      }
    }
    
    result$support[i] <- support_sum
    result$total[i] <- total_sum
    result$total_associated[i] <- associated_sum
    result$percent[i] <- if (!is.na(associated_sum) && associated_sum > 0) (support_sum / associated_sum) * 100 else 0
  }
  
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
  
  # build edges data frame
  edges <- data.frame(
    id = paste0(edge_metrics$seg_src, "_", edge_metrics$seg_tgt),
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
