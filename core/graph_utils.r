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

# finds neighbors connected by edges that pass min_support and min_percent filters
add_filtered_neighbors <- function(segment_id, min_support_val, min_percent_val, lib = NULL) {
  if (is.null(segment_id) || segment_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  count_mat <- load_seg_adj_count(state$assembly)
  total_mat <- load_seg_adj_total(state$assembly)
  associated_mat <- load_seg_adj_associated(state$assembly)
  
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(NULL)
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(NULL)
  }
  
  # find edges where segment is source or target (exclude self-loops)
  edge_rows <- count_mat[(count_mat$seg_src == segment_id | count_mat$seg_tgt == segment_id) &
                        count_mat$seg_src != count_mat$seg_tgt, ]
  
  if (nrow(edge_rows) == 0) {
    return(NULL)
  }
  
  # aggregate metrics for edges
  em <- aggregate_edge_metrics(edge_rows, count_mat, total_mat, associated_mat, lib = lib)
  if (is.null(em) || nrow(em) == 0) {
    return(NULL)
  }
  
  # apply filters
  em <- em[em$support >= min_support_val & em$percent >= min_percent_val, ]
  if (nrow(em) == 0) {
    return(NULL)
  }
  
  # extract neighbor segments (either src or tgt, excluding the current segment)
  neighbors <- unique(c(em$seg_src[em$seg_src != segment_id], 
                       em$seg_tgt[em$seg_tgt != segment_id]))
  
  return(neighbors)
}

# finds neighbor node IDs for selected node IDs, using adjacency matrices and filters
# returns node IDs (with _L/_R suffixes) for neighbors that pass min_support and min_percent
find_neighbor_nodes <- function(selected_node_ids, assembly, min_support_val, min_percent_val, lib = NULL) {
  if (length(selected_node_ids) == 0) {
    return(character())
  }
  
  # extract segment IDs from node IDs (strips _L or _R suffix)
  selected_seg_ids <- sub("_[LR]$", "", selected_node_ids)
  selected_seg_ids <- unique(selected_seg_ids)
  
  if (length(selected_seg_ids) == 0) {
    return(character())
  }
  
  # load adjacency matrices
  count_mat <- load_seg_adj_count(assembly)
  total_mat <- load_seg_adj_total(assembly)
  associated_mat <- load_seg_adj_associated(assembly)
  
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(character())
  }
  
  # find all edges from selected segments in the adjacency matrix
  edge_rows <- count_mat[count_mat$seg_src %in% selected_seg_ids | 
                        count_mat$seg_tgt %in% selected_seg_ids, ]
  
  if (nrow(edge_rows) == 0) {
    return(character())
  }
  
  # aggregate metrics and filter by thresholds
  em <- aggregate_edge_metrics(edge_rows, count_mat, total_mat, associated_mat, lib = lib)
  if (is.null(em) || nrow(em) == 0) {
    return(character())
  }
  
  # filter by min_support and min_percent
  filtered_em <- em[em$support >= min_support_val & em$percent >= min_percent_val, ]
  if (nrow(filtered_em) == 0) {
    return(character())
  }
  
  # extract neighbor segments (either src or tgt, excluding selected segments)
  neighbor_seg_ids <- unique(c(filtered_em$seg_src[!filtered_em$seg_src %in% selected_seg_ids],
                               filtered_em$seg_tgt[!filtered_em$seg_tgt %in% selected_seg_ids]))
  neighbor_seg_ids <- neighbor_seg_ids[!is.na(neighbor_seg_ids)]
  
  if (length(neighbor_seg_ids) == 0) {
    return(character())
  }
  
  # convert segment IDs to node IDs (both L and R for each neighbor segment)
  neighbor_nodes <- c(paste0(neighbor_seg_ids, "_L"), paste0(neighbor_seg_ids, "_R"))
  # remove any that are already selected
  neighbor_nodes <- neighbor_nodes[!neighbor_nodes %in% selected_node_ids]
  
  return(neighbor_nodes)
}

#########################################################################
# edge metrics aggregation
#########################################################################

aggregate_edge_metrics <- function(edge_rows, count_mat, total_mat, associated_mat, lib = NULL) {
  if (is.null(edge_rows) || nrow(edge_rows) == 0) {
    return(NULL)
  }
  
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  all_lib_cols <- names(count_mat)[!names(count_mat) %in% metadata_cols]
  key_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt")
  
  # determine which library columns to use
  if (is.null(lib) || lib == "" || lib == "all") {
    lib_cols <- all_lib_cols
  } else {
    if (lib %in% all_lib_cols) {
      lib_cols <- lib
    } else {
      lib_cols <- all_lib_cols
    }
  }
  browser()
  # row sums per matrix limited to library columns
  if (length(lib_cols) == 1) {
    count_sums <- data.frame(count_mat[, key_cols], 
                             support = count_mat[, lib_cols])
    total_sums <- data.frame(total_mat[, key_cols], 
                             total = total_mat[, lib_cols])
    associated_sums <- data.frame(associated_mat[, key_cols], 
                                  total_associated = associated_mat[, lib_cols])
  } else {
    count_sums <- data.frame(count_mat[, key_cols], 
                             support = rowSums(count_mat[, lib_cols, drop = FALSE], na.rm = TRUE))
    total_sums <- data.frame(total_mat[, key_cols], 
                             total = rowSums(total_mat[, lib_cols, drop = FALSE], na.rm = TRUE))
    associated_sums <- data.frame(associated_mat[, key_cols], 
                                  total_associated = rowSums(associated_mat[, lib_cols, drop = FALSE], na.rm = TRUE))
  }
  
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

# creates two side-nodes per segment (left and right)
build_graph_nodes <- function(segment_ids) {
  if (length(segment_ids) == 0) {
    return(data.frame(id = character(), label = character(), segment = character(), 
                      side = character(), stringsAsFactors = FALSE))
  }
  
  # create left and right nodes for each segment
  node_ids <- c(paste0(segment_ids, "_L"), paste0(segment_ids, "_R"))
  segments <- rep(segment_ids, 2)
  sides <- c(rep("L", length(segment_ids)), rep("R", length(segment_ids)))
  
  nodes <- data.frame(
    id = node_ids,
    label = node_ids,
    title = node_ids,
    segment = segments,
    side = sides,
    stringsAsFactors = FALSE
  )
  
  return(nodes)
}

# maps side values to L/R suffix
normalize_side <- function(side) {
  side_lower <- tolower(side)
  ifelse(side_lower %in% c("left", "l"), "L", 
         ifelse(side_lower %in% c("right", "r"), "R", side))
}

# creates internal edges (within segment) and inter-segment edges
build_graph_edges <- function(segment_ids, count_mat, total_mat, associated_mat, bin_segment_table, min_support_val, min_percent_val, directed = FALSE, lib = NULL) {
  empty_edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
  
  if (length(segment_ids) == 0) {
    return(empty_edges)
  }
  
  # internal edges connecting left and right sides of each segment
  internal_edges <- data.frame(
    id = paste0(segment_ids, "_internal"),
    from = paste0(segment_ids, "_L"),
    to = paste0(segment_ids, "_R"),
    title = segment_ids,
    ..is_internal.. = TRUE,
    ..segment.. = segment_ids,
    arrows = "",
    stringsAsFactors = FALSE
  )
  
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(internal_edges)
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(internal_edges)
  }
  
  # filter to edges where both src and tgt are in graph_segments (exclude self-loops)
  edge_rows <- count_mat[count_mat$seg_src %in% segment_ids & 
                        count_mat$seg_tgt %in% segment_ids &
                        count_mat$seg_src != count_mat$seg_tgt, ]
  
  if (nrow(edge_rows) == 0) {
    return(internal_edges)
  }
  
  # aggregate metrics for all directed edges
  em <- aggregate_edge_metrics(edge_rows, count_mat, total_mat, associated_mat, lib = lib)
  if (is.null(em) || nrow(em) == 0) {
    return(internal_edges)
  }
  
  # add bin information
  if (!is.null(bin_segment_table) && "segment" %in% names(bin_segment_table) && "bin" %in% names(bin_segment_table)) {
    seg_to_bin <- setNames(bin_segment_table$bin, bin_segment_table$segment)
    em$bin_src <- seg_to_bin[em$seg_src]
    em$bin_tgt <- seg_to_bin[em$seg_tgt]
    em$bin_src[is.na(em$bin_src)] <- "N/A"
    em$bin_tgt[is.na(em$bin_tgt)] <- "N/A"
  } else {
    em$bin_src <- "N/A"
    em$bin_tgt <- "N/A"
  }
  
  # normalize sides to L/R and create node IDs
  em$side_src_norm <- normalize_side(em$side_src)
  em$side_tgt_norm <- normalize_side(em$side_tgt)
  em$node_src <- paste0(em$seg_src, "_", em$side_src_norm)
  em$node_tgt <- paste0(em$seg_tgt, "_", em$side_tgt_norm)
  
  # apply filters
  em <- em[em$support >= min_support_val & em$percent >= min_percent_val, ]
  if (nrow(em) == 0) {
    return(internal_edges)
  }
  
  # store per-library data for color/label computation
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  lib_cols <- names(count_mat)[!names(count_mat) %in% metadata_cols]
  key_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt")
  
  # compute per-library percent for each edge in em
  # create lookup keys for matching
  em_keys <- paste(em$seg_src, em$seg_tgt, em$side_src, em$side_tgt, sep = "|")
  count_keys <- paste(count_mat$seg_src, count_mat$seg_tgt, 
                     count_mat$side_src, count_mat$side_tgt, sep = "|")
  assoc_keys <- paste(associated_mat$seg_src, associated_mat$seg_tgt,
                      associated_mat$side_src, associated_mat$side_tgt, sep = "|")
  
  lib_percent_data <- list()
  lib_support_data <- list()
  lib_assoc_data <- list()
  for (lib in lib_cols) {
    if (lib %in% names(count_mat) && lib %in% names(associated_mat)) {
      # match rows by keys
      count_idx <- match(em_keys, count_keys)
      assoc_idx <- match(em_keys, assoc_keys)
      
      lib_count <- count_mat[[lib]][count_idx]
      lib_assoc_val <- associated_mat[[lib]][assoc_idx]
      lib_count[is.na(lib_count)] <- 0
      lib_assoc_val[is.na(lib_assoc_val)] <- 0
      lib_percent <- ifelse(lib_assoc_val > 0, (lib_count / lib_assoc_val) * 100, 0)
      lib_percent_data[[lib]] <- lib_percent
      lib_support_data[[lib]] <- lib_count
      lib_assoc_data[[lib]] <- lib_assoc_val
    }
  }
  
  # build directed edges data frame
  inter_edges <- data.frame(
    id = paste0(em$node_src, "_", em$node_tgt),
    from = em$node_src,
    to = em$node_tgt,
    title = paste0(em$seg_src, " → ", em$seg_tgt),
    ..is_internal.. = FALSE,
    ..support.. = em$support,
    ..seg_from.. = em$seg_src,
    ..seg_to.. = em$seg_tgt,
    ..side_from.. = em$side_src_norm,
    ..side_to.. = em$side_tgt_norm,
    ..bin_from.. = em$bin_src,
    ..bin_to.. = em$bin_tgt,
    ..total.. = em$total,
    ..total_associated.. = em$total_associated,
    ..percent.. = em$percent,
    arrows = "to",
    stringsAsFactors = FALSE
  )
  
  # add per-library percent, support, and associated columns
  for (lib in lib_cols) {
    if (lib %in% names(lib_percent_data)) {
      inter_edges[[paste0("..pct_", lib, "..")]] <- lib_percent_data[[lib]]
    }
    if (lib %in% names(lib_support_data)) {
      inter_edges[[paste0("..sup_", lib, "..")]] <- lib_support_data[[lib]]
    }
    if (lib %in% names(lib_assoc_data)) {
      inter_edges[[paste0("..assoc_", lib, "..")]] <- lib_assoc_data[[lib]]
    }
  }
  
  # if undirected: merge pairs into single edges
  if (!directed) {
    inter_edges$key <- apply(inter_edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
    
    # for undirected, compute mean percent for each library across both directions
    merged_edges <- data.frame()
    for (key in unique(inter_edges$key)) {
      key_rows <- inter_edges[inter_edges$key == key, ]
      if (nrow(key_rows) == 1) {
        merged_edges <- rbind(merged_edges, key_rows)
      } else {
        # take first row as base, average per-library percents, support, and associated
        base_row <- key_rows[1, ]
        for (lib in lib_cols) {
          pct_col <- paste0("..pct_", lib, "..")
          sup_col <- paste0("..sup_", lib, "..")
          assoc_col <- paste0("..assoc_", lib, "..")
          if (pct_col %in% names(key_rows)) {
            base_row[[pct_col]] <- mean(key_rows[[pct_col]], na.rm = TRUE)
          }
          if (sup_col %in% names(key_rows)) {
            base_row[[sup_col]] <- mean(key_rows[[sup_col]], na.rm = TRUE)
          }
          if (assoc_col %in% names(key_rows)) {
            base_row[[assoc_col]] <- mean(key_rows[[assoc_col]], na.rm = TRUE)
          }
        }
        # also average the overall percent
        base_row$..percent.. <- mean(key_rows$..percent.., na.rm = TRUE)
        merged_edges <- rbind(merged_edges, base_row)
      }
    }
    inter_edges <- merged_edges
    inter_edges$id <- inter_edges$key
    inter_edges$title <- gsub(" → ", " - ", inter_edges$title)
    inter_edges$arrows <- ""
    inter_edges$key <- NULL
  }
  
  # combine internal and inter-segment edges
  for (col in names(inter_edges)) {
    if (!col %in% names(internal_edges)) internal_edges[[col]] <- NA
  }
  for (col in names(internal_edges)) {
    if (!col %in% names(inter_edges)) inter_edges[[col]] <- NA
  }
  internal_edges <- internal_edges[, names(inter_edges)]
  
  edges <- rbind(internal_edges, inter_edges)
  return(edges)
}

#########################################################################
# edge coloring and labeling
#########################################################################

# compute edge colors based on metric
compute_edge_colors <- function(edges, metric, lib1 = NULL, lib2 = NULL) {
  if (is.null(edges) || nrow(edges) == 0) {
    return(character())
  }
  
  colors <- rep("#000000", nrow(edges))
  
  # internal edges always light gray
  internal_mask <- !is.na(edges$..is_internal..) & edges$..is_internal..
  colors[internal_mask] <- "#D3D3D3"
  if (metric == "none" || all(internal_mask)) {
    return(colors)
  }
  
  inter_mask <- !internal_mask
  
  if (metric == "support") {
    support_vals <- edges$..support..[inter_mask]
    support_vals[is.na(support_vals)] <- 0
    max_support <- max(support_vals, na.rm = TRUE)
    if (max_support > 0) {
      # scale from light gray (#D3D3D3) to red (#FF0000)
      normalized <- support_vals / max_support
      normalized[normalized > 1] <- 1
      normalized[normalized < 0] <- 0
      for (i in seq_along(normalized)) {
        gray_val <- 211 - (normalized[i] * 211)  # 211 to 0
        red_val <- normalized[i] * 255
        colors[which(inter_mask)[i]] <- sprintf("#%02X%02X%02X", 
                                                 round(red_val), 
                                                 round(gray_val), 
                                                 round(gray_val))
      }
    } else {
      colors[inter_mask] <- "#D3D3D3"
    }
  } else if (metric == "percent") {
    percent_vals <- edges$..percent..[inter_mask]
    percent_vals[is.na(percent_vals)] <- 0
    # scale from 0 to 100, light gray to red
    normalized <- percent_vals / 100
    normalized[normalized > 1] <- 1
    normalized[normalized < 0] <- 0
    for (i in seq_along(normalized)) {
      gray_val <- 211 - (normalized[i] * 211)  # 211 to 0
      red_val <- normalized[i] * 255
      colors[which(inter_mask)[i]] <- sprintf("#%02X%02X%02X", 
                                               round(red_val), 
                                               round(gray_val), 
                                               round(gray_val))
    }
  } else if (metric == "change" && !is.null(lib1) && !is.null(lib2)) {
    pct_col1 <- paste0("..pct_", lib1, "..")
    pct_col2 <- paste0("..pct_", lib2, "..")
    
    if (pct_col1 %in% names(edges) && pct_col2 %in% names(edges)) {
      pct1 <- edges[[pct_col1]][inter_mask]
      pct2 <- edges[[pct_col2]][inter_mask]
      pct1[is.na(pct1)] <- 0
      pct2[is.na(pct2)] <- 0
      
      # add pseudo 1% to both before change calc to track changes with zeros
      pct1 <- pct1 + 1
      pct2 <- pct2 + 1
      
      # compute log2 change
      log2_change <- log2(pct2 / pct1)
      
      # check for infinity or NaN (can happen with very small/large ratios)
      log2_change[!is.finite(log2_change)] <- 0
      
      # clamp to [-2, 2], but preserve 0 for gray
      log2_change[log2_change > 2] <- 2
      log2_change[log2_change < -2] <- -2
      # ensure any remaining non-finite values are set to 0
      log2_change[!is.finite(log2_change)] <- 0
      
      # color: -2 (blue) to 0 (light gray) to +2 (red)
      for (i in seq_along(log2_change)) {
        val <- log2_change[i]
        if (val < 0) {
          # blue gradient: -2 = #0000FF, 0 = #D3D3D3 (light gray)
          normalized <- (val + 2) / 2  # 0 to 1
          r <- round(211 * normalized)  # 0 to 211
          g <- round(211 * normalized)  # 0 to 211
          b <- round(211 + (255 - 211) * normalized)  # 211 to 255
        } else {
          # red gradient: 0 = #D3D3D3 (light gray), +2 = #FF0000
          normalized <- val / 2  # 0 to 1
          r <- round(211 + (255 - 211) * normalized)  # 211 to 255
          g <- round(211 * (1 - normalized))  # 211 to 0
          b <- round(211 * (1 - normalized))  # 211 to 0
        }
        colors[which(inter_mask)[i]] <- sprintf("#%02X%02X%02X", r, g, b)
      }
    }
  }
  
  return(colors)
}

# compute edge labels based on metric
compute_edge_labels <- function(edges, metric, show_labels, lib1 = NULL, lib2 = NULL) {
  if (is.null(edges) || nrow(edges) == 0 || !show_labels) {
    return(rep("", nrow(edges)))
  }
  
  labels <- rep("", nrow(edges))
  internal_mask <- !is.na(edges$..is_internal..) & edges$..is_internal..
  inter_mask <- !internal_mask
  
  if (metric == "support") {
    support_vals <- edges$..support..[inter_mask]
    support_vals[is.na(support_vals)] <- 0
    labels[inter_mask] <- as.character(round(support_vals))
  } else if (metric == "percent") {
    percent_vals <- edges$..percent..[inter_mask]
    percent_vals[is.na(percent_vals)] <- 0
    labels[inter_mask] <- paste0(round(percent_vals), "%")
  } else if (metric == "change" && !is.null(lib1) && !is.null(lib2)) {
    # for change, show the log2 change value
    pct_col1 <- paste0("..pct_", lib1, "..")
    pct_col2 <- paste0("..pct_", lib2, "..")
    
    if (pct_col1 %in% names(edges) && pct_col2 %in% names(edges)) {
      pct1 <- edges[[pct_col1]][inter_mask]
      pct2 <- edges[[pct_col2]][inter_mask]
      pct1[is.na(pct1)] <- 0
      pct2[is.na(pct2)] <- 0
      
      # add pseudo 1% to both before change calc to track changes with zeros
      pct1 <- pct1 + 1
      pct2 <- pct2 + 1
      
      # compute log2 change
      log2_change <- log2(pct2 / pct1)
      
      # check for infinity or NaN
      log2_change[!is.finite(log2_change)] <- 0
      
      # clamp to [-2, 2] for display
      log2_change[log2_change > 2] <- 2
      log2_change[log2_change < -2] <- -2
      
      labels[inter_mask] <- sprintf("%.1f", log2_change)
    }
  }
  
  return(labels)
}

#########################################################################
# table building
#########################################################################

# builds table from inter-segment edges (excludes internal edges)
build_edge_table <- function(edges) {
  empty_table <- data.frame(
    seg_from = character(),
    seg_to = character(),
    bin_from = character(),
    bin_to = character(),
    side_from = character(),
    side_to = character(),
    support = numeric(),
    total_fwd = numeric(),
    total_associated_fwd = numeric(),
    percent_fwd = numeric(),
    total_rev = numeric(),
    total_associated_rev = numeric(),
    percent_rev = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (is.null(edges) || nrow(edges) == 0) {
    return(empty_table)
  }
  
  # filter to inter-segment edges only
  inter_edges <- edges[!edges$..is_internal.., ]
  if (nrow(inter_edges) == 0) {
    return(empty_table)
  }
  
  table_data <- data.frame(
    seg_from = inter_edges$..seg_from..,
    seg_to = inter_edges$..seg_to..,
    bin_from = inter_edges$..bin_from..,
    bin_to = inter_edges$..bin_to..,
    side_from = inter_edges$..side_from..,
    side_to = inter_edges$..side_to..,
    support = inter_edges$..support..,
    total_fwd = inter_edges$..total_fwd..,
    total_associated_fwd = inter_edges$..total_associated_fwd..,
    percent_fwd = inter_edges$..percent_fwd..,
    total_rev = inter_edges$..total_rev..,
    total_associated_rev = inter_edges$..total_associated_rev..,
    percent_rev = inter_edges$..percent_rev..,
    stringsAsFactors = FALSE
  )
  
  return(table_data)
}

