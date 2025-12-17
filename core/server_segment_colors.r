# server-side segment color functions
# observes state changes and registers color mapping functions

#########################################################################
# helper: compute segment degrees from filtered adjacency data
#########################################################################

compute_segment_degrees <- function(assembly, min_support_val, min_percent_val, lib = NULL) {
  # mode-aware loading
  mode <- tryCatch(graph_mode(), error = function(e) "segments")
  
  if (mode == "csegments") {
    count_mat <- load_cseg_adj_count(assembly)
    total_mat <- load_cseg_adj_total(assembly)
    associated_mat <- load_cseg_adj_associated(assembly)
  } else {
   count_mat <- load_seg_adj_count(assembly)
   total_mat <- load_seg_adj_total(assembly)
   associated_mat <- load_seg_adj_associated(assembly)
  }

  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(list())
  }
  
  if (!"src" %in% names(count_mat) || !"tgt" %in% names(count_mat)) {
    return(list())
  }
  
  edge_metrics <- aggregate_edge_metrics(count_mat, count_mat, total_mat, associated_mat, lib = lib)
  if (is.null(edge_metrics) || nrow(edge_metrics) == 0) {
    return(list())
  }
  
  filtered <- edge_metrics[edge_metrics$support >= min_support_val & 
                          edge_metrics$percent >= min_percent_val, ]
  
  if (nrow(filtered) == 0) {
    return(list())
  }
  
  # count unique neighbors (matching hover degree calculation)
  degrees_src <- sapply(split(filtered$tgt, filtered$src), function(x) length(unique(x)))
  degrees_tgt <- sapply(split(filtered$src, filtered$tgt), function(x) length(unique(x)))
  all_segments <- unique(c(filtered$src, filtered$tgt))
  rr <- data.frame(ids = all_segments)
  idx_src <- match(rr$ids, names(degrees_src))
  idx_tgt <- match(rr$ids, names(degrees_tgt))
  degrees <- (ifelse(is.na(idx_src), 0, degrees_src[idx_src]) + 
              ifelse(is.na(idx_tgt), 0, degrees_tgt[idx_tgt]))
  names(degrees) <- all_segments

  return(as.list(degrees))
}

#########################################################################
# helper: compute csegment degrees from filtered adjacency data
#########################################################################

compute_csegment_degrees <- function(assembly, min_support_val, min_percent_val, lib = NULL) {
  count_mat <- load_cseg_adj_count(assembly)
  total_mat <- load_cseg_adj_total(assembly)
  associated_mat <- load_cseg_adj_associated(assembly)
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(numeric())
  }
  if (!"src" %in% names(count_mat) || !"tgt" %in% names(count_mat)) {
    return(numeric())
  }
  edge_metrics <- aggregate_edge_metrics(count_mat, count_mat, total_mat, associated_mat, lib = lib)
  if (is.null(edge_metrics) || nrow(edge_metrics) == 0) {
    return(numeric())
  }
  filtered <- edge_metrics[edge_metrics$support >= min_support_val & 
                          edge_metrics$percent >= min_percent_val, ]
  if (nrow(filtered) == 0) {
    return(numeric())
  }
  # count unique neighbors (matching hover degree calculation)
  degrees_src <- sapply(split(filtered$tgt, filtered$src), function(x) length(unique(x)))
  degrees_tgt <- sapply(split(filtered$src, filtered$tgt), function(x) length(unique(x)))
  all_segments <- unique(c(filtered$src, filtered$tgt))
  rr <- data.frame(ids = all_segments)
  idx_src <- match(rr$ids, names(degrees_src))
  idx_tgt <- match(rr$ids, names(degrees_tgt))
  degrees <- (ifelse(is.na(idx_src), 0, degrees_src[idx_src]) + 
              ifelse(is.na(idx_tgt), 0, degrees_tgt[idx_tgt]))
  names(degrees) <- all_segments
  
  degrees
}

#########################################################################
# helper: degree-based color mapping functions
#########################################################################

# shared degree color thresholds
degree_color_thresholds <- list(
  gray = 2,      # 0-2: grey
  light_green = 6,  # 3-6: light green
  orange = 20    # 7-20: orange, >20: red
)

# shared function to map degree to color
degree_to_color <- function(deg) {
  ifelse(deg <= degree_color_thresholds$gray, "#E0E0E0",  # grey
    ifelse(deg <= degree_color_thresholds$light_green, "#90EE90",  # light green
      ifelse(deg <= degree_color_thresholds$orange, "#FF8C00",  # orange
        "#FF0000")))  # red
}

# color mapping for segment degrees
colors_by_segment_degree <- function(seg_degree_vec, ids) {
  deg <- ifelse(ids %in% names(seg_degree_vec), seg_degree_vec[ids], 0)
  degree_to_color(deg)
}

# color mapping for csegment degrees
colors_by_csegment_degree <- function(cseg_degree_vec, assembly, ids) {
  # IDs can be either segments or csegments - convert segments to csegments if needed
  ids_to_lookup <- ids
  mapping_table <- get_cluster_mapping(assembly)
  if (!is.null(mapping_table)) {
    are_segments <- ids %in% mapping_table$segment
    if (any(are_segments)) {
      ids_to_lookup[are_segments] <- mapping_table$csegment[match(ids[are_segments], mapping_table$segment)]
    }
    ids_to_lookup[is.na(ids_to_lookup)] <- ids[is.na(ids_to_lookup)]
  }
  
  deg <- ifelse(ids_to_lookup %in% names(cseg_degree_vec), cseg_degree_vec[ids_to_lookup], 0)
  degree_to_color(deg)
}

#########################################################################
# register color mappings when state changes
#########################################################################

# observe assembly, thresholds, segments - re-register all mappings
observe({
  assembly <- state$assembly
  segments <- state$segments
  min_support_val <- if (!is.null(input$selectorMinSupport)) input$selectorMinSupport else 0
  min_percent_val <- if (!is.null(input$selectorMinPercent)) input$selectorMinPercent else 0
  
  # get ordered segments for order schemes
  ordered_segs <- if (!is.null(segments) && nrow(segments) > 0) segments$segment else character()
  
  # register "order_rainbow" mapping
  if (length(ordered_segs) > 0) {
    rainbow_colors <- rainbow(length(ordered_segs), s = 0.6, v = 0.85)
    names(rainbow_colors) <- ordered_segs
    register_color_mapping_f("order_rainbow", function(ids) {
      result <- rep("#E0E0E0", length(ids))
      idx <- ids %in% names(rainbow_colors)
      result[idx] <- rainbow_colors[ids[idx]]
      result
    })
  } else {
    register_color_mapping_f("order_rainbow", function(ids) rep("#E0E0E0", length(ids)))
  }
  
  # register "order_grayscale" mapping
  if (length(ordered_segs) > 0) {
    gray_values <- seq(0.85, 0.3, length.out = length(ordered_segs))
    gray_colors <- gray(gray_values)
    names(gray_colors) <- ordered_segs
    register_color_mapping_f("order_grayscale", function(ids) {
      result <- rep("#E0E0E0", length(ids))
      idx <- ids %in% names(gray_colors)
      result[idx] <- gray_colors[ids[idx]]
      result
    })
  } else {
    register_color_mapping_f("order_grayscale", function(ids) rep("#E0E0E0", length(ids)))
  }
  
  # register degree mapping (both segment and csegment)
  if (!is.null(assembly)) {
    edge_lib <- if (!is.null(input$graphEdgeLibs)) input$graphEdgeLibs else "all"
    if (edge_lib == "all") edge_lib <- NULL
    
    # compute both segment and csegment degrees
    seg_degrees <- compute_segment_degrees(assembly, min_support_val, min_percent_val, lib = edge_lib)
    seg_degree_vec <- unlist(seg_degrees)
    cseg_degree_vec <- compute_csegment_degrees(assembly, min_support_val, min_percent_val, lib = edge_lib)
    
    # register the appropriate mapping based on current mode
    current_mode <- tryCatch(graph_mode(), error = function(e) "segments")
    if (current_mode == "csegments") {
      register_color_mapping_f("degree", function(ids) {
        colors_by_csegment_degree(cseg_degree_vec, assembly, ids)
      })
    } else {
      register_color_mapping_f("degree", function(ids) {
        colors_by_segment_degree(seg_degree_vec, ids)
      })
    }
  }
  
  # register user-defined scheme mappings
  user_schemes <- get_user_color_schemes()
  if (!is.null(assembly) && length(user_schemes) > 0) {
    bin_table <- load_bin_segment_table(assembly)
    for (scheme_name in names(user_schemes)) {
      color_field <- get_user_scheme_color_field(scheme_name)
      
      if (!is.null(color_field) && !is.null(bin_table) && color_field %in% names(bin_table)) {
        color_lookup <- setNames(bin_table[[color_field]], bin_table$segment)
        make_map_f <- function(lookup) {
          force(lookup)
          function(ids) {
            result <- rep("#E0E0E0", length(ids))
            idx <- ids %in% names(lookup)
            matched_cols <- lookup[ids[idx]]
            valid <- !is.na(matched_cols) & matched_cols != ""
            result[idx][valid] <- matched_cols[valid]
            result
          }
        }
        register_color_mapping_f(scheme_name, make_map_f(color_lookup))
      } else {
        register_color_mapping_f(scheme_name, function(ids) rep("#E0E0E0", length(ids)))
      }
    }
  }
})

#########################################################################
# direct degree color functions (bypass registry for graph nodes)
#########################################################################

get_csegment_degree_colors_direct <- function(csegment_ids, assembly, min_support_val, min_percent_val, lib = NULL) {
  if (length(csegment_ids) == 0) {
    return(setNames(character(0), character(0)))
  }
  
  cseg_degree_vec <- compute_csegment_degrees(assembly, min_support_val, min_percent_val, lib = lib)
  
  # directly map csegment IDs to colors (no conversion needed)
  deg <- ifelse(csegment_ids %in% names(cseg_degree_vec), cseg_degree_vec[csegment_ids], 0)
  colors <- degree_to_color(deg)
  
  return(setNames(colors, csegment_ids))
}

#########################################################################
# batch color function (for organizer and graph)
#########################################################################

get_segment_colors <- function(segment_ids) {
  if (length(segment_ids) == 0) {
    return(setNames(character(0), character(0)))
  }
  
  map_f <- get_current_color_map()
  colors <- map_f(segment_ids)
  
  return(setNames(colors, segment_ids))
}

#########################################################################
# observers for UI
#########################################################################

# update dropdown choices when assembly changes
observeEvent(state$assembly, {
  schemes <- get_segment_color_schemes()
  
  # load cached color scheme based on current mode
  mode <- tryCatch(graph_mode(), error = function(e) "segments")
  cache_key <- if (mode == "csegments") "csegment.color_scheme" else "segment.color_scheme"
  cached_scheme <- cache_get_if_exists(cache_key, NULL)
  
  current <- get_segment_color_scheme()
  if (!is.null(cached_scheme) && cached_scheme %in% schemes) {
    current <- cached_scheme
    set_segment_color_scheme(cached_scheme)
  } else if (!current %in% schemes) {
    set_segment_color_scheme("order_rainbow")
    current <- "order_rainbow"
  }
  
  updateSelectInput(session, "segmentColorScheme",
                   choices = schemes,
                   selected = current)
})

# set scheme when dropdown changes
observeEvent(input$segmentColorScheme, {
  if (!is.null(input$segmentColorScheme)) {
    set_segment_color_scheme(input$segmentColorScheme)
    
    # cache the color scheme based on current mode
    mode <- tryCatch(graph_mode(), error = function(e) "segments")
    cache_key <- if (mode == "csegments") "csegment.color_scheme" else "segment.color_scheme"
    cache_set(cache_key, input$segmentColorScheme)
    
    current_val <- refresh_trigger()
    refresh_trigger(current_val + 1)
  }
})

# update color scheme when mode changes
observeEvent(input$graphMode, {
  if (!is.null(input$graphMode)) {
    mode <- input$graphMode
    cache_key <- if (mode == "csegments") "csegment.color_scheme" else "segment.color_scheme"
    cached_scheme <- cache_get_if_exists(cache_key, NULL)
    
    if (!is.null(cached_scheme)) {
      schemes <- get_segment_color_schemes()
      if (cached_scheme %in% schemes) {
        set_segment_color_scheme(cached_scheme)
        updateSelectInput(session, "segmentColorScheme",
                         choices = schemes,
                         selected = cached_scheme)
      }
    }
  }
}, ignoreInit = TRUE)

