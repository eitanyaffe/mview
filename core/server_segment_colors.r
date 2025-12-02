# server-side segment color functions
# observes state changes and registers color mapping functions

#########################################################################
# helper: compute segment degrees from filtered adjacency data
#########################################################################

compute_segment_degrees <- function(assembly, min_support_val, min_percent_val) {
  count_mat <- load_seg_adj_count(assembly)
  total_mat <- load_seg_adj_total(assembly)
  associated_mat <- load_seg_adj_associated(assembly)

  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(list())
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(list())
  }
  
  edge_metrics <- aggregate_edge_metrics(count_mat, count_mat, total_mat, associated_mat)
  if (is.null(edge_metrics) || nrow(edge_metrics) == 0) {
    return(list())
  }
  
  filtered <- edge_metrics[edge_metrics$support >= min_support_val & 
                          edge_metrics$percent >= min_percent_val, ]
  
  if (nrow(filtered) == 0) {
    return(list())
  }
  
  degrees <- as.list(table(filtered$seg_src))
  return(degrees)
}

#########################################################################
# register color mappings when state changes
#########################################################################

# observe assembly, thresholds, and segments - re-register all mappings
observe({
  assembly <- state$assembly
  segments <- state$segments
  min_support_val <- if (!is.null(input$selectorMinSupport)) input$selectorMinSupport else 0
  min_percent_val <- if (!is.null(input$selectorMinPercent)) input$selectorMinPercent else 0
  
  # get ordered segments for "order" scheme
  ordered_segs <- if (!is.null(segments) && nrow(segments) > 0) segments$segment else character()
  
  # register "order" mapping
  if (length(ordered_segs) > 0) {
    rainbow_colors <- rainbow(length(ordered_segs), s = 0.6, v = 0.85)
    names(rainbow_colors) <- ordered_segs
    register_color_mapping_f("order", function(ids) {
      result <- rep("#E0E0E0", length(ids))
      idx <- ids %in% names(rainbow_colors)
      result[idx] <- rainbow_colors[ids[idx]]
      result
    })
  } else {
    register_color_mapping_f("order", function(ids) rep("#E0E0E0", length(ids)))
  }
  
  # register "degree" mapping
  if (!is.null(assembly)) {
    degrees <- compute_segment_degrees(assembly, min_support_val, min_percent_val)
    degree_vec <- unlist(degrees)
    register_color_mapping_f("degree", function(ids) {
      deg <- ifelse(ids %in% names(degree_vec), degree_vec[ids], 0)
      ifelse(deg == 0, "#E0E0E0", ifelse(deg == 1, "#4A90E2", ifelse(deg == 2, "#27AE60", "#E67E22")))
    })
  }
  
  # register user-defined scheme mappings
  user_schemes <- get_user_color_schemes()
  if (!is.null(assembly) && length(user_schemes) > 0) {
    bin_table <- load_bin_segment_table(assembly)
    for (scheme_name in names(user_schemes)) {
      field_name <- user_schemes[[scheme_name]]
      if (!is.null(bin_table) && field_name %in% names(bin_table)) {
        color_lookup <- setNames(bin_table[[field_name]], bin_table$segment)
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
# batch color function (for organizer and graph)
#########################################################################

get_segment_colors <- function(segment_ids, assembly, ordered_segments = NULL) {
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
  current <- get_segment_color_scheme()
  
  if (!current %in% schemes) {
    set_segment_color_scheme("order")
    current <- "order"
  }
  
  updateSelectInput(session, "segmentColorScheme",
                   choices = schemes,
                   selected = current)
})

# set scheme when dropdown changes
observeEvent(input$segmentColorScheme, {
  if (!is.null(input$segmentColorScheme)) {
    set_segment_color_scheme(input$segmentColorScheme)
    current_val <- refresh_trigger()
    refresh_trigger(current_val + 1)
  }
})

