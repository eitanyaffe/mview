# server-side segment color functions
# observes state changes and registers color mapping functions

#########################################################################
# helper: convert hex color to grayscale
#########################################################################

color_to_grayscale <- function(hex_color) {
  hex_color <- gsub("^#", "", hex_color)
  if (nchar(hex_color) != 6) return("#888888")
  
  r <- strtoi(substr(hex_color, 1, 2), base = 16)
  g <- strtoi(substr(hex_color, 3, 4), base = 16)
  b <- strtoi(substr(hex_color, 5, 6), base = 16)
  
  gray <- round(0.299 * r + 0.587 * g + 0.114 * b)
  gray_hex <- sprintf("#%02X%02X%02X", gray, gray, gray)
  return(gray_hex)
}

#########################################################################
# helper: compute segment degrees from filtered adjacency data
#########################################################################

compute_segment_degrees <- function(assembly, min_support_val, min_percent_val, lib = NULL) {
  count_mat <- load_seg_adj_count(assembly)
  total_mat <- load_seg_adj_total(assembly)
  associated_mat <- load_seg_adj_associated(assembly)

  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(list())
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
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
  
  degrees <- as.list(table(filtered$seg_src))
  return(degrees)
}

#########################################################################
# register color mappings when state changes
#########################################################################

# observe assembly, thresholds, segments, and grayscale checkbox - re-register all mappings
observe({
  assembly <- state$assembly
  segments <- state$segments
  min_support_val <- if (!is.null(input$selectorMinSupport)) input$selectorMinSupport else 0
  min_percent_val <- if (!is.null(input$selectorMinPercent)) input$selectorMinPercent else 0
  use_grayscale <- if (!is.null(input$segmentGrayscale)) input$segmentGrayscale else FALSE
  
  # get ordered segments for "order" scheme
  ordered_segs <- if (!is.null(segments) && nrow(segments) > 0) segments$segment else character()
  
  # register "order" mapping
  if (length(ordered_segs) > 0) {
    if (use_grayscale) {
      # create light-to-dark grayscale gradient based on order
      gray_values <- seq(0.85, 0.3, length.out = length(ordered_segs))
      gray_colors <- gray(gray_values)
      names(gray_colors) <- ordered_segs
      register_color_mapping_f("order", function(ids) {
        result <- rep("#E0E0E0", length(ids))
        idx <- ids %in% names(gray_colors)
        result[idx] <- gray_colors[ids[idx]]
        result
      })
    } else {
      rainbow_colors <- rainbow(length(ordered_segs), s = 0.6, v = 0.85)
      names(rainbow_colors) <- ordered_segs
      register_color_mapping_f("order", function(ids) {
        result <- rep("#E0E0E0", length(ids))
        idx <- ids %in% names(rainbow_colors)
        result[idx] <- rainbow_colors[ids[idx]]
        result
      })
    }
  } else {
    register_color_mapping_f("order", function(ids) rep("#E0E0E0", length(ids)))
  }
  
  # register "degree" mapping
  if (!is.null(assembly)) {
    edge_lib <- if (!is.null(input$graphEdgeLibs)) input$graphEdgeLibs else "all"
    if (edge_lib == "all") edge_lib <- NULL
    degrees <- compute_segment_degrees(assembly, min_support_val, min_percent_val, lib = edge_lib)
    degree_vec <- unlist(degrees)
    if (use_grayscale) {
      # create light-to-dark grayscale gradient based on degree
      max_deg <- if (length(degree_vec) > 0) max(degree_vec) else 0
      register_color_mapping_f("degree", function(ids) {
        deg <- ifelse(ids %in% names(degree_vec), degree_vec[ids], 0)
        if (max_deg == 0) {
          return(rep("#E0E0E0", length(ids)))
        }
        # map degree 0 to light gray, higher degrees to darker
        gray_levels <- ifelse(deg == 0, 0.9, 0.9 - (deg / max_deg) * 0.6)
        gray_levels <- pmax(0.3, gray_levels)
        colors <- gray(gray_levels)
        colors
      })
    } else {
      register_color_mapping_f("degree", function(ids) {
        deg <- ifelse(ids %in% names(degree_vec), degree_vec[ids], 0)
        ifelse(deg == 0, "#E0E0E0", ifelse(deg == 1, "#4A90E2", ifelse(deg == 2, "#27AE60", "#E67E22")))
      })
    }
  }
  
  # register user-defined scheme mappings
  user_schemes <- get_user_color_schemes()
  if (!is.null(assembly) && length(user_schemes) > 0) {
    bin_table <- load_bin_segment_table(assembly)
    for (scheme_name in names(user_schemes)) {
      scheme_def <- user_schemes[[scheme_name]]
      color_field <- scheme_def$color
      gray_field <- scheme_def$gray
      
      if (use_grayscale && !is.null(gray_field) && !is.null(bin_table) && gray_field %in% names(bin_table)) {
        # use gray field if available and grayscale is checked
        color_lookup <- setNames(bin_table[[gray_field]], bin_table$segment)
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
      } else if (!is.null(color_field) && !is.null(bin_table) && color_field %in% names(bin_table)) {
        # use color field, create grayscale gradient if needed
        color_lookup <- setNames(bin_table[[color_field]], bin_table$segment)
        if (use_grayscale) {
          # create light-to-dark grayscale gradient based on unique color values
          unique_colors <- unique(bin_table[[color_field]])
          unique_colors <- unique_colors[!is.na(unique_colors) & unique_colors != ""]
          if (length(unique_colors) > 0) {
            # sort colors by their luminance to maintain relative ordering
            color_luminance <- function(hex) {
              hex <- gsub("^#", "", hex)
              if (nchar(hex) != 6) return(0.5)
              r <- strtoi(substr(hex, 1, 2), base = 16)
              g <- strtoi(substr(hex, 3, 4), base = 16)
              b <- strtoi(substr(hex, 5, 6), base = 16)
              (0.299 * r + 0.587 * g + 0.114 * b) / 255
            }
            luminances <- vapply(unique_colors, color_luminance, numeric(1))
            sorted_colors <- unique_colors[order(luminances, decreasing = TRUE)]
            gray_values <- seq(0.85, 0.3, length.out = length(sorted_colors))
            gray_colors <- gray(gray_values)
            names(gray_colors) <- sorted_colors
            # create mapping from original colors to grayscale
            gray_lookup <- setNames(gray_colors[bin_table[[color_field]]], bin_table$segment)
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
            register_color_mapping_f(scheme_name, make_map_f(gray_lookup))
          } else {
            register_color_mapping_f(scheme_name, function(ids) rep("#E0E0E0", length(ids)))
          }
        } else {
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
        }
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

# refresh profiles when grayscale checkbox changes
observeEvent(input$segmentGrayscale, {
  if (exists("refresh_trigger")) {
    current_val <- refresh_trigger()
    refresh_trigger(current_val + 1)
  }
}, ignoreInit = TRUE)

