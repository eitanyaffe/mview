interval_profile <- function(id, name, height = 60, is_fixed = TRUE,
                          intervals_f = NULL,
                          color = "#2E86AB",
                          color_field = NULL,
                          color_f = NULL,
                          merge_adjacent = TRUE,
                          auto_register = TRUE) 
{
  plot_f <- function(profile, gg) {
    assembly <- cxt_get_assembly()
    intervals <- NULL
    if (is.function(profile$intervals_f)) {
      intervals <- profile$intervals_f(assembly)
    } else if (is.data.frame(profile$intervals_f)) {
      intervals <- profile$intervals_f
    } else {
      warning(sprintf("invalid intervals_f for interval_profile '%s': %s", id, profile$intervals_f))
      return(list(plot = gg, legends = list()))
    }
    
    if (is.null(intervals) || nrow(intervals) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # filter intervals by current assembly first (if assembly column exists)
    if ("assembly" %in% colnames(intervals)) {
      assembly_intervals <- intervals[intervals$assembly == assembly, ]
      if (is.null(assembly_intervals) || nrow(assembly_intervals) == 0) {
        return(list(plot = gg, legends = list()))
      }
    } else {
      assembly_intervals <- intervals
    }

    # filter intervals using context
    filtered_intervals <- cxt_filter_intervals(assembly_intervals, merge_adjacent = profile$merge_adjacent)
    
    if (is.null(filtered_intervals) || nrow(filtered_intervals) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # determine fill colors - priority: color_f > color_field > color
    fill_colors <- NULL
    
    # try color_f first (function(id) -> color)
    if (!is.null(profile$color_f) && is.function(profile$color_f)) {
      if ("id" %in% names(filtered_intervals)) {
        fill_colors <- sapply(filtered_intervals$id, profile$color_f)
        fill_colors[is.na(fill_colors)] <- profile$color
      }
    }
    
    # fallback to color_field
    if (is.null(fill_colors) && !is.null(profile$color_field)) {
      color_field_name <- profile$color_field
      if (!grepl("_color$", color_field_name)) {
        color_field_name <- paste0(color_field_name, "_color")
      }
      if (color_field_name %in% names(filtered_intervals)) {
        fill_colors <- filtered_intervals[[color_field_name]]
        fill_colors[is.na(fill_colors) | fill_colors == ""] <- profile$color
      }
    }
    
    # fallback to default color
    if (is.null(fill_colors)) {
      fill_colors <- profile$color
    }
    
    # create hover text
    hover_text <- paste0(
      filtered_intervals$contig, ": ", 
      filtered_intervals$start, "-", filtered_intervals$end,
      "\n", filtered_intervals$desc
    )
    
    # add intervals layer as rectangles
    gg <- gg + 
      ggplot2::geom_rect(
        data = filtered_intervals,
        ggplot2::aes(
          xmin = gstart, xmax = gend,
          ymin = -0.15, ymax = 0.15,
          text = hover_text
        ),
        fill = fill_colors,
        color = "black",
        size = 0.5
      )
    xlim <- cxt_get_xlim()
    gg <- gg +
      ggplot2::geom_text(
        data = filtered_intervals,
        ggplot2::aes(
          x = gstart + diff(xlim) * 0.01,
          y = 0.25,
          label = id,
          text = hover_text
        ),
        color = "black",
        size = 2,
        hjust = 1,
        vjust = 0.5
      ) +
      ggplot2::ylim(-0.5, 0.6)
    
    return(list(plot = gg, legends = list()))
  }

  profile_create(
    id = id, name = name, type = "intervals", height = height, is_fixed = is_fixed,
    attr = list(hide_y_ticks = TRUE),
    params = NULL, plot_f = plot_f,
    intervals_f = intervals_f,
    color = color,
    color_field = color_field,
    color_f = color_f,
    merge_adjacent = merge_adjacent,
    auto_register = auto_register
  )
}

