interval_profile <- function(id, name, height = 60, is_fixed = TRUE,
                         intervals_f = NULL,
                           color = "#2E86AB",
                           color_field = NULL,
                           merge_adjacent = TRUE,
                           auto_register = TRUE) {

  plot_f <- function(profile, gg) {
    assembly <- cxt_get_assembly()
    intervals <- NULL
    if (is.function(profile$intervals_f)) {
      intervals <- profile$intervals_f(assembly)
    } else if (is.character(profile$intervals_f)) {
      if (cache_exists(profile$intervals_f)) {
        intervals <- cache_get(profile$intervals_f)
      }
    } else if (is.data.frame(profile$intervals_f)) {
      intervals <- profile$intervals_f
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
    
    # determine color field name if color_field is provided
    color_field_name <- NULL
    if (!is.null(profile$color_field)) {
      color_field_name <- profile$color_field
      if (!grepl("_color$", color_field_name)) {
        color_field_name <- paste0(color_field_name, "_color")
      }
    }
    
    # determine fill colors
    if (!is.null(color_field_name) && color_field_name %in% names(filtered_intervals)) {
      fill_colors <- filtered_intervals[[color_field_name]]
      fill_colors[is.na(fill_colors) | fill_colors == ""] <- profile$color
    } else {
      fill_colors <- profile$color
    }
    
    # create hover text
    hover_text <- paste0(
      filtered_intervals$contig, ": ", 
      filtered_intervals$start, "-", filtered_intervals$end,
      "\n", filtered_intervals$desc
    )
    
    # add intervals layer as rectangles
    if (!is.null(color_field_name) && color_field_name %in% names(filtered_intervals)) {
      gg <- gg + 
        ggplot2::geom_rect(
          data = filtered_intervals,
          ggplot2::aes(
            xmin = gstart, xmax = gend,
            ymin = -0.15, ymax = 0.15,
            text = hover_text,
            fill = .data[[color_field_name]]
          ),
          color = "black",
          size = 0.5
        ) +
        ggplot2::scale_fill_identity()
    } else {
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
    }
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
    merge_adjacent = merge_adjacent,
    auto_register = auto_register
  )
}

