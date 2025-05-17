#' Create an axis profile
#' @param id Unique profile id
#' @param name Display name
#' @param height Relative height weight
#' @param attr List of additional attributes
#' @param auto_register Whether to register immediately
#' @return Profile object
axis_profile <- function(id = "coord_axis",
                         name = "coordinate axis",
                         height = 0.5,
                         attr = list(),
                         auto_register = TRUE) {
  # data_f: no data needed for axis
  data_f <- function(cxt) {
    NULL
  }

  # plot_f: draw axes for each contig
  plot_f <- function(profile, cxt, gg) {
    if (is.null(cxt$mapper) || is.null(cxt$mapper$cdf) || nrow(cxt$mapper$cdf) == 0 || is.null(cxt$mapper$xlim)) {
      return(gg)
    }

    contig_df_all <- cxt$mapper$cdf
    zoom_xlim <- range(cxt$mapper$xlim)
    y_baseline <- 0

    # Filter contigs that are within the current zoom xlim
    contig_df <- contig_df_all[contig_df_all$start < zoom_xlim[2] & contig_df_all$end > zoom_xlim[1], , drop = FALSE]

    if (nrow(contig_df) == 0) {
      return(gg) # No contigs in the current view
    }

    axis_lines_list <- list()
    tick_marks_list <- list()
    tick_labels_list <- list()

    for (i in 1:nrow(contig_df)) {
      contig_info <- contig_df[i, ]
      # Adjust global start/end to be within the zoom window for drawing, but use original for local coord calculation
      draw_start <- max(contig_info$start, zoom_xlim[1])
      draw_end <- min(contig_info$end, zoom_xlim[2])
      contig_length <- contig_info$length
      # contig_name <- contig_info$cid # Not used if no hover

      # Axis line segments (clipped to zoom window)
      axis_lines_list[[i]] <- data.frame(
        x = draw_start, y = y_baseline, xend = draw_end, yend = y_baseline
      )

      # Generate ticks for local coordinates, but only draw those within the zoom window
      # Base tick positions on the original contig length and start for correct local coordinates
      tick_positions_local <- pretty(c(0, contig_length), n = 5) # n=5 is a suggestion, can be dynamic
      tick_positions_local <- tick_positions_local[tick_positions_local >= 0 & tick_positions_local <= contig_length]
      tick_positions_global <- contig_info$start + tick_positions_local

      # Filter ticks to be within the draw_start and draw_end of the current contig segment
      visible_ticks_idx <- tick_positions_global >= draw_start & tick_positions_global <= draw_end
      tick_positions_global_visible <- tick_positions_global[visible_ticks_idx]
      tick_positions_local_visible <- tick_positions_local[visible_ticks_idx]

      if (length(tick_positions_global_visible) > 0) {
        # Determine a sensible tick height, e.g., based on a small fraction of the y-range or a fixed value
        # Since y-range is not explicitly defined for axis, let's use a small absolute value or relative to x-axis range if dynamic
        # For simplicity, using a small fixed proportion of the x-axis zoom range. This might need adjustment.
        tick_height <- (zoom_xlim[2] - zoom_xlim[1]) * 0.005
        if (tick_height == 0) tick_height <- 0.1 # Prevent zero height if xlim range is tiny or zero

        # Tick mark segments
        tick_marks_list[[length(tick_marks_list) + 1]] <- data.frame(
          x = tick_positions_global_visible, y = y_baseline,
          xend = tick_positions_global_visible, yend = y_baseline - tick_height
        )

        # Tick labels
        tick_labels_list[[length(tick_labels_list) + 1]] <- data.frame(
          x = tick_positions_global_visible, y = y_baseline - tick_height * 2,
          label = tick_positions_local_visible
        )
      }
    }

    # Combine list of data.frames into single data.frames
    # Safely rbind only if lists are not empty
    final_axis_lines_df <- if (length(axis_lines_list) > 0) do.call(rbind, axis_lines_list) else NULL
    final_tick_marks_df <- if (length(tick_marks_list) > 0) do.call(rbind, tick_marks_list) else NULL
    final_tick_labels_df <- if (length(tick_labels_list) > 0) do.call(rbind, tick_labels_list) else NULL

    if (!is.null(final_axis_lines_df) && nrow(final_axis_lines_df) > 0) {
      gg <- gg + ggplot2::geom_segment(
        data = final_axis_lines_df,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "black", size = 0.5
      )
    }
    if (!is.null(final_tick_marks_df) && nrow(final_tick_marks_df) > 0) {
      gg <- gg + ggplot2::geom_segment(
        data = final_tick_marks_df,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "black", size = 0.5
      )
    }
    if (!is.null(final_tick_labels_df) && nrow(final_tick_labels_df) > 0) {
      gg <- gg + ggplot2::geom_text(
        data = final_tick_labels_df,
        ggplot2::aes(x = x, y = y, label = label),
        color = "black", size = 3, vjust = 1
      )
    }
    return(gg)
  }

  profile_create(
    id = id,
    name = name,
    type = "axis",
    height = height,
    attr = attr,
    data_f = data_f,
    plot_f = plot_f,
    auto_register = auto_register
  )
}
