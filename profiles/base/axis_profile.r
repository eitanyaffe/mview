#' Create an axis profile
#' @param id Unique profile id
#' @param name Display name
#' @param height Relative height weight
#' @param attr List of additional attributes
#' @param auto_register Whether to register immediately
#' @return Profile object
axis_profile <- function(id = "coord_axis",
                         name = "coordinate axis",
                         height = 0.1,
                         attr = list(),
                         auto_register = TRUE) {
  # Set axis-specific attributes
  axis_attr <- list(
    hide_y_label = TRUE, # Hide y-label for axis profile
    hide_y_ticks = TRUE # Hide y-ticks for axis profile
  )

  # Merge with user-provided attributes, keeping axis-specific overrides
  attr <- c(axis_attr, attr)

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
    contig_names_list <- list() # list to store contig name labels

    for (i in seq_len(nrow(contig_df))) {
      contig_info <- contig_df[i, ]
      # Adjust global start/end to be within the zoom window for drawing, but use original for local coord calculation
      draw_start <- max(contig_info$start, zoom_xlim[1])
      draw_end <- min(contig_info$end, zoom_xlim[2])
      contig_length <- contig_info$length
      contig_name <- contig_info$contig # Get contig name

      # Axis line segments (clipped to zoom window)
      axis_lines_list[[i]] <- data.frame(
        x = draw_start, y = y_baseline, xend = draw_end, yend = y_baseline
      )

      # Calculate local coordinates for the visible portion
      local_start <- max(0, zoom_xlim[1] - contig_info$start)
      local_end <- min(contig_length, zoom_xlim[2] - contig_info$start)

      # Generate ticks based on the visible region in local coordinates
      # Use more ticks for small regions to ensure visibility
      n_ticks <- ifelse(local_end - local_start < 1000, 10, 5)
      tick_positions_local <- pretty(c(local_start, local_end), n = n_ticks)
      tick_positions_local <- tick_positions_local[tick_positions_local >= local_start &
        tick_positions_local <= local_end]

      # Convert back to global coordinates
      tick_positions_global <- contig_info$start + tick_positions_local

      # Filter ticks to be within the draw_start and draw_end
      visible_ticks_idx <- tick_positions_global >= draw_start & tick_positions_global <= draw_end
      tick_positions_global_visible <- tick_positions_global[visible_ticks_idx]
      tick_positions_local_visible <- tick_positions_local[visible_ticks_idx]

      if (length(tick_positions_global_visible) > 0) {
        # Determine a sensible tick height, e.g., based on a small fraction of the y-range or a fixed value
        # Since y-range is not explicitly defined for axis,
        # let's use a small absolute value or relative to x-axis range if dynamic
        # For simplicity, using a small fixed proportion of the x-axis zoom range. This might need adjustment.
        tick_height <- (zoom_xlim[2] - zoom_xlim[1]) * 0.0025
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

      # Define a suitable y_offset for contig names to bring them closer to the axis
      # Base the calculation on a unit equivalent to tick_height logic
      base_y_multiplier <- (zoom_xlim[2] - zoom_xlim[1]) * 0.005
      if (base_y_multiplier == 0) base_y_multiplier <- 0.1 # Ensure a minimum step, similar to tick_height

      if (length(tick_positions_global_visible) > 0) {
        # Ticks are present. Tick labels are at an offset of base_y_multiplier * 2.
        # Place contig names below tick labels, but closer to the axis.
        contig_name_y_offset <- base_y_multiplier * 2.5
      } else {
        # No ticks. Place contig names closer to the axis.
        contig_name_y_offset <- base_y_multiplier * 1.0
      }

      contig_names_list[[length(contig_names_list) + 1]] <- data.frame(
        x = (draw_start + draw_end) / 2,
        y = y_baseline - contig_name_y_offset, # Position below ticks
        label = contig_name,
        stringsAsFactors = FALSE
      )
    }

    # Combine list of data.frames into single data.frames
    # Safely rbind only if lists are not empty
    final_axis_lines_df <- if (length(axis_lines_list) > 0) do.call(rbind, axis_lines_list) else NULL
    final_tick_marks_df <- if (length(tick_marks_list) > 0) do.call(rbind, tick_marks_list) else NULL
    final_tick_labels_df <- if (length(tick_labels_list) > 0) do.call(rbind, tick_labels_list) else NULL
    final_contig_names_df <- if (length(contig_names_list) > 0) do.call(rbind, contig_names_list) else NULL

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
    # Add contig names to the plot
    if (!is.null(final_contig_names_df) && nrow(final_contig_names_df) > 0) {
      # Calculate the minimum y value to determine proper ylim
      min_y_value <- min(final_contig_names_df$y, na.rm = TRUE)
      # Add padding to ensure text isn't clipped
      padding_factor <- 1.5

      gg <- gg + ggplot2::geom_text(
        data = final_contig_names_df,
        ggplot2::aes(x = x, y = y, label = label),
        color = "black", size = 3.5, vjust = 1, fontface = "bold" # Adjust appearance as needed
      ) +
        # Set ylim to ensure enough space for labels
        ggplot2::ylim(min_y_value * padding_factor, y_baseline) +
        # Disable clipping to allow text to extend beyond plot boundaries
        ggplot2::coord_cartesian(clip = "off")
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
