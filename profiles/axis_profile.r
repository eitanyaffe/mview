#' Create a simple axis profile
#' @param id Unique profile id
#' @param name Display name
#' @param height Height in pixels
#' @param auto_register Whether to register immediately
#' @return Profile object
axis_profile <- function(id = "simple_axis",
                         name = "simple axis",
                         height = 100,
                         auto_register = TRUE) {
  data_f <- function(cxt) NULL

  # compute coordinate annotations for the current view
  get_annotations_f <- function(profile, cxt) {
    if (is.null(cxt$mapper) || is.null(cxt$mapper$cdf) || nrow(cxt$mapper$cdf) == 0) return(NULL)

    contig_df_all <- cxt$mapper$cdf
    zoom_xlim <- if (!is.null(cxt$mapper$xlim)) range(cxt$mapper$xlim) else c(min(contig_df_all$start), max(contig_df_all$end))
    contig_df <- contig_df_all[contig_df_all$start < zoom_xlim[2] & contig_df_all$end > zoom_xlim[1], , drop = FALSE]

    if (nrow(contig_df) != 1) return(NULL)

    contig_start <- contig_df$start[1]
    contig_end <- contig_df$end[1]

    # compute tick marks in contig-local coordinates
    local_visible_start <- max(0, zoom_xlim[1] - contig_start)
    local_visible_end <- min(contig_end - contig_start, zoom_xlim[2] - contig_start)
    local_visible_range <- local_visible_end - local_visible_start
    if (local_visible_range <= 0) return(NULL)

    target_ticks <- 10
    raw_interval <- local_visible_range / target_ticks
    magnitude <- 10^floor(log10(raw_interval))
    normalized <- raw_interval / magnitude
    nice_interval <- if (normalized <= 1) {
        1 * magnitude
    } else if (normalized <= 2) {
        2 * magnitude
    } else if (normalized <= 5) {
        5 * magnitude
    } else {
        10 * magnitude
    }
    nice_interval <- max(1, nice_interval)

    start_tick_local <- ceiling(local_visible_start / nice_interval) * nice_interval
    end_tick_local <- floor(local_visible_end / nice_interval) * nice_interval
    local_ticks <- if (start_tick_local <= end_tick_local) {
        seq(start_tick_local, end_tick_local, by = nice_interval)
    } else {
        numeric(0)
    }
    if (length(local_ticks) > 0) {
        local_ticks <- round(local_ticks)
        local_ticks <- unique(local_ticks)
        local_ticks <- local_ticks[local_ticks != 0]
    }
    if (length(local_ticks) == 0) return(NULL)

    global_ticks <- contig_start + local_ticks
    cat("axis_profile: returning ", length(global_ticks), " tick positions\n")
    return(data.frame(
        x = global_ticks,
        y = 0,
        label = format(local_ticks, big.mark = ",", scientific = FALSE)
    ))
  }

  plot_f <- function(profile, cxt, gg) {
    if (is.null(cxt$mapper) || is.null(cxt$mapper$cdf) || nrow(cxt$mapper$cdf) == 0) return(gg)

    contig_df_all <- cxt$mapper$cdf
    zoom_xlim <- if (!is.null(cxt$mapper$xlim)) range(cxt$mapper$xlim) else c(min(contig_df_all$start), max(contig_df_all$end))
    contig_df <- contig_df_all[contig_df_all$start < zoom_xlim[2] & contig_df_all$end > zoom_xlim[1], , drop = FALSE]

    if (nrow(contig_df) == 0) return(gg)

    if (nrow(contig_df) == 1) {
      name_data <- data.frame(x = (zoom_xlim[1] + zoom_xlim[2]) / 2, y = 0.45, label = contig_df$contig)
    } else {
      name_data <- data.frame(x = (contig_df$start + contig_df$end) / 2, y = 0.55, label = contig_df$contig)
    }
    gg <- gg + ggplot2::geom_text(data = name_data, ggplot2::aes(x = x, y = y, label = label), color = "black", size = 3, vjust = 0.5)

    if (nrow(contig_df) == 1) {
      axis_start <- max(min(contig_df$start), zoom_xlim[1])
      axis_end <- min(max(contig_df$end), zoom_xlim[2])
      axis_line_data <- data.frame(x = axis_start, y = 0.8, xend = axis_end, yend = 0.8)
      gg <- gg + ggplot2::geom_segment(data = axis_line_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 0.4)
    }

    if (nrow(contig_df) == 1) {
      contig_start <- contig_df$start[1]
      contig_end <- contig_df$end[1]

      # compute tick marks in contig-local coordinates
      local_visible_start <- max(0, zoom_xlim[1] - contig_start)
      local_visible_end <- min(contig_end - contig_start, zoom_xlim[2] - contig_start)
      local_visible_range <- local_visible_end - local_visible_start
      if (local_visible_range > 0) {
        target_ticks <- 10
        raw_interval <- local_visible_range / target_ticks
        magnitude <- 10^floor(log10(raw_interval))
        normalized <- raw_interval / magnitude
        nice_interval <- if (normalized <= 1) 1 * magnitude else if (normalized <= 2) 2 * magnitude else if (normalized <= 5) 5 * magnitude else 10 * magnitude
        nice_interval <- max(1, nice_interval)

        start_tick_local <- ceiling(local_visible_start / nice_interval) * nice_interval
        end_tick_local <- floor(local_visible_end / nice_interval) * nice_interval
        local_ticks <- if (start_tick_local <= end_tick_local) seq(start_tick_local, end_tick_local, by = nice_interval) else numeric(0)
        if (length(local_ticks) > 0) {
          local_ticks <- round(local_ticks)
          local_ticks <- unique(local_ticks)
          local_ticks <- local_ticks[local_ticks != 0]
        }
        if (length(local_ticks) > 0) {
          tick_height <- 0.05
          global_ticks <- contig_start + local_ticks
          tick_data <- data.frame(x = global_ticks, y = 0.8, xend = global_ticks, yend = 0.8 - tick_height)
          gg <- gg + ggplot2::geom_segment(data = tick_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 0.3)
        }
      }
    }

    gg + ggplot2::coord_cartesian(ylim = c(0.35, 0.9), clip = "off")
  }

  profile_create(
    id = id,
    name = name,
    type = "axis",
    height = height,
    attr = list(hide_y_label = TRUE, hide_y_ticks = TRUE),
    data_f = data_f,
    plot_f = plot_f,
    get_annotations = get_annotations_f,
    auto_register = auto_register
  )
}