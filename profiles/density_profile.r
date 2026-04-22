
################################################################################
# SHARED UTILITIES
################################################################################

BIN_LADDER <- c(10L, 20L, 50L, 100L, 1000L, 10000L)

snap_bin_size <- function(raw) {
  candidates <- BIN_LADDER[BIN_LADDER <= raw]
  if (length(candidates) == 0L) return(BIN_LADDER[1L])
  max(candidates)
}

density_compute_bins <- function(profile, assembly, xlim) {
  pts <- if (is.function(profile$points_f)) profile$points_f(assembly) else NULL
  if (is.null(pts) || nrow(pts) == 0)
    return(NULL)

  suppressWarnings(g <- cxt_contig2global(pts$contig, pts$coord))
  pts$gcoord <- g
  pts <- pts[!is.na(pts$gcoord), , drop = FALSE]
  if (nrow(pts) == 0)
    return(NULL)

  bin_mode    <- if (!is.null(profile$bin_mode))    profile$bin_mode              else "auto"
  max_points  <- if (!is.null(profile$max_points)) as.integer(profile$max_points) else 400L

  view_width <- diff(xlim)
  bin_size   <- if (bin_mode == "auto") {
    snap_bin_size(ceiling(view_width / max_points))
  } else {
    as.integer(bin_mode)
  }

  grid_start <- floor(xlim[1] / bin_size) * bin_size - bin_size
  grid_end   <- floor(xlim[2] / bin_size) * bin_size + bin_size
  bin_starts <- seq(grid_start, grid_end, by = bin_size)

  pts$bin_start <- floor(pts$gcoord / bin_size) * bin_size
  counts        <- tapply(rep(1L, nrow(pts)), pts$bin_start, sum)
  n             <- rep(0L, length(bin_starts))
  idx           <- match(as.integer(names(counts)), bin_starts)
  valid         <- !is.na(idx)
  n[idx[valid]] <- as.integer(counts)[valid]

  density <- n / (bin_size / 1000)

  list(
    bin_starts = bin_starts,
    n          = n,
    density    = density,
    bin_size   = bin_size
  )
}

################################################################################
# LINE PLOT MODE
################################################################################

density_plot_line <- function(profile, gg, bins) {
  show_counts <- isTRUE(profile$show_counts)

  y_vals  <- if (show_counts) as.numeric(bins$n) else as.numeric(bins$density)
  y_label <- if (show_counts) "sites" else "sites / kb"
  hover_y <- if (show_counts)
    paste0(bins$n, " sites")
  else
    paste0(round(bins$density, 2), " sites/kb", "  (n=", bins$n, ")")

  plot_df <- data.frame(
    x     = bins$bin_starts + bins$bin_size / 2,
    y     = y_vals,
    label = paste0("bin ", bins$bin_starts, "-", bins$bin_starts + bins$bin_size - 1L, "\n", hover_y),
    stringsAsFactors = FALSE
  )

  gg <- gg +
    ggplot2::geom_line(
      data = plot_df,
      ggplot2::aes(x = x, y = y),
      color = profile$color, linewidth = 0.6) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = x, y = y, text = label),
      color = profile$color, size = 1) +
    ggplot2::ylab(y_label)

  list(plot = gg, legends = list())
}

################################################################################
# HEATMAP PLOT MODE
################################################################################

density_plot_heatmap <- function(profile, gg, bins) {
  show_counts <- isTRUE(profile$show_counts)

  y_vals  <- if (show_counts) as.numeric(bins$n) else as.numeric(bins$density)
  y_label <- if (show_counts) "sites" else "sites / kb"
  hover_y <- if (show_counts)
    paste0(bins$n, " sites")
  else
    paste0(round(bins$density, 2), " sites/kb", "  (n=", bins$n, ")")

  max_val <- max(y_vals, na.rm = TRUE)
  if (is.na(max_val) || max_val == 0) max_val <- 1

  norm_vals <- y_vals / max_val
  base_rgb  <- grDevices::col2rgb(profile$color)[, 1] / 255
  fill_colors <- grDevices::rgb(
    1 - norm_vals * (1 - base_rgb[1]),
    1 - norm_vals * (1 - base_rgb[2]),
    1 - norm_vals * (1 - base_rgb[3])
  )

  plot_df <- data.frame(
    xmin  = bins$bin_starts,
    xmax  = bins$bin_starts + bins$bin_size,
    ymin  = 0,
    ymax  = 1,
    fill  = fill_colors,
    label = paste0("bin ", bins$bin_starts, "-", bins$bin_starts + bins$bin_size - 1L, "\n", hover_y),
    stringsAsFactors = FALSE
  )

  gg <- gg +
    ggplot2::geom_rect(
      data = plot_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill), text = label)
    ) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::ylab(y_label)

  legend_plot <- density_heatmap_legend(profile$color, max_val, y_label)

  list(plot = gg, legends = list(list(plot = legend_plot, width = 100, height = 200)))
}

density_heatmap_legend <- function(color, max_val, y_label) {
  n_steps <- 100
  gradient_vals <- seq(0, max_val, length.out = n_steps)
  norm_vals <- gradient_vals / max_val
  base_rgb  <- grDevices::col2rgb(color)[, 1] / 255
  fill_colors <- grDevices::rgb(
    1 - norm_vals * (1 - base_rgb[1]),
    1 - norm_vals * (1 - base_rgb[2]),
    1 - norm_vals * (1 - base_rgb[3])
  )

  legend_df <- data.frame(
    ymin = gradient_vals[-length(gradient_vals)],
    ymax = gradient_vals[-1],
    fill = fill_colors[-length(fill_colors)],
    stringsAsFactors = FALSE
  )

  breaks_vals <- c(0, max_val / 2, max_val)
  breaks_labels <- sprintf("%.1f", breaks_vals)

  gg_legend <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = legend_df,
      ggplot2::aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = I(fill))
    ) +
    ggplot2::scale_y_continuous(breaks = breaks_vals, labels = breaks_labels) +
    ggplot2::labs(y = y_label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()
    )

  gg_legend
}

################################################################################
# PROFILE CONSTRUCTOR
################################################################################

density_profile <- function(id, name,
                            points_f,
                            height        = 100,
                            color         = "#2171B5",
                            param_group   = id,
                            auto_register = TRUE) {

  plot_f <- function(profile, gg) {
    assembly   <- cxt_get_assembly()
    xlim       <- cxt_get_xlim()
    plot_mode  <- if (!is.null(profile$plot_mode)) profile$plot_mode else "line"

    bins <- density_compute_bins(profile, assembly, xlim)
    if (is.null(bins))
      return(list(plot = gg, legends = list()))

    if (plot_mode == "heatmap") {
      density_plot_heatmap(profile, gg, bins)
    } else {
      density_plot_line(profile, gg, bins)
    }
  }

  params <- list(
    height = list(
      group_id = param_group,
      type     = "integer",
      default  = height
    ),
    plot_mode = list(
      group_id = param_group,
      type     = "select",
      default  = "line",
      choices  = c("line", "heatmap")
    ),
    bin_mode = list(
      group_id = param_group,
      type     = "select",
      default  = "auto",
      choices  = c("auto", "10", "20", "50", "100", "1000")
    ),
    max_points = list(
      group_id = param_group,
      type     = "integer",
      default  = 400
    ),
    show_counts = list(
      group_id = param_group,
      type     = "boolean",
      default  = FALSE
    )
  )

  profile_create(
    id            = id,
    name          = name,
    type          = "density",
    height        = height,
    attr          = list(),
    params        = params,
    plot_f        = plot_f,
    points_f      = points_f,
    color         = color,
    auto_register = auto_register
  )
}
