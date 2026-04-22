
default_line_profile_params <- function(id, default_height) {
  list(
    height = list(
      group_id = id,
      type     = "integer",
      default  = default_height
    ),
    bin_mode = list(
      group_id = id,
      type     = "select",
      default  = "auto",
      choices  = c("auto", "full", "10", "100", "1000")
    ),
    max_points = list(
      group_id = id,
      type     = "integer",
      default  = 400
    ),
    agg_method = list(
      group_id = id,
      type     = "select",
      default  = "mean",
      choices  = c("mean", "median")
    )
  )
}

line_profile <- function(id, name,
                         points_f,
                         height        = 100,
                         color         = "#2171B5",
                         auto_register = TRUE) {

  BIN_LADDER <- c(1L, 10L, 100L, 1000L, 10000L)

  snap_bin_size <- function(raw) {
    candidates <- BIN_LADDER[BIN_LADDER <= raw]
    if (length(candidates) == 0L) return(BIN_LADDER[1L])
    max(candidates)
  }

  agg_fun <- function(method) {
    if (method == "median") stats::median else mean
  }

  plot_f <- function(profile, gg) {
    assembly   <- cxt_get_assembly()
    xlim       <- cxt_get_xlim()
    bin_mode   <- if (!is.null(profile$bin_mode))   profile$bin_mode   else "auto"
    max_points <- if (!is.null(profile$max_points)) as.integer(profile$max_points) else 400L
    method     <- if (!is.null(profile$agg_method)) profile$agg_method else "mean"
    agg        <- agg_fun(method)

    pts <- if (is.function(profile$points_f)) profile$points_f(assembly) else NULL
    if (is.null(pts) || nrow(pts) == 0)
      return(list(plot = gg, legends = list()))

    filtered <- cxt_filter_coords(pts)
    if (is.null(filtered) || nrow(filtered) == 0)
      return(list(plot = gg, legends = list()))

    view_width <- diff(xlim)

    if (bin_mode == "full") {
      # one point per site, no binning
      plot_df <- data.frame(
        x     = filtered$vcoord,
        y     = filtered$value,
        label = paste0(filtered$contig, ":", filtered$coord, "\n", method, "=", round(filtered$value, 3)),
        stringsAsFactors = FALSE
      )
      gg <- gg +
        ggplot2::geom_point(
          data = plot_df,
          ggplot2::aes(x = x, y = y, text = label),
          color = profile$color, size = 1) +
        ggplot2::geom_line(
          data = plot_df,
          ggplot2::aes(x = x, y = y),
          color = profile$color, linewidth = 0.4)
    } else {
      bin_size <- if (bin_mode == "auto") {
        snap_bin_size(ceiling(view_width / max_points))
      } else {
        as.integer(bin_mode)
      }

      # assign each point to a bin by its genomic coord (on the original contig)
      filtered$bin_start <- floor(filtered$coord / bin_size) * bin_size

      # aggregate per (contig, bin_start) — vectorized via tapply
      keys    <- paste(filtered$contig, filtered$bin_start, sep = "\t")
      agg_y   <- tapply(filtered$value, keys, agg)
      agg_n   <- tapply(filtered$value, keys, length)
      ukeys   <- names(agg_y)
      tab_pos <- regexpr("\t", ukeys, fixed = TRUE)
      bin_df  <- data.frame(
        contig    = substr(ukeys, 1L, tab_pos - 1L),
        bin_start = as.integer(substr(ukeys, tab_pos + 1L, nchar(ukeys))),
        y         = as.numeric(agg_y),
        n         = as.integer(agg_n),
        bin_end   = as.integer(substr(ukeys, tab_pos + 1L, nchar(ukeys))) + bin_size - 1L,
        stringsAsFactors = FALSE
      )

      # map bin midpoints to virtual coords; build a single data.frame to preserve types
      mid_input <- data.frame(
        contig    = bin_df$contig,
        coord     = bin_df$bin_start + as.integer(bin_size / 2L),
        value     = as.numeric(bin_df$y),
        bin_start = as.integer(bin_df$bin_start),
        bin_end   = as.integer(bin_df$bin_end),
        n         = as.integer(bin_df$n),
        stringsAsFactors = FALSE
      )
      mid_mapped <- cxt_filter_coords(mid_input)
      if (is.null(mid_mapped) || nrow(mid_mapped) == 0)
        return(list(plot = gg, legends = list()))

      mid_mapped$value     <- as.numeric(mid_mapped$value)
      mid_mapped$bin_start <- as.integer(mid_mapped$bin_start)
      mid_mapped$bin_end   <- as.integer(mid_mapped$bin_end)
      mid_mapped$n         <- as.integer(mid_mapped$n)
      mid_mapped <- mid_mapped[order(mid_mapped$vcoord), ]

      plot_df <- data.frame(
        x     = mid_mapped$vcoord,
        y     = mid_mapped$value,
        label = paste0(mid_mapped$contig, ":", mid_mapped$bin_start, "-", mid_mapped$bin_end,
                       "\n", method, "=", round(mid_mapped$value, 3),
                       "  n=", mid_mapped$n),
        stringsAsFactors = FALSE
      )

      gg <- gg +
        ggplot2::geom_line(
          data = plot_df,
          ggplot2::aes(x = x, y = y, text = label),
          color = profile$color, linewidth = 0.6) +
        ggplot2::geom_point(
          data = plot_df,
          ggplot2::aes(x = x, y = y, text = label),
          color = profile$color, size = 1)
    }

    gg <- gg + ggplot2::ylab(name)
    list(plot = gg, legends = list())
  }

  params <- default_line_profile_params(id, height)

  profile_create(
    id            = id,
    name          = name,
    type          = "line",
    height        = height,
    attr          = list(),
    params        = params,
    plot_f        = plot_f,
    points_f      = points_f,
    color         = color,
    auto_register = auto_register
  )
}
