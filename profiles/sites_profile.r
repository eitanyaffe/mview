# sites_profile: display CME predicted sites from SITES_PREDICT_ANNOTATE_SITES
#
# data source: per-assembly annotated sites file (real only), loaded lazily per assembly
# rendering:
#   zoomed out (view width > zoom_threshold): vertical segment at site midpoint
#   zoomed in  (view width <= zoom_threshold): filled rectangle spanning start..end
# color: per-quantile of score within the current assembly × CME
#        (stable across zoom/pan), or fixed gray if use_score_color = FALSE
#
# constructor parameters:
#   id, name      - profile identity
#   height        - track height (default 40)
#   sites_f       - function(assembly) returning the full sites data.frame or NULL
#   rv_cme_id     - reactive value returning the selected cme_id string
#   zoom_threshold- bp width below which to switch to rect rendering (default 400)
#   use_score_color - logical; if FALSE all sites are drawn in a fixed gray (default TRUE)
#   color_low / color_high - colour ramp endpoints (defaults: "#DEEBF7" / "#08306B")
#   auto_register - logical (default TRUE)

sites_score_quantile <- function(scores) {
  scores <- as.numeric(scores)
  scores[is.na(scores)] <- min(scores, na.rm = TRUE)
  (rank(scores, ties.method = "average") - 1) / max(1, length(scores) - 1)
}

sites_quantile_colors <- function(q, color_low, color_high) {
  ramp <- grDevices::colorRamp(c(color_low, color_high))
  rgb_m <- ramp(q)
  grDevices::rgb(rgb_m[, 1], rgb_m[, 2], rgb_m[, 3], maxColorValue = 255)
}

sites_build_hover <- function(df) {
  contig <- if ("contig" %in% names(df)) df$contig else rep("", nrow(df))
  start  <- if ("start" %in% names(df)) df$start else rep(NA, nrow(df))
  end    <- if ("end" %in% names(df)) df$end else rep(NA, nrow(df))
  strand <- if ("strand" %in% names(df)) df$strand else rep("", nrow(df))
  score  <- if ("score" %in% names(df)) df$score else rep(NA, nrow(df))
  pval   <- if ("p_value" %in% names(df)) df$p_value else rep(NA, nrow(df))
  host   <- if ("host_bin" %in% names(df)) df$host_bin else rep("", nrow(df))
  cls    <- if ("site_class" %in% names(df)) df$site_class else rep("", nrow(df))

  paste0(
    contig, ":", start, "-", end,
    "\nstrand=", strand,
    "\nscore=", score,
    "\np=", pval,
    "\nhost_bin=", host,
    "\nclass=", cls
  )
}

sites_profile <- function(id, name, height = 40, is_fixed = TRUE,
                          sites_f = NULL,
                          rv_cme_id = NULL,
                          zoom_threshold = 400,
                          use_score_color = TRUE,
                          color_low  = "#DEEBF7",
                          color_high = "#08306B",
                          color_fixed = "#808080",
                          auto_register = TRUE)
{
  plot_f <- function(profile, gg) {
    assembly <- cxt_get_assembly()

    sites_all <- NULL
    if (is.function(profile$sites_f)) {
      sites_all <- profile$sites_f(assembly)
    }

    if (is.null(sites_all) || nrow(sites_all) == 0) {
      return(list(plot = gg, legends = list()))
    }

    # filter to current assembly
    if ("assembly_id" %in% names(sites_all)) {
      sites_all <- sites_all[sites_all$assembly_id == assembly, ]
    }

    # filter to selected CME, only real sites
    cme_id <- if (is.function(profile$rv_cme_id)) profile$rv_cme_id() else NULL
    if (!is.null(cme_id) && nchar(cme_id) > 0 && "cme_id" %in% names(sites_all)) {
      sites_all <- sites_all[sites_all$cme_id == cme_id, ]
    }
    if ("site_class" %in% names(sites_all)) {
      sites_all <- sites_all[sites_all$site_class != "control", ]
    }

    if (is.null(sites_all) || nrow(sites_all) == 0) {
      return(list(plot = gg, legends = list()))
    }

    sites_all$start <- as.integer(sites_all$start)
    sites_all$end   <- as.integer(sites_all$end)

    # quantile over full per-CME set (stable colors across zoom/pan)
    if (profile$use_score_color && "score" %in% names(sites_all)) {
      sites_all$q <- sites_score_quantile(sites_all$score)
    }

    filtered <- cxt_filter_intervals(sites_all, merge_adjacent = TRUE)

    if (is.null(filtered) || nrow(filtered) == 0) {
      return(list(plot = gg, legends = list()))
    }

    if (profile$use_score_color && "q" %in% names(filtered)) {
      fill_colors <- sites_quantile_colors(filtered$q, profile$color_low, profile$color_high)
    } else {
      fill_colors <- rep(profile$color_fixed, nrow(filtered))
    }

    filtered$hover <- sites_build_hover(filtered)

    xlim <- cxt_get_xlim()
    view_width <- if (length(xlim) == 2) diff(as.numeric(xlim)) else Inf

    if (view_width > profile$zoom_threshold) {
      # zoomed out: vertical segment at midpoint
      filtered$vmid <- (filtered$vstart + filtered$vend) / 2
      gg <- gg +
        ggplot2::geom_segment(
          data = filtered,
          ggplot2::aes(x = vmid, xend = vmid, y = -0.45, yend = 0.45, text = hover),
          color = fill_colors,
          linewidth = 0.6
        )
    } else {
      # zoomed in: filled rectangle spanning the full motif
      # vstart/vend are inclusive contig coords; nt n is drawn centered on x=n,
      # so its visual cell is [n-0.5, n+0.5]
      filtered$xmin_r <- filtered$vstart - 0.5
      filtered$xmax_r <- filtered$vend   + 0.5
      gg <- gg +
        ggplot2::geom_rect(
          data = filtered,
          ggplot2::aes(xmin = xmin_r, xmax = xmax_r, ymin = -0.35, ymax = 0.35, text = hover),
          fill  = fill_colors,
          color = "black",
          linewidth = 0.3
        )
    }

    gg <- gg + ggplot2::ylim(-0.5, 0.5)

    return(list(plot = gg, legends = list()))
  }

  profile_create(
    id             = id,
    name           = name,
    type           = "sites",
    height         = height,
    is_fixed       = is_fixed,
    attr           = list(hide_y_ticks = TRUE),
    params         = NULL,
    plot_f         = plot_f,
    sites_f        = sites_f,
    rv_cme_id      = rv_cme_id,
    zoom_threshold = zoom_threshold,
    use_score_color = use_score_color,
    color_low      = color_low,
    color_high     = color_high,
    color_fixed    = color_fixed,
    auto_register  = auto_register
  )
}
