default_differing_variants_params <- list(
  height = list(
    group_id = "differing_variants",
    type = "integer",
    default = 60
  )
)

# profile that shows variants differing between two selected samples.
# expects get_variants_f to return a data frame with columns:
#   contig, coord, variant_id, host_bin, type, is_genic, freq1, freq2
# freq1/freq2 must already be normalized to sample1/sample2 order.
# colors: green = gained in sample2 (freq2 > freq1), red = lost in sample2.
differing_variants_profile <- function(id, name, get_variants_f,
                                       height = 60, is_fixed = TRUE,
                                       params = default_differing_variants_params,
                                       auto_register = TRUE) {
  plot_f <- function(profile, gg) {
    df <- profile$get_variants_f()
    if (is.null(df) || nrow(df) == 0) {
      return(list(plot = gg, legends = list()))
    }

    # filter to visible view range and compute gcoord
    df <- cxt_filter_coords(df)
    if (is.null(df) || nrow(df) == 0) {
      return(list(plot = gg, legends = list()))
    }

    df$hover <- paste0(
      "Variant: ", df$variant_id, "<br>",
      "Bin: ", df$host_bin, "<br>",
      "Type: ", df$type, "<br>",
      "Genic: ", df$is_genic, "<br>",
      "freq1: ", round(df$freq1, 3), "<br>",
      "freq2: ", round(df$freq2, 3)
    )

    # split by direction: gained in sample2 (green) vs lost (red)
    gained <- df[df$freq2 > df$freq1, ]
    lost   <- df[df$freq2 <= df$freq1, ]

    if (nrow(gained) > 0) {
      gg <- gg + ggplot2::geom_segment(
        data = gained,
        ggplot2::aes(x = gcoord, xend = gcoord, y = -0.5, yend = 0.5, text = hover),
        color = "#29a329", linewidth = 1
      )
    }
    if (nrow(lost) > 0) {
      gg <- gg + ggplot2::geom_segment(
        data = lost,
        ggplot2::aes(x = gcoord, xend = gcoord, y = -0.5, yend = 0.5, text = hover),
        color = "#cc2200", linewidth = 1
      )
    }

    gg <- gg + ggplot2::ylim(-0.5, 0.5)

    xlim <- cxt_get_xlim()
    gg <- suppressWarnings(gg + ggplot2::coord_cartesian(xlim = xlim, clip = "off"))

    return(list(plot = gg, legends = list()))
  }

  profile_create(
    id = id, name = name, type = "differing_variants",
    height = height, is_fixed = is_fixed,
    attr = list(hide_y_ticks = TRUE),
    params = params, plot_f = plot_f,
    get_variants_f = get_variants_f,
    auto_register = auto_register
  )
}
