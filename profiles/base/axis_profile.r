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
  plot_f <- function(cxt, gg) {
    if (is.null(cxt$mapper) || is.null(cxt$mapper$contigs)) {
      return(list(plot = gg, info_df = NULL))
    }

    contigs <- cxt$mapper$contigs
    y_baseline <- 0
    info_df <- NULL

    for (i in seq_along(contigs)) {
      contig <- contigs[[i]]
      name <- names(contigs)[i]

      # Get contig coordinates
      global_start <- contig$start
      global_end <- contig$end
      contig_length <- contig$length

      # Draw main axis line
      gg <- gg + ggplot2::geom_segment(
        data = data.frame(x = global_start, y = y_baseline, xend = global_end, yend = y_baseline),
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "black",
        size = 0.5
      )

      # Generate ticks for local coordinates
      tick_positions_local <- pretty(c(0, contig_length), n = 5)
      tick_positions_local <- tick_positions_local[tick_positions_local >= 0 & tick_positions_local <= contig_length]
      tick_positions_global <- global_start + tick_positions_local

      if (length(tick_positions_global) > 0) {
        # Draw tick marks
        gg <- gg + ggplot2::geom_segment(
          data = data.frame(
            x = tick_positions_global,
            y = y_baseline,
            xend = tick_positions_global,
            yend = y_baseline - (cxt$mapper$y_max - cxt$mapper$y_min) * 0.02
          ),
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          color = "black",
          size = 0.5
        )

        # Add tick labels (local coordinates)
        gg <- gg + ggplot2::geom_text(
          data = data.frame(
            x = tick_positions_global,
            y = y_baseline - (cxt$mapper$y_max - cxt$mapper$y_min) * 0.04,
            label = tick_positions_local
          ),
          ggplot2::aes(x = x, y = y, label = label),
          color = "black",
          size = 3,
          vjust = 1
        )

        # Add to info_df for hover
        if (is.null(info_df)) {
          info_df <- data.frame(
            gcoord = tick_positions_global,
            y_value = y_baseline,
            description = paste0("Local coordinate: ", tick_positions_local),
            raw_y = y_baseline,
            stringsAsFactors = FALSE
          )
        } else {
          info_df <- rbind(info_df, data.frame(
            gcoord = tick_positions_global,
            y_value = y_baseline,
            description = paste0("Local coordinate: ", tick_positions_local),
            raw_y = y_baseline,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    list(plot = gg, info_df = info_df)
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
