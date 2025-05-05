# profiles.r

register_profile("randbar", list(
  id = "randbar",
  name = "Random Bar Height",
  plot_f = function(data, mapper) {
    data$value <- runif(nrow(data))
    ggplot(data) +
      geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value), fill = "skyblue", alpha = 0.6) +
      theme_minimal()
  },
  info_f = function(cid, coord, y) {
    paste0("randbar: ", cid, ":", coord, " = ", round(y, 2))
  },
  size_f = function() 1,
  spec_f = function() list(ylabel = "Random Value", grid = TRUE)
))

register_profile("flatline", list(
  id = "flatline",
  name = "Flatline",
  plot_f = function(data, mapper) {
    data$value <- 0.5
    ggplot(data) +
      geom_segment(aes(x = start, xend = end, y = value, yend = value),
        color = "orange", linewidth = 1
      ) +
      theme_minimal()
  },
  info_f = function(cid, coord, y) {
    paste0("flatline: ", cid, ":", coord, " = ", round(y, 2))
  },
  size_f = function() 0.5,
  spec_f = function() list(ylabel = "Flatline", grid = FALSE)
))

register_profile("line", list(
  id = "line",
  name = "Line Profile",
  plot_f = function(data, mapper) {
    set.seed(123)
    profile_df <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) {
      row <- data[i, ]
      x <- seq(row$start, row$end, length.out = 20)
      y <- runif(length(x))
      data.frame(global = x, value = y, cid = row$cid)
    }))

    ggplot(profile_df, aes(x = global, y = value, group = cid)) +
      geom_line(color = "steelblue") +
      theme_minimal()
  },
  info_f = function(cid, coord, y) {
    paste0("line: ", cid, ":", coord, " = ", round(y, 2))
  },
  size_f = function() 0.7,
  spec_f = function() list(ylabel = "Line Profile", grid = TRUE)
))

register_profile("coords", list(
  id = "coords",
  name = "Coordinates",
  plot_f = function(data, mapper) {
    # Create base plot
    p <- ggplot() +
      theme_minimal()

    # Add global coordinate ruler on top
    p <- p + geom_segment(
      data = data.frame(x = 0, xend = max(data$end)),
      aes(x = x, xend = xend, y = 0.75, yend = 0.75),
      color = "black", linewidth = 1
    )

    # Add global ticks and labels
    ticks <- seq(0, max(data$end), length.out = min(10, max(data$end) / 100 + 1))
    p <- p + geom_segment(
      data = data.frame(x = ticks),
      aes(x = x, xend = x, y = 0.7, yend = 0.8),
      color = "black"
    )
    p <- p + geom_text(
      data = data.frame(x = ticks, label = format(round(ticks), big.mark = ",")),
      aes(x = x, y = 0.9, label = label),
      size = 3
    )

    # Add per-contig rulers
    for (i in 1:nrow(data)) {
      row <- data[i, ]
      p <- p + geom_segment(
        data = data.frame(x = row$start, xend = row$end),
        aes(x = x, xend = xend, y = 0.4, yend = 0.4),
        color = "blue", linewidth = 1
      )

      # Contig ticks and labels
      contig_ticks <- seq(1, row$length, length.out = min(5, row$length / 100 + 1))
      p <- p + geom_segment(
        data = data.frame(
          x = mapper$l2g(row$cid, contig_ticks)
        ),
        aes(x = x, xend = x, y = 0.35, yend = 0.45),
        color = "blue"
      )

      # Contig ID and tick labels
      mid_pos <- row$start + (row$end - row$start) / 2
      p <- p + geom_text(
        data = data.frame(x = mid_pos, label = row$cid),
        aes(x = x, y = 0.25, label = label),
        size = 3, color = "blue"
      )

      p <- p + geom_text(
        data = data.frame(
          x = mapper$l2g(row$cid, contig_ticks),
          label = format(round(contig_ticks), big.mark = ",")
        ),
        aes(x = x, y = 0.15, label = label),
        size = 2.5, color = "blue"
      )
    }

    # Set y limits and remove axis elements
    p <- p + ylim(0, 1) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      )

    return(p)
  },
  info_f = function(cid, coord, y) {
    paste0(
      "Position: ", cid, ":", format(coord, big.mark = ",")
    )
  },
  size_f = function() 1.5, # Fixed size factor for 150px (relative to other profiles)
  spec_f = function() list(ylabel = "Coordinates", grid = FALSE)
))
