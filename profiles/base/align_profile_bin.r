# Bin mode for alignment profile

get_current_bin_size <- function(xlim, bin_type) {
  if (bin_type != "auto") {
    bs <- as.numeric(bin_type)
    return(max(1, bs, na.rm = TRUE))
  }

  if (is.null(xlim) || length(xlim) != 2) {
    return(1000)
  }

  range_bp <- (xlim[2] + 1) - xlim[1]
  if (range_bp <= 0) {
    return(100)
  }

  target_bins <- 512
  raw_bin_size <- range_bp / target_bins

  if (raw_bin_size <= 10) {
    return(10)
  }
  if (raw_bin_size <= 100) {
    return(100)
  }
  if (raw_bin_size <= 1000) {
    return(1000)
  }
  if (raw_bin_size <= 10000) {
    return(10000)
  }
  return(100000)
}

align_query_bin_mode <- function(aln, cxt, bin_type) {
  bin_size <- get_current_bin_size(cxt$mapper$xlim, bin_type)

  df <- aln_query_bin(aln, cxt$intervals, bin_size)
  cat(sprintf("bin_size: %d, bin_count: %d\n", bin_size, nrow(df)))
  if (!is.null(df) && nrow(df) > 0) {
    df$start <- df$start + 1
    df$end <- df$end
    return(filter_segments(df, cxt, cxt$mapper$xlim))
  }

  return(NULL)
}

align_profile_bin <- function(profile, cxt, aln, gg) {
  df <- align_query_bin_mode(aln, cxt, profile$bin_type)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }
  # Calculate metrics
  df$cov <- ifelse(df$length > 0, df$read_count / df$length, 0)
  df$mut_density <- ifelse(df$length > 0, df$mutation_count / df$length, 0)

  # Calculate colors based on mutation density
  max_density <- 0.01
  intensity <- pmin(df$mut_density / max_density, 1.0)
  intensity[is.na(intensity) | !is.finite(intensity)] <- 0

  gray_rgb <- grDevices::col2rgb("lightgray") / 255
  red_rgb <- grDevices::col2rgb("red") / 255

  r_val <- gray_rgb[1, ] * (1 - intensity) + red_rgb[1, ] * intensity
  g_val <- gray_rgb[2, ] * (1 - intensity) + red_rgb[2, ] * intensity
  b_val <- gray_rgb[3, ] * (1 - intensity) + red_rgb[3, ] * intensity

  df$fill_color <- grDevices::rgb(r_val, g_val, b_val)

  # Hover text
  df$hover_text <- paste0(
    "Bin Contig: ", df$contig, "\n",
    "Local (1-based): ", df$start, "-", df$end, "\n",
    "Global (0-based): ", df$gstart, "-", df$gend, "\n",
    "Length: ", format_bp(df$length), "\n",
    "Reads: ", df$read_count, " (Cov: ", sprintf("%.2f", df$cov), "x)\n",
    "Muts: ", df$mutation_count, " (Dens: ", sprintf("%.4f", df$mut_density), " muts/bp)"
  )

  # Plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, text = hover_text
    ), size = 0.2
  )

  return(gg)
}
