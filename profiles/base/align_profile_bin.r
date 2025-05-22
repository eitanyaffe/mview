# Bin mode for alignment profile

get_current_bin_size <- function(xlim, bin_type, target_bins = 1024) {
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

  # Calculate raw bin size to get close to target_bins
  raw_bin_size <- range_bp / target_bins

  # Ensure minimum bin size of 1
  return(max(1, round(raw_bin_size)))
}

align_query_bin_mode <- function(aln, cxt, bin_type, target_bins = 1024) {
  bin_size <- get_current_bin_size(cxt$mapper$xlim, bin_type, target_bins = target_bins)

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
  df <- align_query_bin_mode(aln, cxt, profile$bin_type, target_bins = profile$target_bins)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }
  # Calculate metrics
  df$cov <- ifelse(df$length > 0, df$read_count / df$length, 0)
  df$mut_density <- ifelse(df$read_count > 0, df$mutation_count / df$read_count, 0)

  # Calculate colors based on mutation density
  # Gray for 0 mutations per 100bp, red for 1 mutation per 100bp
  mutations_per_100bp <- df$mut_density * 100
  max_mutations_per_100bp <- 1.0

  # Use the general color scale function
  df$fill_color <- get_color_scale(
    values = mutations_per_100bp,
    colors = c("lightgray", "red"),
    min_val = 0,
    max_val = max_mutations_per_100bp,
    num_steps = 20
  )

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
    ), size = 0.2, color = NA
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_minimal()

  return(gg)
}
