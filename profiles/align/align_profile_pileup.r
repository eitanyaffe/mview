# Pileup mode for alignment profile

align_query_pileup_mode <- function(aln, cxt) {
  # Query pileup data
  df <- aln_query_pileup(aln, cxt$intervals, report_mode_str = "all")

  if (!is.null(df) && nrow(df) > 0) {
    # position from alntools is already 1-based, no need to add 1
    df$coord <- df$position
    return(filter_coords(df, cxt, cxt$mapper$xlim))
  }

  return(NULL)
}

align_profile_pileup <- function(profile, cxt, aln, gg) {
  df <- align_query_pileup_mode(aln, cxt)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }
  
  cat(sprintf("plotting pileup with %d variant records\n", nrow(df)))
  
  # ensure gcoord exists (should be added by filter_coords in align_query_pileup_mode)
  if (!"gcoord" %in% colnames(df)) {
    warning("gcoord column missing from pileup data")
    return(gg)
  }
  
  # get shared variant colors (consistent between full and pileup modes)
  variant_colors <- get_variant_colors(aln, cxt)
  
  # calculate cumulative heights using ave function
  df$y_bottom <- ave(df$count, df$coord, FUN = function(x) cumsum(c(0, x[-length(x)])))
  df$y_top <- df$y_bottom + df$count
  
  # add color and position info  
  df$fill_color <- variant_colors[df$variant]
  # position rectangles centered on the coordinate, similar to mutation display in full mode
  df$x_left <- df$gcoord - 0.5
  df$x_right <- df$gcoord + 0.5
  
  # create hover text
  df$hover_text <- paste0(
    "Position: ", df$coord, "\n",
    "Variant: ", df$variant, "\n",
    "Count: ", df$count, " / ", df$coverage, "\n",
    "Fraction: ", sprintf("%.3f", df$count / df$coverage)
  )
  
  # plot rectangles
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = x_left, xmax = x_right,
      ymin = y_bottom, ymax = y_top,
      fill = fill_color,
      text = hover_text
    ),
    color = NA
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
  
  return(gg)
}
