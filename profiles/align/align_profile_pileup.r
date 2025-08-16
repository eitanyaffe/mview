# Pileup mode for alignment profile

align_query_pileup_mode <- function(aln, cxt, report_mode_str = "all") {
  # Create cache key based on all relevant parameters
  # Use address of external pointer as unique identifier for alignment
  aln_id <- if (is(aln, "externalptr")) {
    paste0("ptr_", format(aln))
  } else {
    digest::digest(aln, algo = "md5")
  }
  
  cache_key <- paste0("pileup_query_",
                     digest::digest(list(
                       aln_id = aln_id,
                       xlim = cxt$mapper$xlim,
                       report_mode_str = report_mode_str
                     ), algo = "md5"))

  # Use cache for the pileup query
  df <- cache(cache_key, {
    aln_query_pileup(aln, cxt$intervals, report_mode_str = report_mode_str)
  })

  if (!is.null(df) && nrow(df) > 0) {
    # position from alntools is already 1-based, no need to add 1
    df$coord <- df$position
    return(filter_coords(df, cxt, cxt$mapper$xlim))
  }

  return(NULL)
}

align_profile_pileup <- function(profile, cxt, aln, gg) {
  df <- align_query_pileup_mode(aln, cxt, report_mode_str = "all")
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }
  
  cat(sprintf("plotting pileup with %d variant records\n", nrow(df)))
  
  # ensure gcoord exists (should be added by filter_coords in align_query_pileup_mode)
  if (!"gcoord" %in% colnames(df)) {
    warning("gcoord column missing from pileup data")
    return(gg)
  }
  
  # calculate cumulative heights using ave function
  df$y_bottom <- ave(df$count, df$coord, FUN = function(x) cumsum(c(0, x[-length(x)])))
  df$y_top <- df$y_bottom + df$count
  
  # assign colors based on variant type
  # use different colors for each variant to distinguish them in pileup view
  unique_variants <- unique(df$variant)
  n_variants <- length(unique_variants)
  
  if (n_variants <= 8) {
    # use predefined distinct colors for small numbers of variants
    distinct_colors <- c("#4e79a7", "#f28e2c", "#e15759", "#76b7b2", 
                        "#59a14f", "#edc949", "#af7aa1", "#ff9d9a")
    variant_colors <- setNames(distinct_colors[seq_len(n_variants)], unique_variants)
  } else {
    # use rainbow colors for larger numbers
    variant_colors <- setNames(rainbow(n_variants), unique_variants)
  }
  
  # special color for REF if present
  if ("REF" %in% unique_variants) {
    variant_colors["REF"] <- "#d3d3d3"
  }
  
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
