# Pileup mode for alignment profile

align_query_pileup_mode <- function(aln, report_mode_str = "all", clip_mode = "all", clip_margin = 10, min_mutations_percent = 0.0, max_mutations_percent = 10.0, min_alignment_length = 0, max_alignment_length = 0, min_indel_length = 3) {
  intervals <- cxt_get_zoom_view()
  plotted_segs <- cxt_get_plotted_segments()
  
  # Create cache key based on all relevant parameters
  # Use address of external pointer as unique identifier for alignment
  aln_id <- if (is(aln, "externalptr")) {
    paste0("ptr_", format(aln))
  } else {
    digest::digest(aln, algo = "md5")
  }
  
  # include plotted segments with strands for cache invalidation when segments flip
  seg_key <- if (nrow(plotted_segs) > 0) {
    paste(plotted_segs$segment, plotted_segs$strand, collapse = ",")
  } else {
    ""
  }
  
  cache_key <- paste0("pileup_query_",
                     digest::digest(list(
                       aln_id = aln_id,
                       seg_key = seg_key,
                       report_mode_str = report_mode_str,
                       clip_mode = clip_mode,
                       clip_margin = clip_margin,
                       min_mutations_percent = min_mutations_percent,
                       max_mutations_percent = max_mutations_percent,
                       min_alignment_length = min_alignment_length,
                       max_alignment_length = max_alignment_length,
                       min_indel_length = min_indel_length
                     ), algo = "md5"))

  # Use cache for the pileup query
  df <- cache(cache_key, {
    aln_query_pileup(aln, intervals, report_mode_str = report_mode_str, clip_mode_str = clip_mode, clip_margin = clip_margin, min_mutations_percent = as.numeric(min_mutations_percent), max_mutations_percent = as.numeric(max_mutations_percent), min_alignment_length = as.integer(min_alignment_length), max_alignment_length = as.integer(max_alignment_length), min_indel_length = as.integer(min_indel_length))
  })

  if (!is.null(df) && nrow(df) > 0) {
    # alntools pileup outputs 1-based coordinates
    df$coord <- df$position
    return(cxt_filter_coords(df))
  }

  return(NULL)
}

align_profile_pileup <- function(profile, aln, gg) {
  intervals <- cxt_get_zoom_view()
    df <- align_query_pileup_mode(aln, report_mode_str = "all", 
                                clip_mode = profile$clip_mode, 
                                clip_margin = profile$clip_margin,
                                min_mutations_percent = as.numeric(profile$min_mutations_percent), 
                                max_mutations_percent = as.numeric(profile$max_mutations_percent),
                                min_alignment_length = as.integer(profile$min_alignment_length),
                                max_alignment_length = as.integer(profile$max_alignment_length),
                                min_indel_length = as.integer(profile$min_indel_length))
  if (is.null(df) || nrow(df) == 0) {
    return(list(plot = gg, legends = list()))
  }
  
  cat(sprintf("plotting pileup with %d variant records\n", nrow(df)))
  
  # ensure gcoord exists (should be added by filter_coords in align_query_pileup_mode)
  if (!"gcoord" %in% colnames(df)) {
    warning("gcoord column missing from pileup data")
    return(list(plot = gg, legends = list()))
  }
  
  # calculate cumulative heights using ave function
  df$y_bottom <- ave(df$count, df$coord, FUN = function(x) cumsum(c(0, x[-length(x)])))
  df$y_top <- df$y_bottom + df$count
  
  # choose color function based on mutation_color_mode
  if (profile$mutation_color_mode == "type") {
    df$fill_color <- get_mutation_type_colors(df$variant)
  } else {
    # detailed mode (default)
    df$fill_color <- get_variant_type_colors(df$variant)
  }
  # position rectangles centered on the coordinate, similar to mutation display in full mode
  df$x_left <- df$gcoord - 0.5
  df$x_right <- df$gcoord + 0.5
  
  # create hover text only if enabled
  if (profile$show_hover) {
    df$hover_text <- paste0(
      "Position: ", df$coord, "\n",
      "Variant: ", df$variant, "\n",
      "Count: ", df$count, " / ", df$coverage, "\n",
      "Fraction: ", sprintf("%.3f", df$count / df$coverage)
    )
  } else {
    df$hover_text <- ""
  }
  
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
  
  # apply force_max_y if set
  if (!is.null(profile$force_max_y) && profile$force_max_y > 0) {
    gg <- gg + ggplot2::coord_cartesian(ylim = c(NA, profile$force_max_y))
  }
  
  # build legends based on mutation_color_mode
  legends <- list()
  mutation_color_mode <- if (!is.null(profile$mutation_color_mode)) profile$mutation_color_mode else "detailed"
  if (mutation_color_mode == "detailed") {
    legend_gg <- create_detailed_variant_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 420, width = 350)))
    }
  } else if (mutation_color_mode == "type") {
    legend_gg <- create_simplified_variant_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 140, width = 300)))
    }
  }
  
  return(list(plot = gg, legends = legends))
}
