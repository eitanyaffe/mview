# Bin mode for alignment profile

get_current_bin_size <- function(xlim, bin_type, target_bins = 1024, 
    binsizes = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)) {
  
  if (bin_type != "auto") {
    bs <- as.numeric(bin_type)
    if (bs %in% binsizes) {
      return(bs)
    } else {
      warning(sprintf("binsize %d not available, using closest smaller binsize", bs))
      smaller_binsizes <- binsizes[binsizes <= bs]
      if (length(smaller_binsizes) > 0) {
        return(max(smaller_binsizes))
      } else {
        return(min(binsizes))
      }
    }
  }

  if (is.null(xlim) || length(xlim) != 2) {
    return(5000)  # default binsize
  }

  range_bp <- (xlim[2] + 1) - xlim[1]
  if (range_bp <= 0) {
    return(1000)
  }

  # calculate raw bin size to get close to target_bins
  raw_bin_size <- range_bp / target_bins
  
  # choose the closest smaller binsize to ensure we have enough bins
  smaller_binsizes <- binsizes[binsizes <= raw_bin_size]
  if (length(smaller_binsizes) > 0) {
    return(max(smaller_binsizes))
  } else {
    return(min(binsizes))
  }
}

align_query_bin_mode <- function(aln, bin_type, target_bins = 1024, 
  seg_threshold = 0.2, non_ref_threshold = 0.9, num_threads = 0, clip_mode = "all", 
  clip_margin = 10, min_mutations_percent = 0.0, max_mutations_percent = 10.0, min_alignment_length = 0, 
  max_alignment_length = 0, min_indel_length = 3, 
  binsizes = c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)) 
{
  intervals <- cxt_get_zoom_view()
  xlim <- cxt_get_xlim()
  bin_size <- get_current_bin_size(xlim, bin_type, target_bins = target_bins, binsizes = binsizes)
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
  
  cache_key <- paste0("bin_query_",
                     digest::digest(list(
                       aln_id = aln_id,
                       seg_key = seg_key,
                       xlim = xlim,
                       bin_size = bin_size,
                       seg_threshold = seg_threshold,
                       non_ref_threshold = non_ref_threshold,
                       num_threads = num_threads,
                       clip_mode = clip_mode,
                       clip_margin = clip_margin,
                       min_mutations_percent = min_mutations_percent,
                       max_mutations_percent = max_mutations_percent,
                       min_alignment_length = min_alignment_length,
                       max_alignment_length = max_alignment_length,
                       min_indel_length = min_indel_length
                     ), algo = "md5"))

  
  # Use cache for the bin query
  df <- cache(cache_key, {
    aln_query_bin(aln, intervals, bin_size, seg_threshold, 
      non_ref_threshold, num_threads, 
      clip_mode_str = clip_mode, 
      clip_margin = clip_margin, 
      min_mutations_percent = as.numeric(min_mutations_percent), 
      max_mutations_percent = as.numeric(max_mutations_percent), 
      min_alignment_length = as.integer(min_alignment_length), 
      max_alignment_length = as.integer(max_alignment_length), 
      min_indel_length = as.integer(min_indel_length))
  })

  if (!is.null(df) && nrow(df) > 0) {
    df$start <- df$start + 1
    df$end <- df$end
    return(cxt_filter_intervals(df))
  }

  return(NULL)
}

# Function to plot stacked mutation rates
plot_stacked_mutation_rates <- function(gg, df, profile, normalize = FALSE) {
  if (nrow(df) == 0) return(gg)
  
  # use shared red color scale for discrete mutation categories
  red_colors <- get_shared_red_scale(continuous = FALSE, num_colors = 6)
  
  # define categories in stacking order (bottom to top)
  categories <- c("dist_1_plus", "dist_2", "dist_3", "dist_4", "dist_5", "dist_none")
  colors <- c(red_colors[6], red_colors[5], red_colors[4], red_colors[3], red_colors[2], red_colors[1])
  
  # descriptions for hover text
  descriptions <- c(">10% per bp", "1% to 10% per bp", "0.1% to 1% per bp", 
                   "0.01% to 0.1% per bp", "0.001% to 0.01% per bp", "0 mutations")
  
  # normalize counts to percentages if requested
  if (normalize) {
    # calculate total counts for each bin
    total_counts <- rep(0, nrow(df))
    for (cat in categories) {
      cat_counts <- df[[cat]]
      if (!is.null(cat_counts)) {
        cat_counts[is.na(cat_counts)] <- 0
        total_counts <- total_counts + cat_counts
      }
    }
    
    # normalize each category to percentage of total
    for (cat in categories) {
      if (!is.null(df[[cat]])) {
        df[[cat]] <- ifelse(total_counts > 0, (df[[cat]] / total_counts) * 100, 0)
      }
    }
  }
  
  # initialize bottom y values for each bin (start at 0)
  bottom_y <- rep(0, nrow(df))
  
  # create stacked data by processing each category
  stacked_data <- data.frame()
  
  for (i in seq_along(categories)) {
    cat <- categories[i]
    color <- colors[i]
    desc <- descriptions[i]
    
    # get counts for this category, handle missing/NA values
    counts <- df[[cat]]
    if (is.null(counts)) counts <- rep(0, nrow(df))
    counts[is.na(counts)] <- 0
    
    # only create rectangles for non-zero counts
    has_count <- counts > 0
    if (any(has_count)) {
      # create data for this category's rectangles
      category_data <- data.frame(
        gstart = df$gstart[has_count],
        gend = df$gend[has_count],
        ymin = bottom_y[has_count],
        ymax = bottom_y[has_count] + counts[has_count],
        fill_color = color,
        category = cat,
        description = desc,
        count = counts[has_count],
        contig = df$contig[has_count],
        start = df$start[has_count],
        end = df$end[has_count],
        length = df$length[has_count],
        read_count = df$read_count[has_count],
        sequenced_bp = df$sequenced_bp[has_count],
        stringsAsFactors = FALSE
      )
      
      # create hover text only if enabled
      if (profile$show_hover) {
        if (normalize) {
          category_data$hover_text <- paste0(
            sprintf("%.1f%%", category_data$count), " of reads\n",
            "Category: ", category_data$description, "\n"
          )
        } else {
          category_data$hover_text <- paste0(
            category_data$count, " out of ", category_data$read_count, " reads\n",
            "Category: ", category_data$description, "\n"
          )
        }
      } else {
        category_data$hover_text <- ""
      }
      
      stacked_data <- rbind(stacked_data, category_data)
    }
    
    # update bottom y values for next category
    bottom_y <- bottom_y + counts
  }
  
  if (nrow(stacked_data) > 0) {
    # add stacked rectangles
    gg <- gg + ggplot2::geom_rect(
      data = stacked_data,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = ymin, ymax = ymax,
        fill = fill_color, color = fill_color,
        text = hover_text
      ),
      size = 0.1
    ) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity()
  }
  
  return(gg)
}

# Plot mutation density bins with red color scale
plot_mutation_density_bins <- function(gg, df, profile) {
  # mutation density visualization
  df$fill_color <- get_mutation_colors(df$mut_density)
  
  # create hover text only if enabled
  if (profile$show_hover) {
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Mutations: ", df$mutation_count, " / ", df$sequenced_bp, " bp\n",
      "Mutation density: ", sprintf("%.3f%%", df$mut_density * 100)
    )
  } else {
    df$hover_text <- ""
  }
  
  # plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# Plot median mutation density bins with red color scale (same as mutation density)
plot_median_mutation_density_bins <- function(gg, df, profile) {
  # median mutation density visualization
  df$fill_color <- get_mutation_colors(df$median_mutation_density)
  
  # create hover text only if enabled
  if (profile$show_hover) {
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Mutations: ", df$mutation_count, " / ", df$sequenced_bp, " bp\n",
      "Median mutation density: ", sprintf("%.3f%%", df$median_mutation_density * 100)
    )
  } else {
    df$hover_text <- ""
  }
  
  # plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# Plot segregating sites bins with blue color scale
plot_segregating_sites_bins <- function(gg, df, profile) {
  # segregating sites density visualization using fixed/auto scale
  seg_density_per_100bp <- df$seg_sites_density
  max_density_scale <- profile$max_col_dist_percent
  if (is.null(max_density_scale) || identical(max_density_scale, "auto")) {
    max_density_scale <- max(seg_density_per_100bp, na.rm = TRUE)
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  } else {
    max_density_scale <- as.numeric(max_density_scale) / 100
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  }
  df$fill_color <- get_color_scale(
    values = seg_density_per_100bp,
    colors = c("#c6dbef", "#08306b"),
    min_val = 0,
    max_val = max_density_scale,
    num_steps = 20
  )
  
  # create hover text only if enabled
  if (profile$show_hover) {
    seg_sites_count <- round(df$seg_sites_density * df$length)
    seg_pct <- df$seg_sites_density * 100
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Segregating sites: ", seg_sites_count, " out of ", df$length, " bp\n",
      "Percent: ", sprintf("%.3f%%", seg_pct)
    )
  } else {
    df$hover_text <- ""
  }
  
  # plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# Plot non-reference sites bins with orange color scale
plot_nonref_sites_bins <- function(gg, df, profile) {
  # non-reference sites density visualization using fixed/auto scale
  nonref_density_per_100bp <- df$non_ref_sites_density
  max_density_scale <- profile$max_col_dist_percent
  if (is.null(max_density_scale) || identical(max_density_scale, "auto")) {
    max_density_scale <- max(nonref_density_per_100bp, na.rm = TRUE)
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  } else {
    max_density_scale <- as.numeric(max_density_scale) / 100
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  }
  df$fill_color <- get_color_scale(
    values = nonref_density_per_100bp,
    colors = c("#fff7bc", "#d94701"),
    min_val = 0,
    max_val = max_density_scale,
    num_steps = 20
  )
  
  # create hover text only if enabled
  if (profile$show_hover) {
    nonref_sites_count <- round(df$non_ref_sites_density * df$length)
    nonref_pct <- df$non_ref_sites_density * 100
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Non-ref sites: ", nonref_sites_count, " out of ", df$length, " bp\n",
      "Percent: ", sprintf("%.3f%%", nonref_pct)
    )
  } else {
    df$hover_text <- ""
  }
  
  # plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# Plot segregating clip sites bins with purple color scale
plot_seg_clip_sites_bins <- function(gg, df, profile) {
  # segregating clip sites density visualization using fixed/auto scale
  seg_clip_density_per_100bp <- df$seg_clip_density
  max_density_scale <- profile$max_col_dist_percent
  if (is.null(max_density_scale) || identical(max_density_scale, "auto")) {
    max_density_scale <- max(seg_clip_density_per_100bp, na.rm = TRUE)
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  } else {
    max_density_scale <- as.numeric(max_density_scale) / 100
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  }
  df$fill_color <- get_color_scale(
    values = seg_clip_density_per_100bp,
    colors = c("#f1eef6", "#6a51a3"),
    min_val = 0,
    max_val = max_density_scale,
    num_steps = 20
  )
  
  # create hover text only if enabled
  if (profile$show_hover) {
    seg_clip_sites_count <- round(df$seg_clip_density * df$length)
    seg_clip_pct <- df$seg_clip_density * 100
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Seg clip sites: ", seg_clip_sites_count, " out of ", df$length, " bp\n",
      "Percent: ", sprintf("%.3f%%", seg_clip_pct)
    )
  } else {
    df$hover_text <- ""
  }
  
  # plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# Plot non-reference clip sites bins with green color scale
plot_nonref_clip_sites_bins <- function(gg, df, profile) {
  # non-reference clip sites density visualization using fixed/auto scale
  nonref_clip_density_per_100bp <- df$non_ref_clip_density
  max_density_scale <- profile$max_col_dist_percent
  if (is.null(max_density_scale) || identical(max_density_scale, "auto")) {
    max_density_scale <- max(nonref_clip_density_per_100bp, na.rm = TRUE)
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  } else {
    max_density_scale <- as.numeric(max_density_scale) / 100
    if (!is.finite(max_density_scale) || max_density_scale <= 0) max_density_scale <- 1.0
  }
  df$fill_color <- get_color_scale(
    values = nonref_clip_density_per_100bp,
    colors = c("#edf8e9", "#238b45"),
    min_val = 0,
    max_val = max_density_scale,
    num_steps = 20
  )
  
  # create hover text only if enabled
  if (profile$show_hover) {
    nonref_clip_sites_count <- round(df$non_ref_clip_density * df$length)
    nonref_clip_pct <- df$non_ref_clip_density * 100
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Non-ref clip sites: ", nonref_clip_sites_count, " out of ", df$length, " bp\n",
      "Percent: ", sprintf("%.3f%%", nonref_clip_pct)
    )
  } else {
    df$hover_text <- ""
  }
  
  # plot bins
  gg <- gg + ggplot2::geom_rect(
    data = df,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = cov,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# create legends for bin profile based on bin style
create_bin_profile_legends <- function(bin_style, profile, df) {
  legends <- list()
  
  if (bin_style == "by_mut_density") {
    legend_gg <- create_mutation_density_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 190, width = 300, title = "mutations per bp")))
    }
  } else if (bin_style == "by_median_mutation_density") {
    legend_gg <- create_mutation_density_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 190, width = 300, title = "median mutations per bp")))
    }
  } else if (bin_style == "by_seg_density") {
    max_density <- profile$max_col_dist_percent
    if (is.null(max_density) || identical(max_density, "auto")) {
      max_density <- max(df$seg_sites_density, na.rm = TRUE)
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    } else {
      max_density <- as.numeric(max_density) / 100
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    }
    color_defs <- get_alignment_color_definitions()
    legend_gg <- create_gradient_legend("segregating sites (per 100 bp)", colors = color_defs$seg_gradient, max_val = max_density, n_steps = 10, as_percent = TRUE)
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 190, width = 300, title = "segregating sites density")))
    }
  } else if (bin_style == "by_nonref_density") {
    max_density <- profile$max_col_dist_percent
    if (is.null(max_density) || identical(max_density, "auto")) {
      max_density <- max(df$non_ref_sites_density, na.rm = TRUE)
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    } else {
      max_density <- as.numeric(max_density) / 100
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    }
    color_defs <- get_alignment_color_definitions()
    legend_gg <- create_gradient_legend("non-ref sites (per 100 bp)", colors = color_defs$nonref_gradient, max_val = max_density, n_steps = 10, as_percent = TRUE)
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 190, width = 300, title = "non-ref sites density")))
    }
  } else if (bin_style == "by_seg_clip_density") {
    max_density <- profile$max_col_dist_percent
    if (is.null(max_density) || identical(max_density, "auto")) {
      max_density <- max(df$seg_clip_density, na.rm = TRUE)
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    } else {
      max_density <- as.numeric(max_density) / 100
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    }
    legend_gg <- create_gradient_legend("segregating clip sites (per 100 bp)", colors = c("#f1eef6", "#6a51a3"), max_val = max_density, n_steps = 10, as_percent = TRUE)
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 190, width = 300, title = "segregating clip sites density")))
    }
  } else if (bin_style == "by_non_ref_clip_density") {
    max_density <- profile$max_col_dist_percent
    if (is.null(max_density) || identical(max_density, "auto")) {
      max_density <- max(df$non_ref_clip_density, na.rm = TRUE)
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    } else {
      max_density <- as.numeric(max_density) / 100
      if (!is.finite(max_density) || max_density <= 0) max_density <- 1.0
    }
    legend_gg <- create_gradient_legend("non-ref clip sites (per 100 bp)", colors = c("#edf8e9", "#238b45"), max_val = max_density, n_steps = 10, as_percent = TRUE)
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 190, width = 300, title = "non-ref clip sites density")))
    }
  } else if (bin_style == "by_genomic_distance") {
    legend_gg <- create_stacked_mutation_rates_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 210, width = 320, title = "mutation rate bins")))
    }
  }
  
  return(legends)
}

align_profile_bin <- function(profile, aln, gg) {
  intervals <- cxt_get_zoom_view()
  
  # Get threshold parameters from profile
  seg_threshold <- if (!is.null(profile$seg_threshold)) profile$seg_threshold else 0.2
  non_ref_threshold <- if (!is.null(profile$non_ref_threshold)) profile$non_ref_threshold else 0.9
  num_threads <- if (!is.null(profile$num_threads)) profile$num_threads else 0
  
  df <- align_query_bin_mode(aln, profile$bin_type, target_bins = profile$target_bins, seg_threshold = seg_threshold, non_ref_threshold = non_ref_threshold, num_threads = num_threads, clip_mode = profile$clip_mode, clip_margin = profile$clip_margin, min_mutations_percent = as.numeric(profile$min_mutations_percent), max_mutations_percent = as.numeric(profile$max_mutations_percent), min_alignment_length = as.integer(profile$min_alignment_length), max_alignment_length = as.integer(profile$max_alignment_length), min_indel_length = as.integer(profile$min_indel_length), binsizes = profile$binsizes)
  if (is.null(df) || nrow(df) == 0) {
    return(list(plot = gg, legends = list()))
  }
  
  # eliminate gaps between bins by adjusting boundaries
  if (nrow(df) > 1) {
    # sort by gstart to ensure proper ordering
    df <- df[order(df$gstart), ]
    
    # adjust gend of each bin to touch the start of the next bin
    for (i in seq_len(nrow(df) - 1)) {
      if (df$gend[i] < df$gstart[i + 1]) {
        # extend current bin to touch next bin
        df$gend[i] <- df$gstart[i + 1]
      }
    }
  }
  
  # Get bin style from profile parameters
  bin_style <- if (!is.null(profile$bin_style)) profile$bin_style else "by_mut_density"
  
  # Get normalization parameter
  normalize_distrib_bins <- if (!is.null(profile$normalize_distrib_bins)) profile$normalize_distrib_bins else FALSE
  
  # Calculate metrics
  df$cov <- ifelse(df$length > 0, df$sequenced_bp / df$length, 0)
  df$mut_density <- ifelse(df$sequenced_bp > 0, df$mutation_count / df$sequenced_bp, 0)

  # choose and call appropriate visualization function
  if (bin_style == "by_mut_density") {
    gg <- plot_mutation_density_bins(gg, df, profile)
  } else if (bin_style == "by_median_mutation_density") {
    gg <- plot_median_mutation_density_bins(gg, df, profile)
  } else if (bin_style == "by_seg_density") {
    gg <- plot_segregating_sites_bins(gg, df, profile)
  } else if (bin_style == "by_nonref_density") {
    gg <- plot_nonref_sites_bins(gg, df, profile)
  } else if (bin_style == "by_seg_clip_density") {
    gg <- plot_seg_clip_sites_bins(gg, df, profile)
  } else if (bin_style == "by_non_ref_clip_density") {
    gg <- plot_nonref_clip_sites_bins(gg, df, profile)
  } else if (bin_style == "by_genomic_distance") {
    gg <- plot_stacked_mutation_rates(gg, df, profile, normalize = normalize_distrib_bins)
  }

  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)

  # apply force_max_y if set
  if (!is.null(profile$force_max_y) && profile$force_max_y > 0) {
    gg <- gg + ggplot2::coord_cartesian(ylim = c(NA, profile$force_max_y))
  }

  # create legends for legend tab
  legends <- create_bin_profile_legends(bin_style, profile, df)
  
  return(list(plot = gg, legends = legends))
}
