# Bin mode for alignment profile

get_current_bin_size <- function(xlim, bin_type, target_bins = 1024) {
  if (bin_type != "auto") {
    bs <- as.numeric(bin_type)
    return(max(1, bs, na.rm = TRUE))
  }

  if (is.null(xlim) || length(xlim) != 2) {
    return(500)
  }

  range_bp <- (xlim[2] + 1) - xlim[1]
  if (range_bp <= 0) {
    return(100)
  }

  # Calculate raw bin size to get close to target_bins
  raw_bin_size <- range_bp / target_bins

  # Ensure minimum bin size of 10, but scale down target_bins for very small windows
  bin_size <- max(10, round(raw_bin_size))
  
  # For very small windows, use fewer target bins to ensure reasonable bin size
  if (range_bp < 200) {
    bin_size <- max(10, round(range_bp / 10))  # aim for ~10 bins minimum
  }
  
  return(bin_size)
}

align_query_bin_mode <- function(aln, cxt, bin_type, target_bins = 1024, seg_threshold = 0.2, non_ref_threshold = 0.9, num_threads = 0, clip_mode = "all", clip_margin = 10, min_mutations_percent = 0.0, max_mutations_percent = 10.0, use_gpu = FALSE) {
  bin_size <- get_current_bin_size(cxt$mapper$xlim, bin_type, target_bins = target_bins)

  # Create cache key based on all relevant parameters
  # Use address of external pointer as unique identifier for alignment
  aln_id <- if (is(aln, "externalptr")) {
    paste0("ptr_", format(aln))
  } else {
    digest::digest(aln, algo = "md5")
  }
  
  cache_key <- paste0("bin_query_",
                     digest::digest(list(
                       aln_id = aln_id,
                       xlim = cxt$mapper$xlim,
                       bin_size = bin_size,
                       seg_threshold = seg_threshold,
                       non_ref_threshold = non_ref_threshold,
                       num_threads = num_threads,
                       clip_mode = clip_mode,
                       clip_margin = clip_margin,
                       min_mutations_percent = min_mutations_percent,
                       max_mutations_percent = max_mutations_percent,
                       use_gpu = use_gpu
                     ), algo = "md5"))

  # Use cache for the bin query
  df <- cache(cache_key, {
    aln_query_bin(aln, cxt$intervals, bin_size, seg_threshold, non_ref_threshold, num_threads, clip_mode_str = clip_mode, clip_margin = clip_margin, min_mutations_percent = as.numeric(min_mutations_percent), max_mutations_percent = as.numeric(max_mutations_percent), use_gpu = use_gpu)
  })
  if (!is.null(df) && nrow(df) > 0) {
    df$start <- df$start + 1
    df$end <- df$end
    return(filter_segments(df, cxt, cxt$mapper$xlim))
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

# Plot segregating sites bins with blue color scale
plot_segregating_sites_bins <- function(gg, df, profile) {
  # segregating sites density visualization
  seg_density_per_100bp <- df$seg_sites_density
  max_seg_density <- max(seg_density_per_100bp, na.rm = TRUE)
  if (is.na(max_seg_density) || max_seg_density == 0) {
    max_seg_density <- 1  # set a default if no segregating sites found
  }
  df$fill_color <- get_color_scale(
    values = seg_density_per_100bp,
    colors = c("#c6dbef", "#08306b"),
    min_val = 0,
    max_val = max_seg_density,
    num_steps = 20
  )
  
  # create hover text only if enabled
  if (profile$show_hover) {
    seg_sites_count <- round(df$seg_sites_density * df$length)
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Segregating sites: ", seg_sites_count
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
  # non-reference sites density visualization
  nonref_density_per_100bp <- df$non_ref_sites_density
  max_nonref_density <- max(nonref_density_per_100bp, na.rm = TRUE)
  if (is.na(max_nonref_density) || max_nonref_density == 0) {
    max_nonref_density <- 1  # set a default if no non-ref sites found
  }
  df$fill_color <- get_color_scale(
    values = nonref_density_per_100bp,
    colors = c("#fff7bc", "#d94701"),
    min_val = 0,
    max_val = max_nonref_density,
    num_steps = 20
  )
  
  # create hover text only if enabled
  if (profile$show_hover) {
    nonref_sites_count <- round(df$non_ref_sites_density * df$length)
    df$hover_text <- paste0(
      "Coverage: ", sprintf("%.2f", df$cov), "x\n",
      "Non-ref sites: ", nonref_sites_count
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

align_profile_bin <- function(profile, cxt, aln, gg) {
  # Get threshold parameters from profile
  seg_threshold <- if (!is.null(profile$seg_threshold)) profile$seg_threshold else 0.2
  non_ref_threshold <- if (!is.null(profile$non_ref_threshold)) profile$non_ref_threshold else 0.9
  num_threads <- if (!is.null(profile$num_threads)) profile$num_threads else 0
  use_gpu <- if (!is.null(profile$use_gpu)) profile$use_gpu else FALSE
  
  df <- align_query_bin_mode(aln, cxt, profile$bin_type, target_bins = profile$target_bins, seg_threshold = seg_threshold, non_ref_threshold = non_ref_threshold, num_threads = num_threads, clip_mode = profile$clip_mode, clip_margin = profile$clip_margin, min_mutations_percent = as.numeric(profile$min_mutations_percent), max_mutations_percent = as.numeric(profile$max_mutations_percent), use_gpu = use_gpu)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
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
  } else if (bin_style == "by_seg_density") {
    gg <- plot_segregating_sites_bins(gg, df, profile)
  } else if (bin_style == "by_nonref_density") {
    gg <- plot_nonref_sites_bins(gg, df, profile)
  } else if (bin_style == "by_genomic_distance") {
    gg <- plot_stacked_mutation_rates(gg, df, profile, normalize = normalize_distrib_bins)
  }

  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)

  return(gg)
}
