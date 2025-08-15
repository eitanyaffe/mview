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

align_query_bin_mode <- function(aln, cxt, bin_type, target_bins = 1024, seg_threshold = 0.2, non_ref_threshold = 0.9) {
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
                       non_ref_threshold = non_ref_threshold
                     ), algo = "md5"))

  # Use cache for the bin query
  df <- cache(cache_key, {
    aln_query_bin(aln, cxt$intervals, bin_size, seg_threshold, non_ref_threshold)
  })
  
  if (!is.null(df) && nrow(df) > 0) {
    df$start <- df$start + 1
    df$end <- df$end
    return(filter_segments(df, cxt, cxt$mapper$xlim))
  }

  return(NULL)
}

# Function to plot stacked mutation rates
plot_stacked_mutation_rates <- function(gg, df) {
  # Define colors for each mutation distance category (light to dark red)
  mut_colors <- list(
    dist_1_plus = "#67000d",       # darkest red (>1e-1 per bp)
    dist_2 = "#a50f15",           # dark red (1e-2 to 1e-1 per bp)
    dist_3 = "#cb181d",           # medium red (1e-3 to 1e-2 per bp)
    dist_4 = "#fb6a4a",           # light red (1e-4 to 1e-3 per bp)
    dist_5 = "#fc9272",           # lighter red (1e-5 to 1e-4 per bp)
    dist_none = "#fee0d2"         # lightest red (no mutations)
  )
  
  # Create stacked data
  stacked_data <- data.frame()
  
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    
    # Calculate cumulative heights for stacking
    cumulative_height <- 0
    
    # Add each mutation distance category (bottom to top: highest to lowest distances)
    categories <- c("dist_1_plus", "dist_2", "dist_3", 
                   "dist_4", "dist_5", "dist_none")
    
    for (cat in categories) {
      count <- row[[cat]]
      if (is.null(count) || length(count) == 0 || is.na(count)) {
        count <- 0
      }
      if (count > 0) {
        stacked_data <- rbind(stacked_data, data.frame(
          gstart = row$gstart,
          gend = row$gend,
          ymin = cumulative_height,
          ymax = cumulative_height + count,
          category = cat,
          fill_color = mut_colors[[cat]],
          contig = row$contig,
          start = row$start,
          end = row$end,
          length = row$length,
          read_count = row$read_count,
          cov = row$cov,
          count = count,
          stringsAsFactors = FALSE
        ))
        cumulative_height <- cumulative_height + count
      }
    }
  }
  
  if (nrow(stacked_data) > 0) {
    # Convert category names to descriptive text
    category_descriptions <- sapply(stacked_data$category, function(cat) {
      switch(cat,
        "dist_none" = "0 mutations",
        "dist_5" = "0.001% to 0.01% per bp",
        "dist_4" = "0.01% to 0.1% per bp", 
        "dist_3" = "0.1% to 1% per bp",
        "dist_2" = "1% to 10% per bp",
        "dist_1_plus" = ">10% per bp",
        cat  # fallback to original name
      )
    })
    
    # Create hover text for each stack segment
    stacked_data$hover_text <- paste0(
      "Bin Contig: ", stacked_data$contig, "\n",
      "Local (1-based): ", stacked_data$start, "-", stacked_data$end, "\n",
      "Length: ", format_bp(stacked_data$length), "\n",
      "Total Reads: ", stacked_data$read_count, "\n",
      "Category: ", category_descriptions, "\n",
      "Count: ", stacked_data$count
    )
    
    # Add stacked rectangles that touch each other
    gg <- gg + ggplot2::geom_rect(
      data = stacked_data,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = ymin, ymax = ymax,
        fill = fill_color, color = fill_color, text = hover_text
      ), size = 0.1  # thin border in matching color to eliminate gaps
    ) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity()
  }
  
  return(gg)
}

align_profile_bin <- function(profile, cxt, aln, gg) {
  # Get threshold parameters from profile
  seg_threshold <- if (!is.null(profile$seg_threshold)) profile$seg_threshold else 0.2
  non_ref_threshold <- if (!is.null(profile$non_ref_threshold)) profile$non_ref_threshold else 0.9
  
  df <- align_query_bin_mode(aln, cxt, profile$bin_type, target_bins = profile$target_bins, seg_threshold = seg_threshold, non_ref_threshold = non_ref_threshold)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }
  
  # Get bin style from profile parameters
  bin_style <- if (!is.null(profile$bin_style)) profile$bin_style else "by_mut_density"
  
  # Calculate metrics
  df$cov <- ifelse(df$length > 0, df$read_count / df$length, 0)
  df$mut_density <- ifelse(df$read_count > 0, df$mutation_count / df$read_count, 0)

  if (bin_style == "by_mut_density") {
    # Original mutation density visualization
    mutations_per_100bp <- df$mut_density * 100
    df$fill_color <- get_mutation_colors(mutations_per_100bp)
    
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
      ggplot2::scale_fill_identity()
      
  } else if (bin_style == "by_seg_density") {
    # Segregating sites density visualization with green colors
    seg_density_per_100bp <- df$seg_sites_density
    max_seg_density <- max(seg_density_per_100bp, na.rm = TRUE)
    if (is.na(max_seg_density) || max_seg_density == 0) {
      max_seg_density <- 1  # Set a default if no segregating sites found
    }
    df$fill_color <- get_color_scale(
      values = seg_density_per_100bp,
      colors = c("#c6dbef", "#08306b"),
      min_val = 0,
      max_val = max_seg_density,
      num_steps = 20
    )
    
    # Hover text
    df$hover_text <- paste0(
      "Bin Contig: ", df$contig, "\n",
      "Local (1-based): ", df$start, "-", df$end, "\n",
      "Global (0-based): ", df$gstart, "-", df$gend, "\n",
      "Length: ", format_bp(df$length), "\n",
      "Reads: ", df$read_count, " (Cov: ", sprintf("%.2f", df$cov), "x)\n",
      "Seg Sites: ", sprintf("%.4f", df$seg_sites_density), " sites/bp"
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
      ggplot2::scale_fill_identity()
      
  } else if (bin_style == "by_nonref_density") {
    # Non-reference sites density visualization with blue colors
    nonref_density_per_100bp <- df$non_ref_sites_density
    max_nonref_density <- max(nonref_density_per_100bp, na.rm = TRUE)
    if (is.na(max_nonref_density) || max_nonref_density == 0) {
      max_nonref_density <- 1  # Set a default if no non-ref sites found
    }
    df$fill_color <- get_color_scale(
      values = nonref_density_per_100bp,
      colors = c("#fff7bc", "#d94701"),
      min_val = 0,
      max_val = max_nonref_density,
      num_steps = 20
    )
    
    # Hover text
    df$hover_text <- paste0(
      "Bin Contig: ", df$contig, "\n",
      "Local (1-based): ", df$start, "-", df$end, "\n",
      "Global (0-based): ", df$gstart, "-", df$gend, "\n",
      "Length: ", format_bp(df$length), "\n",
      "Reads: ", df$read_count, " (Cov: ", sprintf("%.2f", df$cov), "x)\n",
      "Non-ref Sites: ", sprintf("%.4f", df$non_ref_sites_density), " sites/bp"
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
      ggplot2::scale_fill_identity()
      
  } else if (bin_style == "by_genomic_distance") {
    # Stacked mutation rate visualization
    gg <- plot_stacked_mutation_rates(gg, df)
  }

  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)

  return(gg)
}
