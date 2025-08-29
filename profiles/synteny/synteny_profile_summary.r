# Summary mode for synteny profile  
# Shows count of libraries per bin that meet coverage threshold

synteny_profile_summary <- function(profile, cxt, data_list, gg, current_binsize) {
  # filter data to current view first using new function
  sequenced_bp <- filter_synteny_matrix(data_list$sequenced_bp, cxt)
  mutations <- filter_synteny_matrix(data_list$mutations, cxt)
  
  if (is.null(sequenced_bp) || is.null(mutations) || nrow(sequenced_bp) == 0) {
    warning("no data in current view for summary mode")
    return(list(plot = gg, legends = list()))
  }
  
  # extract library columns (exclude contig, start, end)
  coord_cols <- c("contig", "start", "end")
  lib_cols <- setdiff(colnames(sequenced_bp), coord_cols)
  
  if (length(lib_cols) == 0) {
    warning("no library columns found in synteny data")
    return(list(plot = gg, legends = list()))
  }
  
  cat(sprintf("found %d libraries for summary: %s\n", length(lib_cols), paste(lib_cols[seq_len(min(5, length(lib_cols)))], collapse = ", ")))
  
  # prepare summary data
  summary_data <- prepare_summary_data(sequenced_bp, mutations, lib_cols, profile, current_binsize)
  
  if (is.null(summary_data) || nrow(summary_data) == 0) {
    warning("no data to plot in summary mode")
    return(list(plot = gg, legends = list()))
  }
  
  # add global coordinates to summary data
  summary_data <- filter_segments(summary_data, cxt, cxt$mapper$xlim)
  
  if (is.null(summary_data) || nrow(summary_data) == 0) {
    warning("no summary data after coordinate mapping")
    return(list(plot = gg, legends = list()))
  }
  
  # create summary plot
  if (profile$color_style == "none") {
    gg <- plot_summary_gray(gg, summary_data, profile)
  } else if (profile$color_style == "mutations") {
    gg <- plot_summary_stacked_mutations(gg, summary_data, profile)
  }
  
  # apply styling
  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    ggplot2::labs(y = "library count")
  
  # create legends
  legends <- list()
  if (profile$color_style == "mutations") {
    legend_gg <- create_stacked_mutation_rates_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 210, width = 320, title = "libraries per mutation rate")))
    }
  }
  
  return(list(plot = gg, legends = legends))
}

# helper function to categorize mutation density
categorize_mutation_density <- function(mut_density) {
  if (mut_density <= 0) return("dist_none")
  if (mut_density <= 0.0001) return("dist_5")     # ≤ 0.01%
  if (mut_density <= 0.001) return("dist_4")      # ≤ 0.1%
  if (mut_density <= 0.01) return("dist_3")       # ≤ 1%
  if (mut_density <= 0.1) return("dist_2")        # ≤ 10%
  return("dist_1_plus")                            # > 10%
}

# prepare summary data by counting libraries per bin and mutation category
prepare_summary_data <- function(sequenced_bp, mutations, lib_cols, profile, current_binsize) {
  summary_data <- data.frame()
  
  for (i in seq_len(nrow(sequenced_bp))) {
    row_data <- sequenced_bp[i, ]
    mut_data <- mutations[i, ]
    
    # initialize mutation category counts
    mutation_counts <- list(
      dist_none = 0,
      dist_5 = 0,
      dist_4 = 0, 
      dist_3 = 0,
      dist_2 = 0,
      dist_1_plus = 0
    )
    
    total_lib_count <- 0
    mutation_densities <- numeric()
    
    for (lib in lib_cols) {
      seq_bp <- row_data[[lib]]
      mut_density <- mut_data[[lib]]
      
      # handle NA values
      if (is.na(seq_bp)) seq_bp <- 0
      if (is.na(mut_density)) mut_density <- 0
      
      # calculate x-coverage
      xcov <- seq_bp / current_binsize
      
      # check if meets minimum coverage threshold
      if (xcov >= profile$min_xcov) {
        total_lib_count <- total_lib_count + 1
        mutation_densities <- c(mutation_densities, mut_density)
        
        # categorize this library's mutation density
        category <- categorize_mutation_density(mut_density)
        mutation_counts[[category]] <- mutation_counts[[category]] + 1
      }
    }
    
    # calculate average mutation density for libraries that meet threshold
    avg_mutation_density <- 0
    if (length(mutation_densities) > 0) {
      avg_mutation_density <- mean(mutation_densities, na.rm = TRUE)
    }
    
    # only include bins with at least one qualifying library
    if (total_lib_count > 0) {
      summary_data <- rbind(summary_data, data.frame(
        contig = row_data$contig,
        start = row_data$start,
        end = row_data$end,
        lib_count = total_lib_count,
        avg_mutation_density = avg_mutation_density,
        dist_none = mutation_counts$dist_none,
        dist_5 = mutation_counts$dist_5,
        dist_4 = mutation_counts$dist_4,
        dist_3 = mutation_counts$dist_3,
        dist_2 = mutation_counts$dist_2,
        dist_1_plus = mutation_counts$dist_1_plus,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(summary_data)
}

# note: genomic coordinates are now added by filter_segments function

# plot summary mode with gray coloring
plot_summary_gray <- function(gg, summary_data, profile) {
  if (nrow(summary_data) == 0) return(gg)
  
  # use gray color from alignment profile
  gray_color <- "#c0c0c0"
  
  # create hover text if enabled
  if (profile$show_hover) {
    summary_data$hover_text <- paste0(
      "Position: ", summary_data$contig, ":", summary_data$start, "-", summary_data$end, "\n",
      "Libraries with coverage: ", summary_data$lib_count, "\n",
      "Avg mutation density: ", sprintf("%.3f%%", summary_data$avg_mutation_density * 100)
    )
  } else {
    summary_data$hover_text <- ""
  }
  
  # eliminate gaps between bins by adjusting boundaries
  if (nrow(summary_data) > 1) {
    # sort by gstart to ensure proper ordering
    summary_data <- summary_data[order(summary_data$gstart), ]
    
    # adjust gend of each bin to touch the start of the next bin
    for (i in seq_len(nrow(summary_data) - 1)) {
      if (summary_data$gend[i] < summary_data$gstart[i + 1]) {
        # extend current bin to touch next bin
        summary_data$gend[i] <- summary_data$gstart[i + 1]
      }
    }
  }
  
  # add rectangles for library counts
  gg <- gg + ggplot2::geom_rect(
    data = summary_data,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = 0, ymax = lib_count,
      fill = gray_color, color = gray_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# plot summary mode with stacked mutation coloring
plot_summary_stacked_mutations <- function(gg, summary_data, profile) {
  if (nrow(summary_data) == 0) return(gg)
  
  # use shared red color scale from alignment profile (same as stacked mutation rates)
  red_colors <- get_shared_red_scale(continuous = FALSE, num_colors = 6)
  
  # define categories in stacking order (bottom to top)
  categories <- c("dist_1_plus", "dist_2", "dist_3", "dist_4", "dist_5", "dist_none")
  colors <- c(red_colors[6], red_colors[5], red_colors[4], red_colors[3], red_colors[2], red_colors[1])
  
  # descriptions for hover text
  descriptions <- c(">10% per bp", "1% to 10% per bp", "0.1% to 1% per bp", 
                   "0.01% to 0.1% per bp", "0.001% to 0.01% per bp", "0 mutations")
  
  # eliminate gaps between bins by adjusting boundaries
  if (nrow(summary_data) > 1) {
    # sort by gstart to ensure proper ordering
    summary_data <- summary_data[order(summary_data$gstart), ]
    
    # adjust gend of each bin to touch the start of the next bin
    for (i in seq_len(nrow(summary_data) - 1)) {
      if (summary_data$gend[i] < summary_data$gstart[i + 1]) {
        # extend current bin to touch next bin
        summary_data$gend[i] <- summary_data$gstart[i + 1]
      }
    }
  }
  
  # initialize bottom y values for each bin (start at 0)
  bottom_y <- rep(0, nrow(summary_data))
  
  # create stacked data by processing each category
  stacked_data <- data.frame()
  
  for (i in seq_along(categories)) {
    cat <- categories[i]
    color <- colors[i]
    desc <- descriptions[i]
    
    # get counts for this category, handle missing/NA values
    counts <- summary_data[[cat]]
    if (is.null(counts)) counts <- rep(0, nrow(summary_data))
    counts[is.na(counts)] <- 0
    
    # only create rectangles for non-zero counts
    has_count <- counts > 0
    if (any(has_count)) {
      # create data for this category's rectangles
      category_data <- data.frame(
        gstart = summary_data$gstart[has_count],
        gend = summary_data$gend[has_count],
        ymin = bottom_y[has_count],
        ymax = bottom_y[has_count] + counts[has_count],
        fill_color = color,
        category = cat,
        description = desc,
        count = counts[has_count],
        contig = summary_data$contig[has_count],
        start = summary_data$start[has_count],
        end = summary_data$end[has_count],
        lib_count = summary_data$lib_count[has_count],
        avg_mutation_density = summary_data$avg_mutation_density[has_count],
        stringsAsFactors = FALSE
      )
      
      # create hover text
      if (profile$show_hover) {
        # both total and segment hover information
        category_data$hover_text <- paste0(
          "Position: ", category_data$contig, ":", category_data$start, "-", category_data$end, "\n",
          "Total libraries: ", category_data$lib_count, "\n",
          "This segment: ", category_data$count, " libraries\n",
          "Category: ", category_data$description
        )
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
