# Detail mode for synteny profile
# Shows heatmap with libraries as rows and genomic bins as columns

# helper function to format library names for display
format_library_name <- function(lib_name) {
  # replace underscore with space and clean up formatting
  # e.g., "EBC_-1" becomes "EBC -1"
  gsub("\\.", "-", gsub("_", " ", lib_name))
}

synteny_profile_detail <- function(profile, cxt, data_list, gg, current_binsize) {
  sequenced_bp <- data_list$sequenced_bp
  mutations <- data_list$mutations
  
  # extract library columns (exclude contig, start, end)
  coord_cols <- c("contig", "start", "end")
  lib_cols <- setdiff(colnames(sequenced_bp), coord_cols)
  
  if (length(lib_cols) == 0) {
    warning("no library columns found in synteny data")
    return(list(plot = gg, legends = list()))
  }
  
  cat(sprintf("found %d libraries: %s\n", length(lib_cols), paste(lib_cols[seq_len(min(5, length(lib_cols)))], collapse = ", ")))
  
  # prepare data for heatmap plotting
  plot_data <- prepare_detail_plot_data(sequenced_bp, mutations, lib_cols, profile, current_binsize)
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("no data to plot in detail mode")
    return(list(plot = gg, legends = list()))
  }
  
  # add genomic coordinates using the standard filter_segments function
  plot_data <- filter_segments(plot_data, cxt, cxt$mapper$xlim)
  
  # check again after filtering
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("no data to plot after filtering in detail mode")
    return(list(plot = gg, legends = list()))
  }
  
  # get the actual libraries that appear in the filtered plot data (in correct order)
  # this ensures Y-axis labels match the rectangle positions
  libraries_in_plot <- unique(plot_data$library)
  
  # create heatmap plot
  if (profile$color_style == "none") {
    gg <- plot_detail_gray(gg, plot_data, profile)
  } else if (profile$color_style == "mutations") {
    gg <- plot_detail_mutations(gg, plot_data, profile)
  }
  
  # add library labels on y-axis using the same order as in the plot
  gg <- add_library_labels(gg, libraries_in_plot)
  
  # apply styling
  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank()
    )
  
  return(list(plot = gg, legends = list()))
}

# prepare data for detail plotting
prepare_detail_plot_data <- function(sequenced_bp, mutations, lib_cols, profile, current_binsize) {
  # reshape data from wide to long format
  plot_data <- data.frame()
  
  for (i in seq_len(nrow(sequenced_bp))) {
    row_data <- sequenced_bp[i, ]
    mut_data <- mutations[i, ]
    
    for (lib in lib_cols) {
      seq_bp <- row_data[[lib]]
      mut_density <- mut_data[[lib]]
      
      # handle NA values
      if (is.na(seq_bp)) seq_bp <- 0
      if (is.na(mut_density)) mut_density <- 0
      
      # calculate x-coverage
      xcov <- seq_bp / current_binsize
      
      # check if meets minimum coverage threshold
      meets_threshold <- xcov >= profile$min_xcov
      
      if (meets_threshold) {
        plot_data <- rbind(plot_data, data.frame(
          contig = row_data$contig,
          start = row_data$start,
          end = row_data$end,
          library = lib,
          sequenced_bp = seq_bp,
          mutation_density = mut_density,
          xcov = xcov,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(plot_data)
}

# note: genomic coordinates are now added by filter_segments function

# plot detail mode with gray coloring
plot_detail_gray <- function(gg, plot_data, profile) {
  if (nrow(plot_data) == 0) return(gg)
  
  # use gray color from alignment profile
  gray_color <- "#c0c0c0"
  
  # create library factor for y-positioning
  plot_data$lib_factor <- as.numeric(factor(plot_data$library, levels = unique(plot_data$library)))
  
  # create hover text if enabled
  if (profile$show_hover) {
    plot_data$hover_text <- paste0(
      "Library: ", format_library_name(plot_data$library), "\n",
      "Position: ", plot_data$contig, ":", plot_data$start, "-", plot_data$end, "\n",
      "X-coverage: ", sprintf("%.2f", plot_data$xcov), "\n",
      "Sequenced bp: ", plot_data$sequenced_bp
    )
  } else {
    plot_data$hover_text <- ""
  }
  
  # add rectangles for each library-bin combination
  gg <- gg + ggplot2::geom_rect(
    data = plot_data,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = lib_factor - 0.4, ymax = lib_factor + 0.4,
      fill = gray_color, color = gray_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# plot detail mode with mutation coloring
plot_detail_mutations <- function(gg, plot_data, profile) {
  if (nrow(plot_data) == 0) return(gg)
  
  # use mutation colors from align profile  
  plot_data$fill_color <- get_mutation_colors(plot_data$mutation_density)
  
  # create library factor for y-positioning
  plot_data$lib_factor <- as.numeric(factor(plot_data$library, levels = unique(plot_data$library)))
  
  # create hover text if enabled
  if (profile$show_hover) {
    plot_data$hover_text <- paste0(
      "Library: ", format_library_name(plot_data$library), "\n",
      "Position: ", plot_data$contig, ":", plot_data$start, "-", plot_data$end, "\n",
      "X-coverage: ", sprintf("%.2f", plot_data$xcov), "\n",
      "Mutation density: ", sprintf("%.3f%%", plot_data$mutation_density * 100), "\n",
      "Sequenced bp: ", plot_data$sequenced_bp
    )
  } else {
    plot_data$hover_text <- ""
  }
  
  # add rectangles for each library-bin combination
  gg <- gg + ggplot2::geom_rect(
    data = plot_data,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = lib_factor - 0.4, ymax = lib_factor + 0.4,
      fill = fill_color, color = fill_color,
      text = hover_text
    ),
    size = 0.1
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity()
  
  return(gg)
}

# add library labels to y-axis
add_library_labels <- function(gg, lib_cols) {
  lib_positions <- seq_along(lib_cols)
  
  # format library names for y-axis labels
  formatted_labels <- format_library_name(lib_cols)
  
  gg <- gg + 
    ggplot2::scale_y_continuous(
      breaks = lib_positions,
      labels = formatted_labels,
      limits = c(0.5, length(lib_cols) + 0.5)
    )
  
  return(gg)
}

# genomic coordinates are now handled by filter_segments function