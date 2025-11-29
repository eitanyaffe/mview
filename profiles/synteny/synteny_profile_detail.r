# Detail mode for synteny profile
# Shows heatmap with libraries as rows and genomic bins as columns

# Load shared utilities for mutation coloring and legends
source("profiles/align/align_utils.r")

# helper function to format library names for display
format_library_name <- function(lib_name) {
  # replace underscore with space and clean up formatting
  # e.g., "EBC_-1" becomes "EBC -1"
  gsub("\\.", "-", gsub("_", " ", lib_name))
}

# create hover text for detail bins
create_detail_hover_text <- function(plot_data, profile, include_mutations = FALSE) {
  if (!profile$show_hover) {
    return("")
  }
  
  # base hover text components
  hover_text <- paste0(
    "Library: ", format_library_name(plot_data$library), "\n",
    "Position: ", plot_data$contig, ":", plot_data$start, "-", plot_data$end, "\n",
    "X-coverage: ", sprintf("%.2f", plot_data$xcov)
  )
  
  # add mutation density if requested
  if (include_mutations) {
    hover_text <- paste0(hover_text, "\n",
      "Mutation density: ", sprintf("%.3f%%", plot_data$mutation_density * 100)
    )
  }
  
  # add sequenced bp at the end
  hover_text <- paste0(hover_text, "\n",
    "Sequenced bp: ", plot_data$sequenced_bp
  )
  
  return(hover_text)
}

synteny_profile_detail <- function(profile, data_list, gg, current_binsize) {
  # filter data to current view first using new function
  sequenced_bp <- filter_synteny_matrix(data_list$sequenced_bp)
  mutations <- filter_synteny_matrix(data_list$mutations)
  
  if (is.null(sequenced_bp) || is.null(mutations) || nrow(sequenced_bp) == 0) {
    warning("no data in current view for detail mode")
    return(list(plot = gg, legends = list()))
  }
  
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
  
  # add global coordinates to plot data
  plot_data <- cxt_filter_intervals(plot_data)
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("no plot data after coordinate mapping")
    return(list(plot = gg, legends = list()))
  }
  
  # get the actual libraries that appear in the filtered plot data (in correct order)
  # this ensures Y-axis labels match the rectangle positions
  libraries_in_plot <- unique(plot_data$library)
  
  bin_gap = 0.1
  # create heatmap plot
  if (profile$color_style == "none") {
    gg <- plot_detail_bins(gg, plot_data, profile, color_mode = "gray", bin_gap = bin_gap)
  } else if (profile$color_style == "mutations") {
    gg <- plot_detail_bins(gg, plot_data, profile, color_mode = "mutations", bin_gap = bin_gap)
  }
  
  # check if we should overlay consensus mutations
      
  should_show_mutations <- should_show_consensus_mutations(profile)

  if (should_show_mutations && profile$show_mutations) {
    # only load consensus data when we actually need it
    consensus_data <- load_consensus_data(profile$consensus_f, current_binsize)
    if (!is.null(consensus_data)) {
        mutation_plot_data <- prepare_consensus_mutation_data(plot_data, consensus_data, profile, current_binsize)
        if (!is.null(mutation_plot_data) && nrow(mutation_plot_data) > 0) {
          # build hover text for synteny-specific format
          if (profile$show_hover) {
            mutation_plot_data$hover_text <- paste0(
              "Library: ", format_library_name(mutation_plot_data$library), "\n",
              "Position: ", mutation_plot_data$contig, ":", mutation_plot_data$coord, "\n",
              "Variant: ", mutation_plot_data$variant_desc, "\n",
              "Frequency: ", sprintf("%.2f%%", mutation_plot_data$frequency * 100), "\n",
              "Coverage: ", mutation_plot_data$coverage
            )
          } else {
            mutation_plot_data$hover_text <- ""
          }
          
          # map columns to align profile format
          mutation_plot_data$desc <- mutation_plot_data$variant_desc
          mutation_plot_data$ybottom <- mutation_plot_data$y_pos - 0.5 + bin_gap
          mutation_plot_data$ytop <- mutation_plot_data$y_pos + 0.5 - bin_gap
          mutation_plot_data$gcoord <- mutation_plot_data$coord
          
          gg <- plot_mutations_unified(gg, mutation_plot_data, profile)
        }
    }
  }
  
  # add library labels on y-axis using the same order as in the plot
  gg <- add_library_labels(gg, libraries_in_plot)
  
  # apply styling
  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank()
    )
  
  # create legends
  legends <- list()
  if (profile$color_style == "mutations") {
    legend_gg <- create_mutation_density_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 210, width = 320, title = "mutation density")))
    }
  }
  
  # add variant legend if mutations are shown
  if (should_show_mutations && profile$show_mutations) {
    # check if consensus data is available (don't reload it, just check)
    consensus_available <- !is.null(profile$consensus_f)
    if (consensus_available) {
      if (profile$mutation_color_mode == "detailed") {
        legend_gg <- create_detailed_variant_legend()
        if (!is.null(legend_gg)) {
          legends <- c(legends, list(list(gg = legend_gg, height = 420, width = 350)))
        }
      } else if (profile$mutation_color_mode == "type") {
        legend_gg <- create_simplified_variant_legend()
        if (!is.null(legend_gg)) {
          legends <- c(legends, list(list(gg = legend_gg, height = 140, width = 300)))
        }
      }
    }
  }
  
  return(list(plot = gg, legends = legends))
}

# prepare data for detail plotting
prepare_detail_plot_data <- function(sequenced_bp, mutations, lib_cols, profile, current_binsize) {
  # apply lib_grepl filtering to library columns early
  if (!is.null(profile$lib_grepl) && profile$lib_grepl != "") {
    lib_cols <- lib_cols[grepl(profile$lib_grepl, lib_cols)]
    if (length(lib_cols) == 0) {
      cat(sprintf("no libraries match grepl pattern '%s'\n", profile$lib_grepl))
      return(data.frame())
    }
    cat(sprintf("filtered to %d libraries matching pattern '%s'\n", length(lib_cols), profile$lib_grepl))
  }
  
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

# note: genomic coordinates are now added by filter_intervals function

# unified function to plot detail bins with configurable gap
plot_detail_bins <- function(gg, plot_data, profile, color_mode = "gray", bin_gap = 0.1) {
  if (nrow(plot_data) == 0) return(gg)
  
  # create library factor for y-positioning
  plot_data$lib_factor <- as.numeric(factor(plot_data$library, levels = unique(plot_data$library)))
  
  # calculate bin boundaries with gap
  plot_data$ybottom <- plot_data$lib_factor - 0.5 + bin_gap
  plot_data$ytop <- plot_data$lib_factor + 0.5 - bin_gap
  
  # set colors and hover text based on mode
  if (color_mode == "gray") {
    plot_data$fill_color <- "#c0c0c0"
    plot_data$hover_text <- create_detail_hover_text(plot_data, profile, include_mutations = FALSE)
  } else { # mutations mode
    plot_data$fill_color <- get_mutation_colors(plot_data$mutation_density)
    plot_data$hover_text <- create_detail_hover_text(plot_data, profile, include_mutations = TRUE)
  }
  
  # add rectangles for each library-bin combination
  gg <- gg + ggplot2::geom_rect(
    data = plot_data,
    ggplot2::aes(
      xmin = gstart, xmax = gend,
      ymin = ybottom, ymax = ytop,
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
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.8))  # 20% smaller
    )
  
  return(gg)
}

# check if consensus mutations should be shown based on zoom level
should_show_consensus_mutations <- function(profile) {
  xlim <- cxt_get_xlim()
  if (is.null(xlim) || length(xlim) != 2) {
    return(FALSE)
  }
  
  range_bp <- (xlim[2] + 1) - xlim[1]
  return(range_bp <= profile$full_threshold)
}



# prepare consensus mutation data for plotting with efficient filtering
prepare_consensus_mutation_data <- function(plot_data, consensus_data, profile, current_binsize) {
  if (is.null(consensus_data) || nrow(consensus_data) == 0) {
    return(NULL)
  }
  
  # prepare bin keys from plot_data (already filtered and final)
  plot_data$library <- gsub("\\.", "-", plot_data$library)
  bin_keys <- paste(plot_data$library, plot_data$contig, plot_data$start, sep="_")
  
  # efficient vectorized filtering using existing keys
  cat(sprintf("filtering %d mutations using %d bin keys\n", nrow(consensus_data), length(bin_keys)))
  filtered_mutations <- consensus_data[consensus_data$key %in% bin_keys, ]
  cat(sprintf("filtered to %d mutations in current view\n", nrow(filtered_mutations)))
  
  if (nrow(filtered_mutations) == 0) {
    return(NULL)
  }
  
  # get actual libraries from filtered plot data for factor levels
  libraries_in_plot <- unique(plot_data$library)
  lib_factor <- as.numeric(factor(filtered_mutations$lib_id, levels = libraries_in_plot))
  
  # prepare final data structure
  mutation_plot_data <- data.frame(
    library = filtered_mutations$lib_id,
    contig = filtered_mutations$contig,
    # alntools consensus mode outputs 1-based coordinates
    coord = filtered_mutations$position,
    variant_type = filtered_mutations$variant_type,
    variant_desc = filtered_mutations$variant_desc,
    count = filtered_mutations$count,
    coverage = filtered_mutations$coverage,
    frequency = filtered_mutations$frequency,
    y_pos = lib_factor,
    stringsAsFactors = FALSE
  )
  
  # map genomic coordinates
  mutation_plot_data <- cxt_filter_coords(mutation_plot_data)
  
  return(mutation_plot_data)
}
