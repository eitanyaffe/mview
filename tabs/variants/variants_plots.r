# tabs/variants/variants_plots.r
# standalone variant plotting functions that return ggplot objects

library(ggplot2)
library(scales)

# create variant frequency plot as ggplot object for PDF export
# extracted and adapted from variants_tab.r output$variantFrequencyPlot
create_variant_frequency_plot <- function(variant_data, library_ids, plot_mode = "frequency", variant_span_filter = 0.5) {
  
  # early return for no data - create empty plot
  if (is.null(variant_data) || is.null(variant_data$variants) || nrow(variant_data$variants) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No variants found", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Variant Frequency Plot")
    
    # calculate dimensions: height = 3in, width = 0.5in per library (minimum 2in)
    width_inches <- max(2, length(library_ids) * 0.5)
    height_inches <- 3
    
    return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
  }
  
  # apply span filter
  filtered_data <- filter_variants_by_span(variant_data, variant_span_filter)
  
  if (is.null(filtered_data) || is.null(filtered_data$variants) || nrow(filtered_data$variants) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No variants pass filter", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Variant Frequency Plot")
    
    # calculate dimensions
    width_inches <- max(2, length(library_ids) * 0.5)
    height_inches <- 3
    
    return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
  }
  
  # prepare data for plotting
  support_matrix <- filtered_data$support
  coverage_matrix <- filtered_data$coverage
  variants_df <- filtered_data$variants
  
  # add color information to variants dataframe
  variants_df <- add_variant_colors(variants_df)
  
  # reorder matrices to match configured library order
  lib_indices <- match(library_ids, colnames(support_matrix))
  support_matrix <- support_matrix[, lib_indices, drop = FALSE]
  coverage_matrix <- coverage_matrix[, lib_indices, drop = FALSE]
  
  # calculate data matrix based on plot mode
  if (plot_mode == "variant_support") {
    data_matrix <- support_matrix
    y_label <- "Variant Support"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else if (plot_mode == "variant_coverage") {
    data_matrix <- coverage_matrix
    y_label <- "Variant Coverage"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else { # frequency mode
    data_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
    y_label <- "Frequency"
    y_limits <- c(0, 1)
    y_format <- scales::percent_format()
  }
  
  # convert to long format for plotting
  plot_data <- data.frame()
  for (i in seq_len(nrow(data_matrix))) {
    for (j in seq_len(ncol(data_matrix))) {
      plot_data <- rbind(plot_data, data.frame(
        variant_id = variants_df$variant_id[i],
        variant_type = variants_df$type[i],
        contig = variants_df$contig[i],
        coord = variants_df$coord[i],
        library = library_ids[j],
        value = data_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$value) & !is.na(plot_data$value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Variant Frequency Plot")
    
    # calculate dimensions
    width_inches <- max(2, length(library_ids) * 0.5)
    height_inches <- 3
    
    return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
  }
  
  # use colors from variants dataframe for plot
  plot_data$variant_color <- variants_df$color[match(plot_data$variant_id, variants_df$variant_id)]
  
  # set library order as factor using our configured library_ids order
  plot_data$library <- factor(plot_data$library, levels = library_ids)
  
  # create the plot (no selection highlighting for PDF export)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = variant_color, group = variant_id))
  
  # draw all variants with consistent styling
  p <- p + 
    ggplot2::geom_line(alpha = 0.7, size = 0.8) +
    ggplot2::geom_point(size = 2, alpha = 0.9)
  
  p <- p + ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Library",
      y = y_label,
      title = "Variant Frequency Plot"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  # add y-axis formatting based on plot mode
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = y_format)
  }
  
  # calculate dimensions: height = 3in, width = 0.5in per library (minimum 2in)
  width_inches <- max(2, length(library_ids) * 0.5)
  height_inches <- 3
  
  return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
}
