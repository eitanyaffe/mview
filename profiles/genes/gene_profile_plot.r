# gene profile plot function

# helper function to set gene colors based on color field
assign_gene_colors <- function(df, profile) {
  # determine color field name - add _color suffix if not already present
  color_field_name <- profile$color_field
  if (!grepl("_color$", color_field_name)) {
    color_field_name <- paste0(color_field_name, "_color")
  }
  
  # check if color field exists
  if (!color_field_name %in% names(df)) {
    stop(sprintf("error: color field '%s' not found in gene data. available fields: %s", 
                 color_field_name, paste(names(df), collapse = ", ")))
  }
  
  # check if color field contains valid data
  color_values <- df[[color_field_name]]
  if (all(is.na(color_values))) {
    cat(sprintf("warning: color field '%s' contains only NA values - using default light gray color for all genes\n", color_field_name))
    df$color <- "#E8E8E8"
    return(df)
  }
  
  # check if field contains a color mapping instead of individual colors
  if (length(unique(color_values[!is.na(color_values)])) == 1 && 
      is.list(color_values[!is.na(color_values)][1])) {
    stop(sprintf("error: color field '%s' appears to contain a color mapping rather than individual color values. each gene should have its own color value (e.g., '#FF0000'), not a mapping object", color_field_name))
  }
  
  # set default light gray color for all genes
  df$color <- "#E8E8E8"
  
  # use color from specified field
  valid_colors <- !is.na(color_values) & color_values != ""
  if (sum(valid_colors) == 0) {
    cat(sprintf("warning: no valid colors found in field '%s' - all genes will be light gray\n", color_field_name))
  } else {
    df$color[valid_colors] <- color_values[valid_colors]
    cat(sprintf("colored %d out of %d genes using field '%s'\n", 
                sum(valid_colors), nrow(df), color_field_name))
  }
  
  return(df)
}

plot_gene_profile <- function(profile, cxt, genes, gg, mode) {
  if (is.null(genes) || nrow(genes) == 0) {
    return(gg)
  }

  # Process genes for plotting
  df <- genes
  df$contig <- as.character(df$contig)
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)

  # Filter to visible range
  df <- filter_segments(df, cxt, cxt$mapper$xlim)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }

  # clip so that genes are not plotted outside of the visible range
  df$gstart <- pmax(df$gstart, cxt$mapper$xlim[1])
  df$gend <- pmin(df$gend, cxt$mapper$xlim[2])

  # assign colors based on color field
  df <- assign_gene_colors(df, profile)

  cat(sprintf("plotting %d genes in %s mode\n", nrow(df), mode))

  # set hover text from label_field when provided, else use default
  if (profile$show_hover) {
    if (!is.null(profile$label_field) && profile$label_field %in% names(df)) {
      df$hover_text <- as.character(df[[profile$label_field]])
    } else {
      df$hover_text <- paste0(
        "Gene: ", df$gene, "\n",
        if (!is.null(df$prot_desc)) paste0("Description: ", df$prot_desc, "\n") else "",
        if (!is.null(df$tax)) paste0("Taxonomy: ", df$tax) else ""
      )
    }
  } else {
    df$hover_text <- ""
  }

  if (mode == "simple") {
    # Simple mode - efficient rendering for large ranges

    # Sample genes if there are too many
    if (nrow(df) > 5000) {
      set.seed(42) # For consistent sampling
      df <- df[sample(nrow(df), 5000), ]
      cat("sampled to 5000 genes for efficiency\n")
    }

    # Plot genes as simple horizontal lines without strand indicators
    # draw genes as vertical segments with color based on grouping variable
    gg <- gg + ggplot2::geom_segment(
      data = df,
      ggplot2::aes(
        x = gstart, xend = gstart,
        y = 0.45, yend = 0.55,
        text = hover_text
      ),
      color = df$color,
      size = 1
    )
  } else {
    # Full mode - detailed view for zoomed in ranges

    # Plot genes as rectangles
    rect_height <- 0.3
    df$rect_ymin <- 0.5 - rect_height / 2
    df$rect_ymax <- 0.5 + rect_height / 2

    gg <- gg + ggplot2::geom_rect(
      data = df,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = rect_ymin, ymax = rect_ymax,
        fill = I(color),
        text = hover_text
      ),
      size = 0.5
    ) +
      ggplot2::theme(legend.position = "none")

    # Create data for promoter indicators
    promoter_df <- df

    # For + strand, promoter is at start; for - strand, promoter is at end
    promoter_df$x_pos <- ifelse(promoter_df$strand == "+",
      promoter_df$gstart,
      promoter_df$gend
    )

    # Add vertical segments at promoters to indicate strand
    gg <- gg + ggplot2::geom_segment(
      data = promoter_df,
      ggplot2::aes(
        x = x_pos, xend = x_pos,
        y = rect_ymin, yend = rect_ymax,
        text = hover_text
      ),
      color = "black",
      size = 0.5
    )
  }

  return(gg)
}
