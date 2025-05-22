# gene profile plot function

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

  cat(sprintf("plotting %d genes in %s mode\n", nrow(df), mode))

  # Generate colors based on taxonomy
  if (!is.null(df$tax)) {
    # Get unique taxonomies and assign colors
    unique_tax <- sort(unique(df$tax))
    color_palette <- grDevices::rainbow(length(unique_tax))
    tax_colors <- setNames(color_palette, unique_tax)

    # Assign colors to genes
    df$color <- tax_colors[df$tax]
  } else {
    # Default color if taxonomy not available
    df$color <- "steelblue"
  }

  # Create hover text
  df$hover_text <- paste0(
    "Gene: ", df$gene, "\n",
    if (!is.null(df$prot_desc)) paste0("Description: ", df$prot_desc, "\n") else "",
    if (!is.null(df$tax)) paste0("Taxonomy: ", df$tax) else ""
  )

  if (mode == "simple") {
    # Simple mode - efficient rendering for large ranges

    # Sample genes if there are too many
    if (nrow(df) > 5000) {
      set.seed(42) # For consistent sampling
      df <- df[sample(nrow(df), 5000), ]
      cat("sampled to 5000 genes for efficiency\n")
    }

    # Plot genes as simple horizontal lines without strand indicators
    # draw genes as vertical segments with color based on taxonomy
    gg <- gg + ggplot2::geom_segment(
      data = df,
      ggplot2::aes(
        x = gstart, xend = gstart,
        y = 0.45, yend = 0.55,
        color = tax,
        text = hover_text
      ),
      size = 3
    ) +
      # scale_color_discrete maps taxonomy values to distinct colors
      ggplot2::scale_color_discrete(name = "Taxonomy") +
      ggplot2::theme(legend.position = "none")
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
        fill = tax,
        text = hover_text
      ),
      color = NA,
      size = 0.5
    ) +
      ggplot2::scale_fill_discrete(name = "Taxonomy") +
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
