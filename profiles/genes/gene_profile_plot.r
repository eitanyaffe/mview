# gene profile plot function

# Helper function to classify and color genes based on profile settings
classify_and_color_genes <- function(df, profile, unique_tax) {

  # Handle coloring based on color_style
  if (profile$color_style == "by_regex") {
      # Initialize group assignment
    df$group <- NA
      
    # Check each gene against regex patterns
    for (group_name in names(profile$select_groups)) {
      patterns <- profile$select_groups[[group_name]]
        
      for (pattern in patterns) {
        # Find genes matching this pattern that haven't been assigned yet
        matches <- grepl(pattern, df$prot_desc, ignore.case = TRUE) & is.na(df$group)
        df$group[matches] <- group_name
      }
    }
      
    # Filter to only keep genes that matched a pattern
    df <- df[!is.na(df$group), , drop = FALSE]
      
    if (nrow(df) == 0) {
      return(list(df = df, filtered_out = TRUE))
    }

    # Assign colors based on select_colors
    if (is.null(profile$select_colors)) {
        cat("no select_colors provided, skipping rainbow colors\n")
        return(list(df = df, filtered_out = TRUE))
    }
    scolors = unlist(profile$select_colors)
    ix = match(df$group, names(profile$select_groups))
    df$color <- scolors[ix]
    
  } else if (profile$color_style == "by_taxonomy") {
    # Default taxonomy-based coloring
    if (is.null(df$tax)) {
      cat("warning: no tax provided\n")
      return(list(df = df, filtered_out = TRUE))
    }
    ucolors = rainbow(length(unique_tax))
    # Assign colors to genes
    df$color <- ucolors[match(df$tax, unique_tax)]
  } else {
    df$color <- "steelblue"
  }
  
  return(list(df = df, filtered_out = FALSE))
}

plot_gene_profile <- function(profile, cxt, genes, gg, mode) {
  if (is.null(genes) || nrow(genes) == 0) {
    return(gg)
  }

  unique_tax <- cache(paste("gene_profile_unique_tax", cxt$assembly, sep = "_"), {
    sort(unique(genes$tax))
  })

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

  # Classify and color genes
  classification_result <- classify_and_color_genes(df, profile, unique_tax)
  df <- classification_result$df
  
  # Return if all genes were filtered out
  if (classification_result$filtered_out || nrow(df) == 0) {
    return(gg)
  }

  cat(sprintf("plotting %d genes in %s mode\n", nrow(df), mode))

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
    # draw genes as vertical segments with color based on grouping variable
    gg <- gg + ggplot2::geom_segment(
      data = df,
      ggplot2::aes(
        x = gstart, xend = gstart,
        y = 0.45, yend = 0.55,
        color = .data$color,
        text = hover_text
      ),
      size = 0.5
    ) +
      ggplot2::scale_color_discrete(name = "Group") +
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
