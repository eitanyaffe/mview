# Gene profile legends

# gene category legend with counts
create_gene_legend <- function(legend_data, field_name) {
  # add labels with counts
  legend_data$label <- paste0(legend_data$category, " (", legend_data$count, ")")
  
  ggplot2::ggplot(legend_data, ggplot2::aes(x = 1, y = seq_len(nrow(legend_data)))) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 0.65, xmax = 1.35, ymin = seq_len(nrow(legend_data)) - 0.45, ymax = seq_len(nrow(legend_data)) + 0.45, fill = color),
      color = "black", size = 0.3
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.4) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = paste("top", field_name, "categories")) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::coord_cartesian(xlim = c(0.5, 4), ylim = c(0.5, nrow(legend_data) + 0.5))
}

# create gene legend from gene data frame
create_gene_legend_from_data <- function(df, profile) {
  # determine display field and color source
  base_field <- profile$color_field
  if (!base_field %in% names(df)) {
    return(NULL)
  }

  # prefer explicit color column if present, else use df$color assigned upstream
  color_field_name <- if (grepl("_color$", base_field)) base_field else paste0(base_field, "_color")
  if (color_field_name %in% names(df)) {
    df$legend_color <- df[[color_field_name]]
  } else if ("color" %in% names(df)) {
    df$legend_color <- df$color
  } else {
    return(NULL)
  }

  # valid entries
  valid_rows <- !is.na(df[[base_field]]) & df[[base_field]] != "" & !is.na(df$legend_color) & df$legend_color != ""
  if (!any(valid_rows)) return(NULL)

  # if all colors are identical (e.g., default gray), skip legend
  unique_colors <- unique(df$legend_color[valid_rows])
  if (length(unique_colors) == 1) return(NULL)

  max_items <- if (!is.null(profile$max_items_in_legend)) profile$max_items_in_legend else 10

  # counts per category
  category_counts <- sort(table(df[[base_field]][valid_rows]), decreasing = TRUE)
  top_categories <- head(category_counts, max_items)
  if (length(top_categories) == 0) return(NULL)

  legend_data <- data.frame(
    category = names(top_categories),
    count = as.numeric(top_categories),
    stringsAsFactors = FALSE
  )

  # choose representative color per category (first occurrence)
  legend_data$color <- sapply(legend_data$category, function(cat) {
    idx <- which(df[[base_field]] == cat & valid_rows)[1]
    if (!is.na(idx)) df$legend_color[idx] else "#cccccc"
  })

  # reverse order so highest counts appear at bottom
  legend_data <- legend_data[nrow(legend_data):1, ]

  legend_gg <- create_gene_legend(legend_data, base_field)
  if (is.null(legend_gg)) return(NULL)

  title <- paste("Top", nrow(legend_data), base_field)
  height <- 60 + nrow(legend_data) * 25
  list(gg = legend_gg, height = height, width = 750, title = title)
}
