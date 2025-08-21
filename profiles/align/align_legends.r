# Alignment legends (bin modes)

# SINGLE SOURCE OF TRUTH FOR ALL COLORS
get_alignment_color_definitions <- function() {
  list(
    # red palette for mutation density
    red_palette = c("#fee0d2", "#fc9272", "#fb6a4a", "#cb181d", "#a50f15", "#67000d"),
    
    # 12 nucleotide substitutions (detailed mode)
    substitution_colors = c(
      "A:T" = "#ffcc00", "A:C" = "#00ccff", "A:G" = "#cc00ff",
      "T:A" = "#ff6600", "T:C" = "#00ff66", "T:G" = "#6600ff",
      "C:A" = "#ff0066", "C:T" = "#66ff00", "C:G" = "#ff3399",
      "G:A" = "#ff9900", "G:T" = "#0099ff", "G:C" = "#9900ff"
    ),
    
    # indel colors (detailed mode)
    indel_colors = c(
      "+1 bp" = "#c0c0c0", "+2-3 bp" = "#b30000", "+4+ bp" = "#800000",
      "-1 bp" = "#c0c0c0", "-2-3 bp" = "#0000b3", "-4+ bp" = "#000080",
      "REF" = "#e0dede"
    ),
    
    # simplified variant colors (type mode)
    simplified_colors = c(
      "substitution" = "#ff8c00", "deletion" = "#0066cc", 
      "addition" = "#cc0000", "REF" = "#e0dede"
    ),
    
    # segregating sites gradient
    seg_gradient = c("#c6dbef", "#08306b"),
    
    # non-ref sites gradient  
    nonref_gradient = c("#fff7bc", "#d94701")
  )
}

# safe helpers for backward compatibility
if (!exists("get_red_palette")) {
  get_red_palette <- function() {
    get_alignment_color_definitions()$red_palette
  }
}
if (!exists("get_shared_red_scale")) {
  get_shared_red_scale <- function(continuous = TRUE, num_colors = 6) {
    red_palette <- get_alignment_color_definitions()$red_palette
    if (continuous) c(red_palette[1], red_palette[6]) else red_palette[seq_len(min(num_colors, length(red_palette)))]
  }
}

# mutation density legend (discrete thresholds)
create_mutation_density_legend <- function() {
  colors <- get_alignment_color_definitions()
  labels <- c("0", "≤ 0.01%", "≤ 0.1%", "≤ 1%", "≤ 10%", "> 10%")
  red_palette <- colors$red_palette
  legend_data <- data.frame(
    y = seq_along(labels), x = 1,
    color = red_palette[seq_along(labels)], label = labels,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35, ymin = y - 0.45, ymax = y + 0.45, fill = color),
      color = "black", size = 0.4
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.6) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "mutations per bp") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(labels) + 0.5))
}

# fixed range gradient legend (0 to max_val)
create_gradient_legend <- function(title, colors, max_val, n_steps = 10, as_percent = FALSE) {
  cols <- grDevices::colorRampPalette(colors)(n_steps)
  values <- seq(0, max_val, length.out = n_steps)
  legend_data <- data.frame(
    y = seq_len(n_steps), x = 1, color = cols, value = values,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35, ymin = y - 0.45, ymax = y + 0.45, fill = color),
      color = "black", size = 0.3
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::annotate(
      "text", x = 1.7, y = 1,
      label = if (as_percent) sprintf("%.2f%%", 0.0) else sprintf("%.1f", 0.0),
      hjust = 0, vjust = 0.5, size = 3.2
    ) +
    ggplot2::annotate(
      "text", x = 1.7, y = n_steps,
      label = if (as_percent) sprintf("%.2f%%", max_val * 100) else sprintf("%.1f", max_val),
      hjust = 0, vjust = 0.5, size = 3.2
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, n_steps + 0.5))
}

# stacked mutation rates categories legend (discrete)
create_stacked_mutation_rates_legend <- function() {
  colors <- get_alignment_color_definitions()
  red_colors <- colors$red_palette
  categories <- c(">10%", "1%-10%", "0.1%-1%", "0.01%-0.1%", "0.001%-0.01%", "0")
  colors <- c(red_colors[6], red_colors[5], red_colors[4], red_colors[3], red_colors[2], red_colors[1])
  legend_data <- data.frame(
    y = seq_along(categories), x = 1,
    color = colors, label = categories,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35, ymin = y - 0.45, ymax = y + 0.45, fill = color),
      color = "black", size = 0.3
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.4) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "mutation rate bins") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(categories) + 0.5))
}

# detailed variant colors legend (12 substitutions + indels)
create_detailed_variant_legend <- function() {
  colors <- get_alignment_color_definitions()
  all_colors <- c(colors$substitution_colors, colors$indel_colors)
  all_labels <- names(all_colors)
  
  legend_data <- data.frame(
    y = seq_along(all_labels), x = 1,
    color = all_colors, label = all_labels,
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35, ymin = y - 0.45, ymax = y + 0.45, fill = color),
      color = "black", size = 0.3
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.2) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "variant types (detailed)") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0.5, length(all_labels) + 0.5))
}

# simplified variant colors legend (3 types)
create_simplified_variant_legend <- function() {
  color_defs <- get_alignment_color_definitions()
  colors <- color_defs$simplified_colors
  labels <- names(colors)
  
  legend_data <- data.frame(
    y = seq_along(labels), x = 1,
    color = colors, label = labels,
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35, ymin = y - 0.45, ymax = y + 0.45, fill = color),
      color = "black", size = 0.4
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.6) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "variant types (simplified)") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(labels) + 0.5))
}


