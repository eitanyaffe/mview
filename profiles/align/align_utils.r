# Shared utilities for alignment and synteny profiles
# Contains mutation color functions, legend functions, and unified mutation plotting

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

# ============================================================================
# COLOR FUNCTIONS
# ============================================================================

# define the shared red color palette used across all mutation visualizations
get_red_palette <- function() {
  return(get_alignment_color_definitions()$red_palette)
}

# shared red color scale function for both continuous and discrete red colors
get_shared_red_scale <- function(continuous = TRUE, num_colors = 6) {
  red_palette <- get_red_palette()
  
  if (continuous) {
    # return the base colors for continuous scale
    return(c(red_palette[1], red_palette[6]))
  } else {
    # return discrete red colors from light to dark
    return(red_palette[seq_len(min(num_colors, length(red_palette)))])
  }
}

# get mutation density colors (shared between bin and full modes)
get_mutation_colors <- function(mutations_per_bp) {
  # use shared red palette
  red_palette <- get_red_palette()
  
  # apply same thresholds as discrete scale using vectorized ifelse for efficiency
  colors <- ifelse(mutations_per_bp <= 0, red_palette[1],
            ifelse(mutations_per_bp <= 0.0001, red_palette[2],
            ifelse(mutations_per_bp <= 0.001, red_palette[3],
            ifelse(mutations_per_bp <= 0.01, red_palette[4],
            ifelse(mutations_per_bp <= 0.1, red_palette[5],
                   red_palette[6])))))
  
  return(colors)
}

# get fixed colors for variant types (shared between full and pileup modes)
get_variant_type_colors <- function(variants) {
  color_defs <- get_alignment_color_definitions()
  
  # colors for gains (3 red shades, darker range, squeezed)
  gain_light <- "#c0c0c0"   # gray color slightly darker than ref for +1
  gain_medium <- "#b30000"  # darker medium red for +2 to +3
  gain_dark <- "#800000"    # darkest red for +4 or more
  
  # colors for losses (3 blue shades, darker range, squeezed)  
  loss_light <- "#c0c0c0"   # gray color slightly darker than ref for -1
  loss_medium <- "#0000b3"  # darker medium blue for -2 to -3
  loss_dark <- "#000080"    # darkest blue for -4 or more
  
  # ref and unknown colors
  ref_color <- "#e0dede"   # light gray
  unknown_color <- "#000000" # black
  
  # assign colors based on variant pattern
  colors <- character(length(variants))
  
  for (i in seq_along(variants)) {
    variant <- variants[i]
    
    if (variant == "REF") {
      colors[i] <- ref_color
    } else if (variant %in% names(color_defs$substitution_colors)) {
      # direct substitution match
      colors[i] <- color_defs$substitution_colors[variant]
    } else if (grepl("^\\+", variant)) {
      # gain pattern: +XXX
      gain_length <- nchar(gsub("^\\+", "", variant))
      colors[i] <- if (gain_length == 1) {
        gain_light
      } else if (gain_length <= 3) {
        gain_medium
      } else {
        gain_dark
      }
    } else if (grepl("^-", variant)) {
      # loss pattern: -XXX
      loss_length <- nchar(gsub("^-", "", variant))
      colors[i] <- if (loss_length == 1) {
        loss_light
      } else if (loss_length <= 3) {
        loss_medium
      } else {
        loss_dark
      }
    } else {
      # unknown variant type
      colors[i] <- unknown_color
    }
  }
  
  return(colors)
}

# get simplified colors for variant types (type mode: sub=orange, deletion=blue, addition=red)
get_mutation_type_colors <- function(variants) {
  # define simplified color scheme
  substitution_color <- "#ff8c00"  # orange for all substitutions
  deletion_color <- "#0066cc"      # blue for all deletions  
  addition_color <- "#cc0000"      # red for all additions
  ref_color <- "#e0dede"           # light gray
  unknown_color <- "#000000"       # black
  
  # assign colors based on variant pattern
  colors <- character(length(variants))
  
  for (i in seq_along(variants)) {
    variant <- variants[i]
    
    if (variant == "REF") {
      colors[i] <- ref_color
    } else if (grepl(":", variant)) {
      # substitution pattern: X:Y
      colors[i] <- substitution_color
    } else if (grepl("^\\+", variant)) {
      # addition pattern: +XXX
      colors[i] <- addition_color
    } else if (grepl("^-", variant)) {
      # deletion pattern: -XXX
      colors[i] <- deletion_color
    } else {
      # unknown variant type
      colors[i] <- unknown_color
    }
  }
  
  return(colors)
}

# ============================================================================
# LEGEND FUNCTIONS
# ============================================================================

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

# ============================================================================
# UNIFIED MUTATION PLOTTING FUNCTION
# ============================================================================

# unified function to plot mutations for both align and synteny profiles
# based on align_profile_full.r mutation plotting logic
plot_mutations_unified <- function(gg, mutation_data, profile) {
  if (is.null(mutation_data) || nrow(mutation_data) == 0) {
    return(gg)
  }
  
  cat(sprintf("plotting %d mutations\n", nrow(mutation_data)))
  
  # sample if too many mutations (using profile's max_mutations if available)
  max_mutations <- profile$max_mutations %||% 10000
  if (nrow(mutation_data) > max_mutations) {
    mutation_data <- mutation_data[sample(nrow(mutation_data), max_mutations), ]
    cat(sprintf("sampled to %d mutations\n", max_mutations))
  }
  
  # choose color function based on mutation_color_mode
  if (profile$mutation_color_mode == "type") {
    mutation_data$fill_color <- get_mutation_type_colors(mutation_data$desc)
  } else {
    # detailed mode (default)
    mutation_data$fill_color <- get_variant_type_colors(mutation_data$desc)
  }
  
  # expect hover_text to be pre-built in the mutation_data
  if (!"hover_text" %in% colnames(mutation_data)) {
    mutation_data$hover_text <- ""
  }
  
  # determine plot style based on zoom level (from align_profile_full logic)
  xlim <- cxt_get_xlim()
  xlim_range <- xlim[2] - xlim[1]
  use_rectangles <- xlim_range < 1000
  
  if (use_rectangles) {
    # plot mutation rectangles when zoomed in
    mutation_border_size <- (profile$full_mutation_lwd %||% 0.5) * 0.8  # slightly thinner than mutation lines
    if (profile$show_hover) {
      gg <- gg + ggplot2::geom_rect(
        data = mutation_data,
        ggplot2::aes(
          xmin = gcoord - 0.5, xmax = gcoord + 0.5,
          ymin = ybottom, ymax = ytop,
          fill = fill_color,
          text = hover_text
        ),
        color = mutation_data$fill_color, size = mutation_border_size
      ) +
        ggplot2::scale_fill_identity()
    } else {
      gg <- gg + ggplot2::geom_rect(
        data = mutation_data,
        ggplot2::aes(
          xmin = gcoord - 0.5, xmax = gcoord + 0.5,
          ymin = ybottom, ymax = ytop,
          fill = fill_color
        ),
        color = mutation_data$fill_color, size = mutation_border_size
      ) +
        ggplot2::scale_fill_identity()
    }
  } else {
    # plot mutation segments when zoomed out
    if (profile$show_hover) {
      gg <- gg + ggplot2::geom_segment(
        data = mutation_data,
        ggplot2::aes(
          x = gcoord, xend = gcoord,
          y = ybottom, yend = ytop,
          color = fill_color,
          text = hover_text
        ),
        size = profile$full_mutation_lwd %||% 0.5
      ) +
        ggplot2::scale_color_identity()
    } else {
      gg <- gg + ggplot2::geom_segment(
        data = mutation_data,
        ggplot2::aes(
          x = gcoord, xend = gcoord,
          y = ybottom, yend = ytop,
          color = fill_color
        ),
        size = profile$full_mutation_lwd %||% 0.5
      ) +
        ggplot2::scale_color_identity()
    }
  }
  
  return(gg)
}
