# default parameters for variants profile
default_variants_params <- list(
  height = list(
    group_id = "variants",
    type = "integer",
    default = 40
  )
)

variants_profile <- function(id, name, height = 60, is_fixed = TRUE,
                            params = default_variants_params,
                            auto_register = TRUE) {

  plot_f <- function(profile, gg) {
    # get variants data from cache
    variants <- NULL
    if (cache_exists("variants.current")) {
      variants <- cache_get("variants.current")
    }

    if (is.null(variants) || nrow(variants) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # filter variants using context coordinates (similar to filter_coords)
    filtered_variants <- cxt_filter_coords(variants)
    if (is.null(filtered_variants) || nrow(filtered_variants) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # ensure color column exists
    if (!"color" %in% colnames(filtered_variants)) {
      filtered_variants$color <- "#000000"  # default black
    }
    
    # check for selected variants for highlighting
    # try to get current table selection directly from cache first
    selected_variants <- NULL
    if (cache_exists("variants.selected")) {
      selected_variants <- cache_get("variants.selected")
    }
    
    # mark selected variants for highlighting
    filtered_variants$is_selected <- FALSE
    if (!is.null(selected_variants) && nrow(selected_variants) > 0) {
      filtered_variants$is_selected <- filtered_variants$variant_id %in% selected_variants$variant_id
    }
    
    # create hover text (same format as frequency plot)
    position_label <- paste0(filtered_variants$contig, ":", filtered_variants$coord)
    
    # truncate descriptions for hover text (same logic as frequency plot)
    description_display <- sapply(filtered_variants$desc, function(desc) {
      if (is.na(desc) || nchar(desc) <= 15) {
        return(desc)
      } else {
        # truncate middle for hover (keep first and last 5 characters)
        first_5 <- substr(desc, 1, 5)
        last_5 <- substr(desc, nchar(desc) - 4, nchar(desc))
        return(paste0(first_5, "...", last_5))
      }
    })
    
    hover_text <- paste0(
      "Variant: ", filtered_variants$variant_id, "<br>",
      "Position: ", position_label, "<br>",
      "Description: ", description_display, "<br>",
      ifelse(!is.na(filtered_variants$mutation_desc) & filtered_variants$mutation_desc != "", 
             paste0("AA Change: ", filtered_variants$mutation_desc, "<br>"), ""),
      "Type: ", filtered_variants$type, "<br>"
    )
    
    # add vertical lines for variants spanning the full profile height
    # draw non-selected variants first
    non_selected <- filtered_variants[!filtered_variants$is_selected, ]
    selected <- filtered_variants[filtered_variants$is_selected, ]
    
    # add non-selected variants if any exist
    if (nrow(non_selected) > 0) {
      # create hover text for non-selected variants
      non_selected_hover <- hover_text[!filtered_variants$is_selected]
      
      gg <- gg + 
        ggplot2::geom_segment(
          data = non_selected,
          ggplot2::aes(
            x = gcoord, xend = gcoord,
            y = -0.5, yend = 0.5,
            color = color,
            text = non_selected_hover
          ),
          size = 1
        )
    }
    
    # add selected variant with highlighting if any exist
    if (nrow(selected) > 0) {
      # create hover text for selected variants
      selected_hover <- hover_text[filtered_variants$is_selected]
      
      gg <- gg + 
        ggplot2::geom_segment(
          data = selected,
          ggplot2::aes(
            x = gcoord, xend = gcoord,
            y = -0.5, yend = 0.5,
            color = color,
            text = selected_hover
          ),
          size = 3,
          alpha = 1
        )
    }
    
    gg <- gg + ggplot2::scale_color_identity()
    
    # add variant ID labels to the right of lines
    xlim <- cxt_get_xlim()
    label_offset <- (xlim[2] - xlim[1]) * 0.017
    
    # non-selected labels
    if (nrow(non_selected) > 0) {
      non_selected_hover <- hover_text[!filtered_variants$is_selected]
      
      gg <- gg + 
        ggplot2::geom_text(
          data = non_selected,
          ggplot2::aes(
            x = gcoord + label_offset,
            y = -0.05,
            label = variant_id,
            text = non_selected_hover
          ),
          color = "black",
          size = 3,
          hjust = 0,
          vjust = 0.5
        )
    }
    
    # selected labels with bold styling
    if (nrow(selected) > 0) {
      selected_hover <- hover_text[filtered_variants$is_selected]
      
      gg <- gg + 
        ggplot2::geom_text(
          data = selected,
          ggplot2::aes(
            x = gcoord + label_offset,
            y = -0.05,
            label = variant_id,
            text = selected_hover
          ),
          color = "black",
          size = 4,
          fontface = "bold",
          hjust = 0,
          vjust = 0.5
        )
    }
    
    gg <- gg +
      ggplot2::ylim(-0.5, 0.5)
    
    return(list(plot = gg, legends = list()))
  }

  # create profile
  profile_create(
    id = id, name = name, type = "variants", height = height, is_fixed = is_fixed,
    attr = list(hide_y_ticks = TRUE),
    params = params, plot_f = plot_f,
    auto_register = auto_register
  )
}
