# default parameters for rearrangements profile
default_rearrangements_params <- list(
  height = list(
    group_id = "rearrangements",
    type = "integer",
    default = 40
  )
)

rearrangements_profile <- function(id, name, height = 100, is_fixed = TRUE,
                                  params = default_rearrangements_params,
                                  auto_register = TRUE) {

  # function to calculate heights for rearrangements to avoid collisions
  calculate_rearrangement_heights <- function(rearrangements_df) {
    if (is.null(rearrangements_df) || nrow(rearrangements_df) == 0) {
      return(rearrangements_df)
    }
    
    # sort rearrangements by left gcoord position
    left_gcoord <- pmin(
      ifelse(!is.na(rearrangements_df$gcoord_out), rearrangements_df$gcoord_out, Inf),
      ifelse(!is.na(rearrangements_df$gcoord_in), rearrangements_df$gcoord_in, Inf))
    sorted_indices <- order(left_gcoord)
    sorted_rearrangements <- rearrangements_df[sorted_indices, ]
    
    height_coord <- c()
    gap <- 10
    
    for (i in seq_len(nrow(sorted_rearrangements))) {
      event <- sorted_rearrangements[i, ]
      coords <- c(event$gcoord_out, event$gcoord_in)
      coords <- coords[!is.na(coords)]
      if (length(coords) == 0) next
      
      left_coord <- min(coords)
      right_coord <- max(coords)
      
      height <- 0
      while (height < length(height_coord)) {
        if (left_coord > (height_coord[height + 1] + gap)) {
          break
        }
        height <- height + 1
      }
      
      if (height >= length(height_coord)) {
        height_coord <- c(height_coord, 0)
      }
      
      sorted_rearrangements$plot_height[i] <- height
      height_coord[height + 1] <- right_coord
    }

    rearrangements_df$plot_height <- 0
    rearrangements_df$plot_height[sorted_indices] <- sorted_rearrangements$plot_height
    
    return(rearrangements_df)
  }

  # helper function to add vertical clip segment
  add_clip_segment <- function(gg, df, gcoord_col, hover_text, line_size) {
    if (nrow(df) == 0) return(gg)
    
    plot_data <- df
    plot_data$gcoord <- plot_data[[gcoord_col]]
    plot_data$hover <- hover_text
    
    gg + ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(
        x = gcoord, xend = gcoord,
        y = plot_height - 0.4, yend = plot_height + 0.4,
        color = color,
        text = hover
      ),
      size = line_size
    )
  }

  # helper function to add horizontal connecting segment
  add_horizontal_segment <- function(gg, df, hover_text, line_size) {
    if (nrow(df) == 0) return(gg)
    
    plot_data <- df
    plot_data$hover <- hover_text
    
    gg + ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(
        x = gcoord_out, xend = gcoord_in,
        y = plot_height, yend = plot_height,
        color = color,
        text = hover
      ),
      size = line_size
    )
  }

  plot_f <- function(profile, gg) {
    rearrangements <- NULL
    if (cache_exists("rearrangements.current")) {
      rearrangements <- cache_get("rearrangements.current")
    }

    if (is.null(rearrangements) || nrow(rearrangements) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    xlim <- cxt_get_xlim()
    
    # check which clips are in view (plotted segment) and in zoom (within xlim)
    in_view_out <- cxt_coords_in_view(rearrangements$contig, rearrangements$out_clip, limit_to_zoom = FALSE)
    in_view_in <- cxt_coords_in_view(rearrangements$contig, rearrangements$in_clip, limit_to_zoom = FALSE)
    in_zoom_out <- cxt_coords_in_view(rearrangements$contig, rearrangements$out_clip, limit_to_zoom = TRUE)
    in_zoom_in <- cxt_coords_in_view(rearrangements$contig, rearrangements$in_clip, limit_to_zoom = TRUE)
    
    # get gcoords for coords that are in view (safe to convert)
    gcoord_out <- rep(NA_real_, nrow(rearrangements))
    gcoord_in <- rep(NA_real_, nrow(rearrangements))
    if (any(in_view_out)) {
      gcoord_out[in_view_out] <- cxt_contig2global(
        rearrangements$contig[in_view_out], 
        rearrangements$out_clip[in_view_out])
    }
    if (any(in_view_in)) {
      gcoord_in[in_view_in] <- cxt_contig2global(
        rearrangements$contig[in_view_in], 
        rearrangements$in_clip[in_view_in])
    }
    
    # check if both in view but crossing xlim (one on each side)
    crosses_xlim <- in_view_out & in_view_in & 
      ((gcoord_out < xlim[1] & gcoord_in > xlim[2]) | (gcoord_in < xlim[1] & gcoord_out > xlim[2]))
    
    # include event if: at least one in zoom, or crosses xlim
    include_event <- in_zoom_out | in_zoom_in | crosses_xlim
    
    if (!any(include_event)) {
      return(list(plot = gg, legends = list()))
    }
    
    filtered <- rearrangements[include_event, ]
    filtered$gcoord_out <- gcoord_out[include_event]
    filtered$gcoord_in <- gcoord_in[include_event]
    
    # determine which clips to show
    filtered$show_out <- in_zoom_out[include_event] | 
      (in_view_out[include_event] & in_zoom_in[include_event]) |
      crosses_xlim[include_event]
    filtered$show_in <- in_zoom_in[include_event] | 
      (in_view_in[include_event] & in_zoom_out[include_event]) |
      crosses_xlim[include_event]

    # set gcoord to NA if not showing
    filtered$gcoord_out[!filtered$show_out] <- NA
    filtered$gcoord_in[!filtered$show_in] <- NA

    # calculate rightmost visible gcoord for label placement
    filtered$gcoord_right <- pmax(
      ifelse(filtered$show_out, filtered$gcoord_out, -Inf),
      ifelse(filtered$show_in, filtered$gcoord_in, -Inf))

    # set colors based on rearrangement type
    type_colors <- c(
      "large_insert" = "red",
      "large_delete" = "blue", 
      "large_invert" = "orange"
    )
    filtered$color <- type_colors[filtered$type]
    filtered$color[is.na(filtered$color)] <- "black"
    
    filtered <- calculate_rearrangement_heights(filtered)
    
    # check for selected rearrangements
    selected_rearrangements <- NULL
    if (cache_exists("rearrangements.selected")) {
      selected_rearrangements <- cache_get("rearrangements.selected")
    }
    
    filtered$is_selected <- FALSE
    if (!is.null(selected_rearrangements) && nrow(selected_rearrangements) > 0) {
      filtered$is_selected <- filtered$event_id %in% selected_rearrangements$event_id
    }
    
    # create hover text
    position_label <- paste0(filtered$contig, ":", filtered$out_clip)
    hover_text <- paste0(
      "Event: ", filtered$event_id, "<br>",
      "Type: ", filtered$type, "<br>",
      "Position: ", position_label, "<br>",
      "Element: ", filtered$element_contig, ":", 
      filtered$element_start, "-", filtered$element_end, "<br>"
    )
    
    max_height <- if (nrow(filtered) > 0) max(filtered$plot_height) else 0
    
    # split into non-selected and selected
    non_selected <- filtered[!filtered$is_selected, ]
    selected <- filtered[filtered$is_selected, ]
    non_selected_hover <- hover_text[!filtered$is_selected]
    selected_hover <- hover_text[filtered$is_selected]
    
    # plot non-selected rearrangements
    if (nrow(non_selected) > 0) {
      # out clips
      show_out_ns <- non_selected[non_selected$show_out, ]
      if (nrow(show_out_ns) > 0) {
        gg <- add_clip_segment(gg, show_out_ns, "gcoord_out", 
          non_selected_hover[non_selected$show_out], 1)
      }
      
      # in clips
      show_in_ns <- non_selected[non_selected$show_in, ]
      if (nrow(show_in_ns) > 0) {
        gg <- add_clip_segment(gg, show_in_ns, "gcoord_in",
          non_selected_hover[non_selected$show_in], 1)
      }
      
      # horizontal segments (only if both clips shown)
      show_both_ns <- non_selected[non_selected$show_out & non_selected$show_in, ]
      if (nrow(show_both_ns) > 0) {
        gg <- add_horizontal_segment(gg, show_both_ns,
          non_selected_hover[non_selected$show_out & non_selected$show_in], 1)
      }
    }
    
    # plot selected rearrangements
    if (nrow(selected) > 0) {
      show_out_sel <- selected[selected$show_out, ]
      if (nrow(show_out_sel) > 0) {
        gg <- add_clip_segment(gg, show_out_sel, "gcoord_out",
          selected_hover[selected$show_out], 3)
      }
      
      show_in_sel <- selected[selected$show_in, ]
      if (nrow(show_in_sel) > 0) {
        gg <- add_clip_segment(gg, show_in_sel, "gcoord_in",
          selected_hover[selected$show_in], 3)
      }
      
      show_both_sel <- selected[selected$show_out & selected$show_in, ]
      if (nrow(show_both_sel) > 0) {
        gg <- add_horizontal_segment(gg, show_both_sel,
          selected_hover[selected$show_out & selected$show_in], 3)
      }
    }
    
    gg <- gg + ggplot2::scale_color_identity()
    
    # add labels
    if (nrow(non_selected) > 0) {
      gg <- gg + ggplot2::geom_text(
        data = non_selected,
        ggplot2::aes(
          x = gcoord_right + (xlim[2] - xlim[1]) * 0.017,
          y = plot_height,
          label = event_id,
          text = non_selected_hover
        ),
        color = "black", size = 3, hjust = 0, vjust = 0.5
      )
    }
    
    if (nrow(selected) > 0) {
      gg <- gg + ggplot2::geom_text(
        data = selected,
        ggplot2::aes(
          x = gcoord_right + (xlim[2] - xlim[1]) * 0.017,
          y = plot_height,
          label = event_id,
          text = selected_hover
        ),
        color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0.5
      )
    }
    
    y_min <- -0.5
    y_max <- max_height + 0.5
    gg <- gg + ggplot2::ylim(y_min, y_max)
    
    return(list(plot = gg, legends = list()))
  }

  # create profile
  profile_create(
    id = id, name = name, type = "rearrangements", height = height, is_fixed = is_fixed,
    attr = list(hide_y_ticks = TRUE),
    params = params, plot_f = plot_f,
    auto_register = auto_register
  )
}
