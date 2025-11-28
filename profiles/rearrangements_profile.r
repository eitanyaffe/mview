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
    
    # sort rearrangements by clip_out position (left to right)
    sorted_indices <- order(rearrangements_df$gcoord_out)
    sorted_rearrangements <- rearrangements_df[sorted_indices, ]
    
    # tracks the rightmost coordinate at each height level
    height_coord <- c()

    # gap between rearrangements
    gap <- 10
    
    # assign heights to avoid overlaps
    for (i in seq_len(nrow(sorted_rearrangements))) {
      event <- sorted_rearrangements[i, ]
      left_coord <- min(event$gcoord_out, event$gcoord_in)
      right_coord <- max(event$gcoord_out, event$gcoord_in)
      
      # find the lowest available height
      height <- 0
      while (height < length(height_coord)) {
        if (left_coord > (height_coord[height + 1] + gap)) {  # +1 because R is 1-indexed
          break
        }
        height <- height + 1
      }
      
      # if we need a new height level
      if (height >= length(height_coord)) {
        height_coord <- c(height_coord, 0)
      }
      
      # assign the height and update the rightmost position at this height
      sorted_rearrangements$plot_height[i] <- height
      height_coord[height + 1] <- right_coord  # +1 because R is 1-indexed
    }

    # restore original order and add heights back to original dataframe
    rearrangements_df$plot_height <- 0  # initialize
    rearrangements_df$plot_height[sorted_indices] <- sorted_rearrangements$plot_height
    
    return(rearrangements_df)
  }

  plot_f <- function(profile, gg) {
    # get rearrangements data from cache
    rearrangements <- NULL
    if (cache_exists("rearrangements.current")) {
      rearrangements <- cache_get("rearrangements.current")
    }

    if (is.null(rearrangements) || nrow(rearrangements) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    xlim <- cxt_get_xlim()
    
    # filter rearrangements using context coordinates
    # we need to filter by both clip_out and clip_in positions
    # first filter by clip_out
    rearrangements$coord <- rearrangements$out_clip
    filtered_by_out <- cxt_filter_coords(rearrangements)
    
    # then filter by clip_in
    rearrangements$coord <- rearrangements$in_clip
    filtered_by_in <- cxt_filter_coords(rearrangements)
    
    # combine both filters (events that have either clip_out or clip_in in view)
    if (is.null(filtered_by_out) && is.null(filtered_by_in)) {
      return(list(plot = gg, legends = list()))
    }
    
    # merge the results, keeping unique events
    if (!is.null(filtered_by_out) && !is.null(filtered_by_in)) {
      filtered_rearrangements <- rbind(filtered_by_out, filtered_by_in)
      filtered_rearrangements <- filtered_rearrangements[!duplicated(filtered_rearrangements$event_id), ]
    } else if (!is.null(filtered_by_out)) {
      filtered_rearrangements <- filtered_by_out
    } else {
      filtered_rearrangements <- filtered_by_in
    }
    
    if (is.null(filtered_rearrangements) || nrow(filtered_rearrangements) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # add global coordinates for both clip positions
    filtered_rearrangements$gcoord_out <- cxt_contig2global(filtered_rearrangements$contig, filtered_rearrangements$out_clip)
    filtered_rearrangements$gcoord_in <- cxt_contig2global(filtered_rearrangements$contig, filtered_rearrangements$in_clip)

    filtered_rearrangements$gcoord_right <- pmax(filtered_rearrangements$gcoord_out, filtered_rearrangements$gcoord_in)

    # set colors based on rearrangement type
    type_colors <- c(
      "large_insert" = "red",
      "large_delete" = "blue", 
      "large_invert" = "orange"
    )
    filtered_rearrangements$color <- type_colors[filtered_rearrangements$type]
    filtered_rearrangements$color[is.na(filtered_rearrangements$color)] <- "black"  # default for unknown types
    
    # calculate heights to avoid collisions
    filtered_rearrangements <- calculate_rearrangement_heights(filtered_rearrangements)
    
    # check for selected rearrangements for highlighting
    selected_rearrangements <- NULL
    if (cache_exists("rearrangements.selected")) {
      selected_rearrangements <- cache_get("rearrangements.selected")
    }
    
    # mark selected rearrangements for highlighting
    filtered_rearrangements$is_selected <- FALSE
    if (!is.null(selected_rearrangements) && nrow(selected_rearrangements) > 0) {
      filtered_rearrangements$is_selected <- filtered_rearrangements$event_id %in% selected_rearrangements$event_id
    }
    
    # create hover text
    position_label <- paste0(filtered_rearrangements$contig, ":", filtered_rearrangements$out_clip)
    
    hover_text <- paste0(
      "Event: ", filtered_rearrangements$event_id, "<br>",
      "Type: ", filtered_rearrangements$type, "<br>",
      "Position: ", position_label, "<br>",
      "Element: ", filtered_rearrangements$element_contig, ":", 
      filtered_rearrangements$element_start, "-", filtered_rearrangements$element_end, "<br>"
    )
    
    # calculate y positions based on heights
    max_height <- if (nrow(filtered_rearrangements) > 0) max(filtered_rearrangements$plot_height) else 0
    
    # draw rearrangement segments
    # draw non-selected rearrangements first
    non_selected <- filtered_rearrangements[!filtered_rearrangements$is_selected, ]
    selected <- filtered_rearrangements[filtered_rearrangements$is_selected, ]
    
    # add non-selected rearrangements if any exist
    if (nrow(non_selected) > 0) {
      # create hover text for non-selected rearrangements
      non_selected_hover <- hover_text[!filtered_rearrangements$is_selected]
      
      # vertical segment at clip_out position
      gg <- gg + 
        ggplot2::geom_segment(
          data = non_selected,
          ggplot2::aes(
            x = gcoord_out, xend = gcoord_out,
            y = plot_height - 0.4, yend = plot_height + 0.4,
            color = color,
            text = non_selected_hover
          ),
          size = 1
        )
      
      # vertical segment at clip_in position
      gg <- gg + 
        ggplot2::geom_segment(
          data = non_selected,
          ggplot2::aes(
            x = gcoord_in, xend = gcoord_in,
            y = plot_height - 0.4, yend = plot_height + 0.4,
            color = color,
            text = non_selected_hover
          ),
          size = 1
        )
      
      # horizontal segment connecting clip_out to clip_in
      gg <- gg + 
        ggplot2::geom_segment(
          data = non_selected,
          ggplot2::aes(
            x = gcoord_out, xend = gcoord_in,
            y = plot_height, yend = plot_height,
            color = color,
            text = non_selected_hover
          ),
          size = 1
        )
    }
    
    # add selected rearrangements with highlighting if any exist
    if (nrow(selected) > 0) {
      # create hover text for selected rearrangements
      selected_hover <- hover_text[filtered_rearrangements$is_selected]
      
      # vertical segment at clip_out position (highlighted)
      gg <- gg + 
        ggplot2::geom_segment(
          data = selected,
          ggplot2::aes(
            x = gcoord_out, xend = gcoord_out,
            y = plot_height - 0.4, yend = plot_height + 0.4,
            color = color,
            text = selected_hover
          ),
          size = 3,
          alpha = 1
        )
      
      # vertical segment at clip_in position (highlighted)
      gg <- gg + 
        ggplot2::geom_segment(
          data = selected,
          ggplot2::aes(
            x = gcoord_in, xend = gcoord_in,
            y = plot_height - 0.4, yend = plot_height + 0.4,
            color = color,
            text = selected_hover
          ),
          size = 3,
          alpha = 1
        )
      
      # horizontal segment connecting clip_out to clip_in (highlighted)
      gg <- gg + 
        ggplot2::geom_segment(
          data = selected,
          ggplot2::aes(
            x = gcoord_out, xend = gcoord_in,
            y = plot_height, yend = plot_height,
            color = color,
            text = selected_hover
          ),
          size = 3,
          alpha = 1
        )
    }
    
    gg <- gg + ggplot2::scale_color_identity()
    
    # add event ID labels to the right of clip_in lines
    # non-selected labels
    if (nrow(non_selected) > 0) {
      non_selected_hover <- hover_text[!filtered_rearrangements$is_selected]
      
      gg <- gg + 
        ggplot2::geom_text(
          data = non_selected,
          ggplot2::aes(
            x = gcoord_right + (xlim[2] - xlim[1]) * 0.017,
            y = plot_height,
            label = event_id,
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
      selected_hover <- hover_text[filtered_rearrangements$is_selected]
      
      gg <- gg + 
        ggplot2::geom_text(
          data = selected,
          ggplot2::aes(
            x = gcoord_right + (xlim[2] - xlim[1]) * 0.017,
            y = plot_height,
            label = event_id,
            text = selected_hover
          ),
          color = "black",
          size = 4,
          fontface = "bold",
          hjust = 0,
          vjust = 0.5
        )
    }
    
    # set y-axis limits based on the number of height levels
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
