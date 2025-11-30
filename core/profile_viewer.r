# core/profile_viewer.r

# This file handles the plotting logic for all registered profiles.
# For each profile, it calls pre_plot, profile$plot_f, and post_plot, and collects the results.

library(plotly)
library(patchwork)

#' Pre-plot function: sets up the base ggplot object for a profile
#' @param profile The profile object
#' @return A ggplot object with base layers (grid, y-label, etc.)
pre_plot <- function(profile) {
  gg <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = cxt_get_xlim()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )

  # Only add y-label if not hidden
  if (!isTRUE(profile$attr$hide_y_label)) {
    gg <- gg + ggplot2::labs(y = profile$attr$title)
  } else {
    gg <- gg + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  }

  # Hide y-ticks if specified
  if (isTRUE(profile$attr$hide_y_ticks)) {
    gg <- gg + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
  } else {
    # Show the y-axis line when showing ticks
    gg <- gg + ggplot2::theme(
      axis.line.y = ggplot2::element_line(color = "black", size = 0.5)
    )
  }

  if (isTRUE(profile$attr$should_plot_grid)) {
    gg <- gg + ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(),
      panel.grid.major.y = ggplot2::element_line()
    )
  }
  
  gg
}

#' Post-plot function: can add overlays (e.g. vlines) on top of the profile plot
#' @param gg The ggplot object after profile$plot_f
#' @param profile The profile object
#' @return The ggplot object with overlays added
post_plot <- function(gg, profile) {
  # Add vertical lines between contigs (at contig ends, not starts)
  cdf <- cxt_get_entire_view()
  
  if (!is.null(cdf) && nrow(cdf) > 1) {
    # Only draw lines at the END of each contig (except the last one)
    contig_ends <- cdf$vend[-nrow(cdf)]  # all except last
    gg <- gg + ggplot2::geom_vline(xintercept = contig_ends, color = "gray")
  }
  gg
}

#' Update cached ggplot objects for all registered profiles
#' This function does the expensive ggplot creation and ggplotly conversion
#' @param profile_plotly_objects Reactive value to store results
update_gg_objects <- function(profile_plotly_objects) {
  profiles <- profiles_get_all()
  
  if (length(profiles) == 0) {
    profile_plotly_objects(NULL)
    return()
  }

  plotly_list <- list()
  legends <- list()

  for (id in names(profiles)) {
    profile <- profiles[[id]]

    # apply default height fields if missing
    if (is.null(profile$is_fixed)) profile$is_fixed <- FALSE
    if (is.null(profile$height)) profile$height <- 150

    if (length(profile$params) > 0) {
      for (i in seq_along(profile$params)) {
        param_id <- names(profile$params)[i]
        group_id <- profile$params[[i]]$group
        if (is.null(group_id)) {
          next
        }
        param_value <- get_param(group_id, param_id)()
        profile[[param_id]] <- param_value
      }
    }
    
    # update height from parameters if it exists
    if ("height" %in% names(profile$params)) {
      height_param <- profile$params[["height"]]
      if (!is.null(height_param) && !is.null(height_param$group)) {
        height_value <- get_param(height_param$group, "height")()
        profile$height <- height_value
      }
    }

    # Create base ggplot
    gg <- pre_plot(profile)

    # Get the profile's plot result
    pf_result <- profile$plot_f(profile, gg)

    # profiles must return list(plot = gg, legends = list_of_legends)
    if (is.list(pf_result)) {
      gg_plot <- pf_result$plot
      profile_legends <- pf_result$legends
      if (!is.null(profile_legends) && length(profile_legends) > 0) {
        legends[[id]] <- profile_legends
      }
    } else {
      gg_plot <- pf_result
    }

    # validate that we have a plot object (not NULL)
    if (is.null(gg_plot)) {
      warning(sprintf("Profile '%s' returned NULL plot object, skipping", id))
      next
    }

    # Add overlays
    gg_final <- post_plot( gg_plot, profile)

    # Convert to plotly, using 'text' aesthetic when hover is enabled
    if (isFALSE(profile$show_hover)) {
      p_ly <- plotly::ggplotly(gg_final, tooltip = NULL) %>%
         plotly::style(hoverinfo = "none")
    } else {
      p_ly <- plotly::ggplotly(gg_final, tooltip = "text")
    }

    # Store in list with profile metadata
    plotly_list[[id]] <- list(
      plotly = p_ly,
      profile = profile,
      is_fixed = profile$is_fixed,
      height = profile$height
    )
  }

  if (length(plotly_list) == 0) {
    profile_plotly_objects(NULL)
    return()
  }

  # Store plotly objects and legends
  result <- list(
    plotly_list = plotly_list,
    legends = legends
  )
  
  profile_plotly_objects(result)
}

#' Calculate height fractions for profiles based on container height
#' @param plotly_list List of plotly objects with profile metadata
#' @param container_height_px Container height in pixels
#' @return List with height fractions and total height
calculate_height_fractions <- function(plotly_list, container_height_px) {
  if (length(plotly_list) == 0) {
    return(numeric())
  }

  # Extract profile height info
  fixed_profiles <- sapply(plotly_list, function(p) p$is_fixed)
  profile_heights <- sapply(plotly_list, function(p) p$height)
  profile_names <- names(plotly_list)
  
  # Calculate minimum container height needed
  fixed_height_total <- sum(profile_heights[fixed_profiles])
  
  # If no container height specified, use profile heights
  if (is.null(container_height_px)) {
    container_height_px <- as.integer(min(1500, sum(profile_heights)))  # max 1500px
  }
  
  # Calculate height distribution
  remaining_height <- container_height_px - fixed_height_total
  dynamic_profiles <- !fixed_profiles
  n_dynamic <- sum(dynamic_profiles)
  
  final_heights <- numeric(length(plotly_list))
  
  # Fixed profiles keep their height
  final_heights[fixed_profiles] <- profile_heights[fixed_profiles]
  
  # Dynamic profiles share remaining space proportionally based on their heights
  if (n_dynamic > 0) {
    if (remaining_height > 0) {
      # Calculate total height of dynamic profiles for proportional distribution
      dynamic_indices <- which(dynamic_profiles)
      dynamic_heights <- profile_heights[dynamic_indices]
      total_dynamic_height <- sum(dynamic_heights)
      
      if (total_dynamic_height > 0) {
        # Distribute remaining height proportionally based on original heights
        for (i in dynamic_indices) {
          proportion <- profile_heights[i] / total_dynamic_height
          allocated_height <- remaining_height * proportion
          final_heights[i] <- allocated_height
        }
      } else {
        # Fallback to equal distribution if all dynamic heights are 0
        height_per_dynamic <- remaining_height / n_dynamic
        for (i in dynamic_indices) {
          final_heights[i] <- height_per_dynamic
        }
      }
    } else {
      # If no remaining height, give dynamic profiles zero height
      for (i in which(dynamic_profiles)) {
        final_heights[i] <- 0
      }
    }
  }
  
  # Calculate fractions directly from final heights
  total_height <- sum(final_heights)
  
  # Safety check: ensure total_height > 0
  if (total_height <= 0) {
    # Fallback: equal fractions
    fractions <- rep(1/length(plotly_list), length(plotly_list))
  } else {
    fractions <- final_heights / total_height
    # Normalize fractions to ensure they sum to exactly 1.0 (avoid floating point errors)
    fractions <- fractions / sum(fractions)
  }

  fractions
}

#' Plot all registered profiles using cached plotly objects
#' @param plotly_objects_result Result from update_gg_objects containing plotly_list and legends
#' @param container_height_px Container height in pixels (NULL = calculate from profile heights)
#' @return A list with plotly object and total height: list(plot = plotly_obj)
plot_profiles_cached <- function(plotly_objects_result, container_height_px = NULL) {
  if (is.null(plotly_objects_result) || is.null(plotly_objects_result$plotly_list)) {
    return(NULL)
  }
  
  plotly_list <- plotly_objects_result$plotly_list
  legends <- plotly_objects_result$legends
  
  # Calculate height fractions
  fractions <- calculate_height_fractions(plotly_list, container_height_px)
  if (length(fractions) == 0) {
    return(NULL)
  }
  
  # Extract just the plotly objects for subplot
  plotly_objects <- lapply(plotly_list, function(p) p$plotly)
  

  # Combine plots
  if (length(plotly_objects) > 1) {
    combined_plot <- plotly::subplot(
      plotly_objects,
      nrows = length(plotly_objects),
      shareX = TRUE,
      titleY = FALSE,
      margin = 0.01,
      heights = fractions
    )
  } else if (length(plotly_objects) == 1) {
    combined_plot <- plotly_objects[[1]]
  } else {
    return(NULL)
  }

  cat(sprintf("configuring combined plot\n"))

  # Configure the combined plot
  layout_args <- list(
    xaxis = list(title = "", range = cxt_get_xlim()),
    margin = list(l = 150, r = 20, t = 30, b = 30)
  )
  
  # Build horizontal y-labels using plotly annotations aligned to each subplot row
  profiles <- profiles_get_all()
  annotations <- list()
  
  # collect coordinate annotations via profile hooks
  cat("collecting coordinate annotations\n")
  for (i in seq_along(profiles)) {
    profile <- profiles[[i]]
    get_ann_f <- profile$get_annotations
    if (is.null(get_ann_f) || !is.function(get_ann_f)) next

    coord_data <- get_ann_f(profile)
    if (is.null(coord_data) || nrow(coord_data) == 0) next

    # target this subplot's y-axis for annotation placement
    y_axis_ref <- if (i == 1) "y" else paste0("y", i)
    y_val <- 0.52  # place below ticks (ticks at ~0.70)

    cat("profile ", profile$id, " adding ", nrow(coord_data), " annotations\n")
    for (j in seq_len(nrow(coord_data))) {
      annotations[[length(annotations) + 1]] <- list(
        x = coord_data$x[j], xref = "x", xanchor = "center",
        y = y_val, yref = y_axis_ref, yanchor = "middle",
        text = coord_data$label[j], textangle = 90,
        showarrow = FALSE, font = list(size = 10)
      )
    }
  }
  cat("adding y-axis labels\n")
  current_layout <- combined_plot$x$layout
  for (i in seq_along(profiles)) {
    profile <- profiles[[i]]
    # Always clear default vertical y-axis title text
    yaxis_name <- if (i == 1) "yaxis" else paste0("yaxis", i)
    layout_args[[yaxis_name]] <- list(title = list(text = ""))

    # Add horizontal annotation for y-label if not hidden
    if (!isTRUE(profile$attr$hide_y_label) && !is.null(profile$attr$title)) {
      axis_layout <- current_layout[[yaxis_name]]
      domain <- NULL
      if (!is.null(axis_layout) && !is.null(axis_layout$domain)) {
        domain <- axis_layout$domain
      }
      if (!is.null(domain) && length(domain) == 2) {
        y_mid <- mean(domain)
        annotations[[length(annotations) + 1]] <- list(
          x = 0.01, xref = "paper", xanchor = "right", xshift = -30,
          y = y_mid, yref = "paper", yanchor = "middle",
          text = profile$attr$title, textangle = 0,
          showarrow = FALSE, align = "right", font = list(size = 12)
        )
      }
    }
  }
  if (length(annotations) > 0) {
    layout_args$annotations <- annotations
    cat("adding annotations, length: ", length(annotations), "\n")
  }
  
  combined_plot <- do.call(plotly::layout, c(list(combined_plot), layout_args))
  return(list(plot = combined_plot, legends = legends))
}

#' Plot all registered profiles
#' @return A list with plotly object and total height in pixels: list(plot = plotly_obj)
plot_profiles <- function() {
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }

  plotly_list <- list()
  heights <- numeric()
  legends <- list()

  # clear alignment object list before processing profiles
  cache_set("aln_obj", list())

  for (id in names(profiles)) {
    cat(sprintf("plotting profile: %s\n", id))
    profile <- profiles[[id]]

    # apply default height fields if missing
    if (is.null(profile$is_fixed)) profile$is_fixed <- FALSE
    if (is.null(profile$height)) profile$height <- 150

    if (length(profile$params) > 0) {
      for (i in seq_along(profile$params)) {
        param_id <- names(profile$params)[i]
        group_id <- profile$params[[i]]$group
        if (is.null(group_id)) {
          next
        }
        param_value <- get_param(group_id, param_id)()
        profile[[param_id]] <- param_value
      }
    }
    
    # update height from parameters if it exists
    if ("height" %in% names(profile$params)) {
      height_param <- profile$params[["height"]]
      if (!is.null(height_param) && !is.null(height_param$group)) {
        height_value <- get_param(height_param$group, "height")()
        profile$height <- height_value
      }
    }

    # Create base ggplot
    gg <- pre_plot( profile)

    # Get the profile's plot result
    pf_result <- profile$plot_f(profile, gg)

    # profiles must return list(plot = gg, legends = list_of_legends)
    if (is.list(pf_result)) {
      gg_plot <- pf_result$plot
      profile_legends <- pf_result$legends
      if (!is.null(profile_legends) && length(profile_legends) > 0) {
        legends[[id]] <- profile_legends
      }
    } else {
      gg_plot <- pf_result
    }

    # validate that we have a plot object (not NULL)
    if (is.null(gg_plot)) {
      warning(sprintf("Profile '%s' returned NULL plot object, skipping", id))
      next
    }

    # Add overlays
    gg_final <- post_plot( gg_plot, profile)

    # Convert to plotly, using 'text' aesthetic when hover is enabled
    if (isFALSE(profile$show_hover)) {
      p_ly <- plotly::ggplotly(gg_final, tooltip = NULL) %>%
         plotly::style(hoverinfo = "none")
    } else {
      p_ly <- plotly::ggplotly(gg_final, tooltip = "text")
    }

    # Store in list
    plotly_list[[id]] <- p_ly
    heights <- c(heights, profile$height)
  }

  if (length(plotly_list) == 0) {
    return(NULL) # Return NULL if no plots were generated
  }

  # Heights are now in pixels - calculate proportions for plotly subplot
  total_height <- as.integer(sum(heights))
  if (total_height > 0) {
    heights <- heights / total_height
  } else {
    # This should not happen with validation, but fallback
    heights <- rep(1 / length(heights), length(heights))
  }

  cat(sprintf("combining %d profiles using plotly\n", length(plotly_list)))

  # Combine plots
  # Need to handle the case where plotly_list is empty or has only one plot
  if (length(plotly_list) > 1) {
    combined_plot <- plotly::subplot(
      plotly_list,
      nrows = length(plotly_list),
      shareX = TRUE,
      titleY = FALSE,
      margin = 0.01,
      heights = if (length(heights) == length(plotly_list)) heights else NULL
    )
  } else if (length(plotly_list) == 1) {
    combined_plot <- plotly_list[[1]]
  } else {
    return(NULL) # Should not happen if profiles list was not empty initially
  }

  cat(sprintf("configuring combined plot\n"))

  # Configure the combined plot
  layout_args <- list(
    xaxis = list(title = "", range = cxt_get_xlim()),
    margin = list(l = 150, r = 20, t = 30, b = 30)
  )
  
  # Build horizontal y-labels using plotly annotations aligned to each subplot row
  profiles <- profiles_get_all()
  annotations <- list()
  
  # collect coordinate annotations via profile hooks
  cat("collecting coordinate annotations\n")
  for (i in seq_along(profiles)) {
    profile <- profiles[[i]]
    get_ann_f <- profile$get_annotations
    if (is.null(get_ann_f) || !is.function(get_ann_f)) next

    coord_data <- get_ann_f(profile)
    if (is.null(coord_data) || nrow(coord_data) == 0) next

    # target this subplot's y-axis for annotation placement
    y_axis_ref <- if (i == 1) "y" else paste0("y", i)
    y_val <- 0.52  # place below ticks (ticks at ~0.70)

    cat("profile ", profile$id, " adding ", nrow(coord_data), " annotations\n")
    for (j in seq_len(nrow(coord_data))) {
      annotations[[length(annotations) + 1]] <- list(
        x = coord_data$x[j], xref = "x", xanchor = "center",
        y = y_val, yref = y_axis_ref, yanchor = "middle",
        text = coord_data$label[j], textangle = 90,
        showarrow = FALSE, font = list(size = 10)
      )
    }
  }
  cat("adding y-axis labels\n")
  current_layout <- combined_plot$x$layout
  for (i in seq_along(profiles)) {
    profile <- profiles[[i]]
    # Always clear default vertical y-axis title text
    yaxis_name <- if (i == 1) "yaxis" else paste0("yaxis", i)
    layout_args[[yaxis_name]] <- list(title = list(text = ""))

    # Add horizontal annotation for y-label if not hidden
    if (!isTRUE(profile$attr$hide_y_label) && !is.null(profile$attr$title)) {
      axis_layout <- current_layout[[yaxis_name]]
      domain <- NULL
      if (!is.null(axis_layout) && !is.null(axis_layout$domain)) {
        domain <- axis_layout$domain
      }
      if (!is.null(domain) && length(domain) == 2) {
        y_mid <- mean(domain)
        annotations[[length(annotations) + 1]] <- list(
          x = 0.01, xref = "paper", xanchor = "right", xshift = -30,
          y = y_mid, yref = "paper", yanchor = "middle",
          text = profile$attr$title, textangle = 0,
          showarrow = FALSE, align = "right", font = list(size = 12)
        )
      }
    }
  }
  if (length(annotations) > 0) {
    layout_args$annotations <- annotations
    cat("adding annotations, length: ", length(annotations), "\n")
  }
  
  combined_plot <- do.call(plotly::layout, c(list(combined_plot), layout_args))
  return(list(plot = combined_plot, legends = legends))
}

#' Plot all registered profiles as ggplot objects for PDF export
#' @return A list with combined ggplot object and total height: list(plot = ggplot_obj, total_height = numeric)
plot_profiles_ggplot <- function() {
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }

  ggplot_list <- list()
  heights <- numeric()
  legends <- list()

  for (id in names(profiles)) {
    cat(sprintf("plotting profile for PDF: %s\n", id))
    profile <- profiles[[id]]

    # apply default height fields if missing
    if (is.null(profile$is_fixed)) profile$is_fixed <- FALSE
    if (is.null(profile$height)) profile$height <- 150

    if (length(profile$params) > 0) {
      for (i in seq_along(profile$params)) {
        param_id <- names(profile$params)[i]
        group_id <- profile$params[[i]]$group
        if (is.null(group_id)) {
          next
        }
        param_value <- get_param(group_id, param_id)()
        profile[[param_id]] <- param_value
      }
    }
    
    # update height from parameters if it exists
    if ("height" %in% names(profile$params)) {
      height_param <- profile$params[["height"]]
      if (!is.null(height_param) && !is.null(height_param$group)) {
        height_value <- get_param(height_param$group, "height")()
        profile$height <- height_value
      }
    }

    # Create base ggplot
    gg <- pre_plot( profile)

    # Get the profile's plot result
    pf_result <- profile$plot_f(profile, gg)

    # profiles must return list(plot = gg, legends = list_of_legends)
    if (is.list(pf_result)) {
      gg_plot <- pf_result$plot
      profile_legends <- pf_result$legends
      if (!is.null(profile_legends) && length(profile_legends) > 0) {
        legends[[id]] <- profile_legends
      }
    } else {
      gg_plot <- pf_result
    }

    # validate that we have a plot object (not NULL)
    if (is.null(gg_plot)) {
      warning(sprintf("Profile '%s' returned NULL plot object, skipping", id))
      next
    }

    # Add overlays
    gg_final <- post_plot( gg_plot, profile)
    
    # Apply styling for PDF export (ensure it's not overridden)
    if (!isTRUE(profile$attr$hide_y_label)) {
      gg_final <- gg_final + 
        ggplot2::labs(y = profile$attr$title) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text(angle = 0, hjust = 1, vjust = 0.5),
          axis.text.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    } else {
      gg_final <- gg_final + 
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }

    # Store the ggplot object directly (no plotly conversion)
    ggplot_list[[id]] <- gg_final
    heights <- c(heights, profile$height)
  }

  if (length(ggplot_list) == 0) {
    return(NULL)
  }

  # Calculate height ratios for patchwork
  total_height <- as.integer(sum(heights))
  if (total_height > 0) {
    height_ratios <- heights / total_height
  } else {
    height_ratios <- rep(1 / length(heights), length(heights))
  }

  cat(sprintf("combining %d profiles using patchwork for PDF\n", length(ggplot_list)))

  # Combine plots using patchwork
  if (length(ggplot_list) > 1) {
    # Use patchwork to stack plots vertically with height ratios
    combined_plot <- Reduce(`/`, ggplot_list) + 
      patchwork::plot_layout(heights = height_ratios)
  } else {
    combined_plot <- ggplot_list[[1]]
  }
  
  # Add coordinate annotations to the bottom plot (similar to plotly version)
  if (length(ggplot_list) > 0) {
    profiles <- profiles_get_all()
    coord_annotations <- list()
    
    # collect coordinate annotations from profiles
    for (i in seq_along(profiles)) {
      profile <- profiles[[i]]
      get_ann_f <- profile$get_annotations
      if (is.null(get_ann_f) || !is.function(get_ann_f)) next
      
      coord_data <- get_ann_f(profile)
      if (is.null(coord_data) || nrow(coord_data) == 0) next
      
      coord_annotations <- rbind(coord_annotations, coord_data)
    }
    
    # Add coordinate annotations to the last (bottom) plot if we have any
    if (length(coord_annotations) > 0 && nrow(coord_annotations) > 0) {
      bottom_plot_name <- names(ggplot_list)[length(ggplot_list)]
      bottom_plot <- ggplot_list[[bottom_plot_name]]
      
      # Add coordinate labels as text annotations
      bottom_plot <- bottom_plot + 
        ggplot2::annotate("text", 
                         x = coord_annotations$x, 
                         y = -Inf, 
                         label = coord_annotations$label,
                         angle = 90, 
                         hjust = 0, 
                         vjust = 0.5,
                         size = 3)
      
      # Update the plot in the list
      ggplot_list[[bottom_plot_name]] <- bottom_plot
      
      # Recreate the combined plot with annotations
      if (length(ggplot_list) > 1) {
        combined_plot <- Reduce(`/`, ggplot_list) + 
          patchwork::plot_layout(heights = height_ratios)
      } else {
        combined_plot <- ggplot_list[[1]]
      }
    }
  }

  return(list(plot = combined_plot, total_height = total_height, legends = legends))
}
