# core/profile_viewer.r

# This file handles the plotting logic for all registered profiles.
# For each profile, it calls pre_plot, profile$plot_f, and post_plot, and collects the results.

library(plotly)

#' Pre-plot function: sets up the base ggplot object for a profile
#' @param cxt Context object
#' @param profile The profile object
#' @return A ggplot object with base layers (grid, y-label, etc.)
pre_plot <- function(cxt, profile) {
  # Minimal base plot with xlim from cxt$mapper, y-label from profile$attr$title
  gg <- ggplot2::ggplot() +
    ggplot2::xlim(cxt$mapper$xlim) +
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
#' @param cxt Context object
#' @param gg The ggplot object after profile$plot_f
#' @param profile The profile object
#' @return The ggplot object with overlays added
post_plot <- function(cxt, gg, profile) {
  # Add vertical lines between contigs (at contig ends, not starts)
  if (!is.null(cxt$mapper$cdf) && nrow(cxt$mapper$cdf) > 1) {
    # Only draw lines at the END of each contig (except the last one)
    contig_ends <- cxt$mapper$cdf$end[-nrow(cxt$mapper$cdf)]  # all except last
    gg <- gg + ggplot2::geom_vline(xintercept = contig_ends, color = "gray")
  }
  gg
}

#' Plot all registered profiles
#' @param cxt The context object containing mapper and other dynamic info
#' @return A list with plotly object and total height in pixels: list(plot = plotly_obj, total_height = numeric)
plot_profiles <- function(cxt) {
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }

  plotly_list <- list()
  heights <- numeric()

  for (id in names(profiles)) {
    cat(sprintf("plotting profile: %s\n", id))
    profile <- profiles[[id]]

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
    gg <- pre_plot(cxt, profile)

    # Get the profile's plot result
    pf_result <- profile$plot_f(profile, cxt, gg)

    # Add overlays
    gg_final <- post_plot(cxt, pf_result, profile)

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
  total_height <- sum(heights)
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
    xaxis = list(title = ""),
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

    coord_data <- get_ann_f(profile, cxt)
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
  cat(sprintf("plotting done, total height: %dpx\n", total_height))
  return(list(plot = combined_plot, total_height = total_height))
}
