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
  # Example: add vertical lines at contig starts
  if (!is.null(cxt$mapper$cdf$start)) {
    gg <- gg + ggplot2::geom_vline(xintercept = cxt$mapper$cdf$start, color = "gray")
  }
  gg
}

#' Plot all registered profiles
#' @param cxt The context object containing mapper and other dynamic info
#' @return A plotly object combining all profile plots
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
        param_value <- get_param(id, param_id)()
        cat(sprintf("profile %s, param: %s, value: %s\n", id, param_id, param_value))
        profile[[param_id]] <- param_value
      }
    }

    # Create base ggplot
    gg <- pre_plot(cxt, profile)

    # Get the profile's plot result
    pf_result <- profile$plot_f(profile, cxt, gg)

    # Add overlays
    gg_final <- post_plot(cxt, pf_result, profile)

    # Convert to plotly, ensuring 'text' aesthetic is available for hover
    # The actual 'text' aesthetic should be defined in the profile's plot_f
    p_ly <- plotly::ggplotly(gg_final, tooltip = "text")

    # Store in list
    plotly_list[[id]] <- p_ly
    heights <- c(heights, profile$height)
  }

  if (length(plotly_list) == 0) {
    return(NULL) # Return NULL if no plots were generated
  }

  # Normalize heights to add up to 1 if there are any heights
  if (sum(heights) > 0) {
    heights <- heights / sum(heights)
  } else if (length(heights) > 0) {
    # if all heights are 0, distribute equally
    heights <- rep(1 / length(heights), length(heights))
  }


  # Combine plots
  # Need to handle the case where plotly_list is empty or has only one plot
  if (length(plotly_list) > 1) {
    combined_plot <- plotly::subplot(
      plotly_list,
      nrows = length(plotly_list),
      shareX = TRUE,
      titleY = TRUE,
      margin = 0.01,
      heights = if (length(heights) == length(plotly_list)) heights else NULL
    )
  } else if (length(plotly_list) == 1) {
    combined_plot <- plotly_list[[1]]
  } else {
    return(NULL) # Should not happen if profiles list was not empty initially
  }

  # Configure the combined plot
  combined_plot <- plotly::layout(combined_plot,
    xaxis = list(title = ""),
    margin = list(l = 50, r = 20, t = 30, b = 30) # Adjust margins as needed
  )
  cat(sprintf("plotting done\n"))
  return(combined_plot)
}
