# core/profile_viewer.r

# This file handles the plotting logic for all registered profiles.
# For each profile, it calls pre_plot, profile$plot_f, and post_plot, and collects the results.

#' Pre-plot function: sets up the base ggplot object for a profile
#' @param cxt Context object
#' @param profile The profile object
#' @return A ggplot object with base layers (grid, y-label, etc.)
pre_plot <- function(cxt, profile) {
  # Minimal base plot with xlim from cxt$mapper, y-label from profile$attr$title
  gg <- ggplot2::ggplot() +
    ggplot2::xlim(cxt$mapper$xlim) +
    ggplot2::labs(y = profile$attr$title) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
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
#' @return A list of lists: each with plot (ggplot) and info_df (for hover)
plot_profiles <- function(cxt) {
  profiles <- profiles_get_all()
  results <- list()
  for (id in names(profiles)) {
    cat(sprintf("Plotting profile: %s\n", id))
    profile <- profiles[[id]]
    gg <- pre_plot(cxt, profile)
    pf_result <- profile$plot_f(cxt, gg)
    gg_final <- post_plot(cxt, pf_result$plot, profile)
    results[[id]] <- list(plot = gg_final, info_df = pf_result$info_df, profile = profile)
  }
  results
}
