gene_profile <- function(id, name, height = 0.1,
                         gene_f = NULL,
                         color_style = "by_taxonomy",
                         select_groups = NULL,
                         select_colors = NULL,
                         threshold = 200000,
                         params = list(),
                         auto_register = TRUE) {
  # Determine display mode based on view range
  get_display_mode <- function(xlim, threshold) {
    if (is.null(xlim) || length(xlim) != 2) {
      return("simple")
    }

    range_bp <- (xlim[2] + 1) - xlim[1]

    if (range_bp <= threshold) {
      return("full")
    } else {
      return("simple")
    }
  }

  plot_f <- function(profile, cxt, gg) {
    genes <- gene_f(cxt)
    if (is.null(genes)) {
      return(NULL)
    }
    mode <- get_display_mode(cxt$mapper$xlim, profile$threshold)
    cat(sprintf("mode: %s\n", mode))
    plot_gene_profile(profile, cxt, genes, gg, mode)
  }

  # Create profile
  profile_create(
    id = id, name = name, type = "gene", height = height,
    params = params, plot_f = plot_f,
    threshold = threshold,
    gene_f = gene_f,
    color_style = color_style,
    select_groups = select_groups,
    select_colors = select_colors,
    auto_register = auto_register
  )
}
