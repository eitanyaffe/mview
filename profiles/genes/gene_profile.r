# default parameters for gene profile
default_gene_params <- list(
  color_field = list(
    group_id = "gene",
    type = "select",
    choices = c("tax", "mge"),
    default = "tax"
  ),
  height = list(
    group_id = "gene",
    type = "integer",
    default = 50
  ),
  show_hover = list(
    group_id = "gene",
    type = "boolean",
    default = TRUE
  ),
  max_items_in_legend = list(
    group_id = "gene",
    type = "integer",
    default = 10
  )
)

gene_profile <- function(id, name, height = 100,
                         gene_f = NULL,
                         color_field = "tax_color",
                         threshold = 200000,
                         label_field = NULL,
                         show_hover = TRUE,
                         params = default_gene_params,
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
      return(list(plot = gg, legends = list()))
    }
    mode <- get_display_mode(cxt$mapper$xlim, profile$threshold)
    cat(sprintf("mode: %s\n", mode))
    return(plot_gene_profile(profile, cxt, genes, gg, mode))
  }

  # Create profile
  profile_create(
    id = id, name = name, type = "gene", height = height,
    attr = list(hide_y_ticks = TRUE),
    params = params, plot_f = plot_f,
    threshold = threshold,
    gene_f = gene_f,
    color_field = color_field,
    label_field = label_field,
    show_hover = show_hover,
    auto_register = auto_register
  )
}
