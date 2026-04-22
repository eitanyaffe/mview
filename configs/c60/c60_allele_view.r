source("profiles/interval_profile.r")
source("profiles/variants_profile.r")
source("profiles/rearrangements_profile.r")
source("profiles/rrna_profile.r")
source("profiles/allele_matrix_profile.r")
source("profiles/line_profile.r")
source("profiles/density_profile.r")

########################################################
# malign allele matrix profile
########################################################

allele_matrix_profile(
  id      = "allele_matrix",
  name    = "Allele Matrix",
  assoc_f = get_malign_assoc_f,
  sites_f = get_malign_transform_sites_f
)

########################################################
# malign hotspots interval profile (with window param)
########################################################

local({
  params <- list(
    window = list(
      group_id = "hotspots",
      type     = "integer",
      default  = 100
    )
  )

  plot_f <- function(profile, gg) {
    assembly <- cxt_get_assembly()
    hs <- get_malign_hotspots_f(assembly)
    if (is.null(hs) || nrow(hs) == 0)
      return(list(plot = gg, legends = list()))

    w  <- as.integer(if (!is.null(profile$window)) profile$window else min(hs$window))
    hs <- hs[hs$window == w, ]
    if (nrow(hs) == 0)
      return(list(plot = gg, legends = list()))

    filtered <- cxt_filter_intervals(hs, merge_adjacent = FALSE)
    if (is.null(filtered) || nrow(filtered) == 0)
      return(list(plot = gg, legends = list()))

    xlim       <- cxt_get_xlim()
    hover_text <- paste0(filtered$contig, ": ", filtered$start, "-", filtered$end,
                         "\n", filtered$desc)

    gg <- gg +
      ggplot2::geom_rect(
        data = filtered,
        ggplot2::aes(xmin = gstart, xmax = gend, ymin = -0.15, ymax = 0.15, text = hover_text),
        fill = "#CC3333", color = "black", size = 0.5
      ) +
      ggplot2::geom_text(
        data = filtered,
        ggplot2::aes(x = gstart + diff(xlim) * 0.01, y = 0.25, label = id, text = hover_text),
        color = "black", size = 2, hjust = 1, vjust = 0.5
      ) +
      ggplot2::ylim(-0.5, 0.6)

    list(plot = gg, legends = list())
  }

  profile_create(
    id            = "hotspots",
    name          = "Hotspots",
    type          = "intervals",
    height        = 60,
    is_fixed      = TRUE,
    attr          = list(hide_y_ticks = TRUE),
    params        = params,
    plot_f        = plot_f,
    auto_register = TRUE
  )
})

########################################################
# site density profile (sites per kb, binned)
########################################################

density_profile(
  id          = "site_density",
  name        = "Site Count",
  points_f    = get_allele_density_f,
  height      = 80,
  color       = "#2171B5",
  param_group = "malign_density"
)

########################################################
# malign allele sites profile (with color_mode param)
########################################################

local({
  params <- list(
    height = list(
      group_id = "malign_sites",
      type     = "integer",
      default  = 60
    ),
    color_mode = list(
      group_id = "malign_sites",
      type     = "select",
      default  = "gene",
      choices  = c("gene", "n_alleles")
    ),
    show_hover = list(
      group_id = "malign_sites",
      type     = "boolean",
      default  = FALSE
    )
  )

  plot_f <- function(profile, gg) {
    assembly  <- cxt_get_assembly()
    intervals <- profile$intervals_f(assembly)

    if (is.null(intervals) || nrow(intervals) == 0)
      return(list(plot = gg, legends = list()))

    if ("assembly" %in% names(intervals))
      intervals <- intervals[intervals$assembly == assembly, ]
    if (is.null(intervals) || nrow(intervals) == 0)
      return(list(plot = gg, legends = list()))

    filtered <- cxt_filter_intervals(intervals, merge_adjacent = FALSE)
    if (is.null(filtered) || nrow(filtered) == 0)
      return(list(plot = gg, legends = list()))

    color_mode  <- if (!is.null(profile$color_mode)) profile$color_mode else "gene"
    col_name    <- if (color_mode == "n_alleles") "fill_nalleles_color" else "fill_gene_color"
    fill_colors <- filtered[[col_name]]
    if (is.null(fill_colors)) fill_colors <- rep("#AAAAAA", nrow(filtered))
    fill_colors[is.na(fill_colors) | fill_colors == ""] <- "#AAAAAA"

    use_hover <- isTRUE(profile$show_hover)
    if (use_hover && "site_id" %in% names(filtered)) {
      ann_sites   <- get_malign_annotate_sites_f()
      ann_alleles <- get_malign_annotate_alleles_f()
      has_ann     <- !is.null(ann_sites) && nrow(ann_sites) > 0 &&
                     !is.null(ann_alleles) && nrow(ann_alleles) > 0
      hover_text <- if (has_ann) {
        mapply(build_malign_site_hover,
               filtered$site_id, filtered$n_alleles,
               MoreArgs = list(ann_sites = ann_sites, ann_alleles = ann_alleles),
               SIMPLIFY = TRUE)
      } else {
        paste0("site: ", filtered$site_id, "\nalleles: ", filtered$n_alleles)
      }
    } else {
      hover_text <- ""
    }

    filtered$fill_color <- fill_colors

    gg <- gg +
      ggplot2::geom_rect(
        data = filtered,
        ggplot2::aes(xmin = gstart, xmax = gend, ymin = -0.15, ymax = 0.15,
                     fill = fill_color, color = fill_color, text = hover_text),
        size = 0.5) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity() +
      ggplot2::ylim(-0.5, 0.6)

    legends <- list()
    if (color_mode == "gene") {
      labels <- c("intergenic", "syn / identical", "non-syn", "aa indel", "frameshift")
      colors <- c("#AAAAAA", "#ADD8E6", "#FFA500", "#CC88FF", "#E74C3C")
    } else {
      labels <- c("2 alleles", "3 alleles", ">3 alleles")
      colors <- c("#ADD8E6", "#FFA500", "#E74C3C")
    }
    legend_data <- data.frame(
      y = seq_along(labels), x = 1,
      color = colors, label = labels,
      stringsAsFactors = FALSE
    )
    legend_gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35, ymin = y - 0.45, ymax = y + 0.45, fill = color),
        color = "black", size = 0.3
      ) +
      ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.4) +
      ggplot2::scale_fill_identity() +
      ggplot2::labs(title = paste0("site colors (", color_mode, ")")) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                     plot.margin = ggplot2::margin(8, 8, 8, 8)) +
      ggplot2::coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0.5, length(labels) + 0.5))
    legends <- list(list(gg = legend_gg, height = 40 + length(labels) * 30, width = 280,
                         title = "site colors"))

    list(plot = gg, legends = legends)
  }

  profile_create(
    id            = "malign_sites",
    name          = "Allele Sites",
    type          = "intervals",
    height        = 60,
    is_fixed      = TRUE,
    attr          = list(hide_y_ticks = TRUE),
    params        = params,
    plot_f        = plot_f,
    intervals_f   = get_malign_sites_f,
    auto_register = TRUE
  )
})

########################################################
# fancy gene profile
########################################################

source("profiles/genes/fancy_gene_profile.r")

fancy_gene_profile(
  id     = "fancy_genes",
  name   = "Genes (fancy)",
  gene_f = get_genes_f
)

########################################################
# rRNA gene profile
########################################################

rrna_profile(
  id = "rrna",
  name = "rRNA",
  get_gff_f = function() {
    get_data("BARRNAP_TABLE", tag = cxt_get_assembly(), null.on.missing = TRUE, read_f = read_gff_f)
  }
)

########################################################
# regions interval profile
########################################################

interval_profile(
  id = "regions",
  name = "Regions",
  intervals_f = "segments.current_regions",
  merge_adjacent = TRUE,
  height = 60
)

########################################################
# bin segments interval profile
########################################################

interval_profile(
  id = "bin_segments",
  name = "Segments",
  intervals_f = get_bin_segments_f,
  color_f = function(ids) { get_current_color_map()(ids) },
  merge_adjacent = FALSE,
  height = 60
)

########################################################
# variants profile
########################################################

variants_profile(
  id = "variants",
  name = "Variants"
)

########################################################
# rearrangements profile
########################################################

rearrangements_profile(
  id = "rearrangements",
  name = "Rearrangements"
)

########################################################
# axis profile
########################################################

axis_profile(height = 80)
