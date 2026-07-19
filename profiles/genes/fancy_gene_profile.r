# fancy gene profile - three zoom levels, three color modes

fancy_gene_profile <- function(id, name, height = 40, is_fixed = TRUE,
                                gene_f = NULL,
                                threshold_simple = 500000,
                                threshold_full   = 50000,
                                auto_register = TRUE) {

  params <- list(
    color_mode = list(
      group_id = id,
      type     = "select",
      choices  = c("basic", "strand", "uniref"),
      default  = "basic"
    )
  )

  plot_f <- function(profile, gg) {
    genes <- gene_f(cxt_get_assembly())
    if (is.null(genes) || nrow(genes) == 0)
      return(list(plot = gg, legends = list()))

    df <- genes
    df$contig <- as.character(df$contig)
    df$start  <- as.numeric(df$start)
    df$end    <- as.numeric(df$end)

    visible_contigs <- cxt_get_contigs()
    df <- df[df$contig %in% visible_contigs, ]
    if (nrow(df) == 0)
      return(list(plot = gg, legends = list()))

    valid <- !is.na(df$contig) & !is.na(df$start) & !is.na(df$end)
    df <- df[valid, ]
    if (nrow(df) == 0)
      return(list(plot = gg, legends = list()))

    tryCatch({
      df <- cxt_filter_intervals(df, merge_adjacent = TRUE)
    }, error = function(e) {
      df <<- NULL
    })
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
      return(list(plot = gg, legends = list()))

    xlim <- cxt_get_xlim()

    df$gstart <- df$vstart
    df$gend   <- df$vend

    color_mode <- if (!is.null(profile$color_mode)) profile$color_mode else "basic"
    xspan      <- diff(xlim)

    # determine zoom mode
    zoom_mode <- if (xspan > threshold_simple) "simple" else if (xspan > threshold_full) "medium" else "full"

    col_basic      <- "#ADD8E6"
    col_strand_neg <- "#F08080"

    if (color_mode == "basic") {
      df$color <- col_basic
    } else if (color_mode == "strand") {
      df$color <- ifelse(df$strand == "+", col_basic, col_strand_neg)
    } else {
      # uniref: color by uniref-derived taxonomy
      if ("tax_color" %in% names(df) && !all(is.na(df$tax_color))) {
        df$color <- ifelse(is.na(df$tax_color) | df$tax_color == "", col_basic, df$tax_color)
      } else {
        df$color <- col_basic
      }
    }

    # strand geometry: TSS is at gstart for + genes (adjusted for segment flip)
    seg_strand    <- if ("segment_strand" %in% names(df)) df$segment_strand else "+"
    tss_at_gstart <- (df$strand == "+") != (seg_strand == "-")
    df$tss_pos      <- ifelse(tss_at_gstart, df$gstart, df$gend)
    df$gene_end_pos <- ifelse(tss_at_gstart, df$gend, df$gstart)
    df$dir          <- ifelse(tss_at_gstart, 1L, -1L)

    df$hover_text <- paste0(
      "Gene: ", df$gene, "\n",
      if ("prot_desc" %in% names(df)) paste0("Description: ", df$prot_desc, "\n") else "",
      if ("tax" %in% names(df)) paste0("Taxonomy: ", df$tax) else ""
    )

    legends <- list()

    if (zoom_mode == "simple") {
      # TSS vertical segments only
      if (nrow(df) > 5000) {
        set.seed(42)
        df <- df[sample(nrow(df), 5000), ]
      }
      gg <- gg + ggplot2::geom_segment(
        data = df,
        ggplot2::aes(x = tss_pos, xend = tss_pos, y = 0.45, yend = 0.55, text = hover_text),
        color = df$color,
        size  = 1
      )

    } else if (zoom_mode == "medium") {
      # colored rects, no arrows or L-shapes
      df$rect_ymin <- 0.35
      df$rect_ymax <- 0.65

      gg <- gg +
        ggplot2::geom_rect(
          data = df,
          ggplot2::aes(xmin = gstart, xmax = gend,
                       ymin = rect_ymin, ymax = rect_ymax,
                       fill = I(color), text = hover_text),
          alpha = 0.6
        ) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::coord_cartesian(xlim = cxt_get_xlim(), ylim = c(0.3, 0.7))

    } else {
      # full: rects + vertical lines + arrows + L-shapes
      arr_len  <- xspan * 0.010
      tri_len  <- arr_len * 0.65
      tri_h    <- 0.11
      y_lo     <- 0.35
      y_hi     <- 0.65
      arr_y    <- y_hi + 0.03
      tick_y   <- y_lo - 0.03

      df$rect_ymin  <- y_lo
      df$rect_ymax  <- y_hi
      df$arr_tip_x  <- df$tss_pos + df$dir * arr_len
      df$arr_base_x <- df$arr_tip_x - df$dir * tri_len
      df$l_xend     <- df$gene_end_pos - df$dir * arr_len * 0.6

      gg <- gg +
        ggplot2::geom_rect(
          data = df,
          ggplot2::aes(xmin = gstart, xmax = gend,
                       ymin = rect_ymin, ymax = rect_ymax,
                       fill = I(color), text = hover_text),
          alpha = 0.6
        ) +
        ggplot2::theme(legend.position = "none")

      # vertical line at TSS: crosses through rect and continues up to arrow row
      gg <- gg + ggplot2::geom_segment(
        data = df,
        ggplot2::aes(x = tss_pos, xend = tss_pos, y = y_lo, yend = arr_y, text = hover_text),
        color = "black", size = 0.5, lineend = "round"
      )

      # horizontal arrow shaft (stops at triangle base)
      gg <- gg + ggplot2::geom_segment(
        data = df,
        ggplot2::aes(x = tss_pos, xend = arr_base_x, y = arr_y, yend = arr_y, text = hover_text),
        color = "black", size = 0.5, lineend = "round", linejoin = "round"
      )

      # arrowhead as filled triangle polygon (survives ggplotly conversion)
      n_genes <- nrow(df)
      tri_df  <- data.frame(
        x     = c(rbind(df$arr_tip_x, df$arr_base_x, df$arr_base_x)),
        y     = c(rbind(rep(arr_y, n_genes), rep(arr_y + tri_h/2, n_genes), rep(arr_y - tri_h/2, n_genes))),
        group = rep(seq_len(n_genes), each = 3),
        stringsAsFactors = FALSE
      )
      gg <- gg + ggplot2::geom_polygon(
        data = tri_df,
        ggplot2::aes(x = x, y = y, group = group),
        fill = "black", color = NA
      )

      # vertical line at gene end: crosses through rect and continues down to L tick
      gg <- gg + ggplot2::geom_segment(
        data = df,
        ggplot2::aes(x = gene_end_pos, xend = gene_end_pos, y = y_hi, yend = tick_y, text = hover_text),
        color = "gray40", size = 0.4, lineend = "round"
      )

      # L horizontal pointing back toward gene interior
      gg <- gg + ggplot2::geom_segment(
        data = df,
        ggplot2::aes(x = gene_end_pos, xend = l_xend, y = tick_y, yend = tick_y, text = hover_text),
        color = "gray40", size = 0.4, lineend = "round"
      )

      gg <- gg + ggplot2::coord_cartesian(xlim = cxt_get_xlim(), ylim = c(0.3, 0.7))
    }

    if (color_mode == "strand") {
      legend_data <- data.frame(
        y     = c(1, 2),
        x     = 1,
        color = c(col_basic, col_strand_neg),
        label = c("+ strand", "- strand"),
        stringsAsFactors = FALSE
      )
      legend_gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35,
                       ymin = y - 0.45, ymax = y + 0.45, fill = color),
          color = "black", size = 0.3
        ) +
        ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.4) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "strand") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                       plot.margin = ggplot2::margin(8, 8, 8, 8)) +
        ggplot2::coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0.5, 2.5))
      legends <- list(list(gg = legend_gg, height = 100, width = 200, title = "strand"))

    } else if (color_mode == "uniref" && "tax" %in% names(df)) {
      valid_tax <- df[!is.na(df$tax) & df$tax != "none" & df$color != col_basic, ]
      if (nrow(valid_tax) > 0) {
        tax_counts  <- sort(table(valid_tax$tax), decreasing = TRUE)
        top_n       <- min(10, length(tax_counts))
        top_tax     <- head(tax_counts, top_n)
        legend_data <- data.frame(
          category = names(top_tax),
          count    = as.numeric(top_tax),
          stringsAsFactors = FALSE
        )
        legend_data$color <- sapply(legend_data$category, function(t) {
          idx <- which(valid_tax$tax == t)[1]
          if (!is.na(idx)) valid_tax$color[idx] else col_basic
        })
        legend_data$label <- paste0(legend_data$category, " (", legend_data$count, ")")
        legend_data <- legend_data[nrow(legend_data):1, ]

        n <- nrow(legend_data)
        legend_gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = 1, y = seq_len(n))) +
          ggplot2::geom_rect(
            ggplot2::aes(xmin = 0.65, xmax = 1.35,
                         ymin = seq_len(n) - 0.45, ymax = seq_len(n) + 0.45,
                         fill = color),
            color = "black", size = 0.3
          ) +
          ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.4) +
          ggplot2::scale_fill_identity() +
          ggplot2::labs(title = "top taxa (uniref)") +
          ggplot2::theme_void() +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                         plot.margin = ggplot2::margin(8, 8, 8, 8)) +
          ggplot2::coord_cartesian(xlim = c(0.5, 4), ylim = c(0.5, n + 0.5))
        legends <- list(list(gg = legend_gg, height = 60 + n * 25, width = 750,
                             title = "taxa (uniref)"))
      }
    }

    list(plot = gg, legends = legends)
  }

  profile_create(
    id               = id,
    name             = name,
    type             = "gene",
    height           = height,
    is_fixed         = is_fixed,
    attr             = list(hide_y_ticks = TRUE),
    params           = params,
    plot_f           = plot_f,
    gene_f           = gene_f,
    threshold_simple = threshold_simple,
    threshold_full   = threshold_full,
    auto_register    = auto_register
  )
}
