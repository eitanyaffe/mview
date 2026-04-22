default_allele_matrix_params <- list(
  orientation = list(
    group_id = "allele_matrix",
    type = "select",
    default = "triangle",
    choices = c("triangle", "square")
  ),
  style = list(
    group_id = "allele_matrix",
    type = "select",
    default = "cramers_v",
    choices = c("cramers_v", "fgt_excess", "fgt_compat", "n_hap", "p_adj")
  ),
  binsize = list(
    group_id = "allele_matrix",
    type = "select",
    default = "auto",
    choices = c("auto", "full", "10", "50", "100", "200", "500", "1000")
  ),
  show_hover = list(
    group_id = "allele_matrix",
    type = "boolean",
    default = FALSE
  ),
  triangle_aspect = list(
    group_id = "allele_matrix",
    type = "string",
    default = "0.35"
  ),
  height = list(
    group_id = "allele_matrix",
    type = "integer",
    default = 700
  ),
  site_width = list(
    group_id = "allele_matrix",
    type = "integer",
    default = 10
  )
)

allele_matrix_profile <- function(id, name, height = 700, is_fixed = FALSE,
                                  assoc_f = NULL, sites_f = NULL,
                                  params = default_allele_matrix_params,
                                  auto_register = TRUE) {

  cramers_pal <- grDevices::colorRampPalette(c("white", "#2171B5", "#08306B"))(256)
  fgt_pal     <- grDevices::colorRampPalette(c("#2166AC", "#D1E5F0", "white", "#FDDBC7", "#B2182B"))(257)
  nhap_cols   <- c("white", "#9ECAE1", "#2171B5", "#08306B")
  padj_cols   <- c("#CB181D", "#FC4E2A", "#FCBBA1", "white")
  padj_breaks <- c(0, 0.001, 0.01, 0.05, Inf)

  fgt_compat_cols <- c("#ADD8E6", "#FFA500")

  get_fill_color <- function(vals, style) {
    switch(style,
      cramers_v = {
        ix <- as.integer(pmin(pmax(vals * 255 + 1, 1), 256))
        cramers_pal[ix]
      },
      fgt_excess = {
        abs_max <- max(abs(vals), 1L, na.rm = TRUE)
        scaled <- (vals / abs_max + 1) * 0.5
        ix <- as.integer(pmin(pmax(scaled * 256 + 1, 1), 257))
        fgt_pal[ix]
      },
      fgt_compat = ifelse(vals <= 0, fgt_compat_cols[1], fgt_compat_cols[2]),
      n_hap = nhap_cols[pmin(as.integer(pmax(vals, 1)), 4L)],
      p_adj = padj_cols[findInterval(vals, padj_breaks)],
      {
        ix <- as.integer(pmin(pmax(vals * 255 + 1, 1), 256))
        cramers_pal[ix]
      })
  }

  style_col <- function(df, style) {
    if ("mean_cramers_v" %in% names(df)) {
      switch(style,
        cramers_v  = df$mean_cramers_v,
        fgt_excess = df$mean_fgt_excess,
        fgt_compat = df$mean_fgt_excess,
        n_hap      = df$mean_n_hap,
        p_adj      = df$min_p_adj,
        df$mean_cramers_v)
    } else {
      switch(style,
        cramers_v  = df$cramers_v,
        fgt_excess = df$fgt_excess,
        fgt_compat = df$fgt_excess,
        n_hap      = df$n_hap,
        p_adj      = df$p_adj,
        df$cramers_v)
    }
  }

  hover_full <- function(assoc) {
    paste0(
      assoc$site1, " \u2194 ", assoc$site2,
      "<br>Cram\u00e9r's V: ", round(as.numeric(assoc$cramers_v), 3),
      "<br>p_adj: ",           signif(as.numeric(assoc$p_adj), 3),
      if ("n_seq" %in% names(assoc)) paste0("<br>n_seq: ", assoc$n_seq) else "",
      "<br>n_hap: ",           assoc$n_hap,
      "<br>FGT excess: ",      round(as.numeric(assoc$fgt_excess), 3))
  }

  hover_binned <- function(df) {
    paste0(
      df$bin1, " \u2194 ", df$bin2,
      "<br>pairs: ",        df$n_pairs,
      "<br>mean V: ",       round(as.numeric(df$mean_cramers_v), 3),
      "<br>min p_adj: ",    signif(as.numeric(df$min_p_adj), 3),
      if ("mean_n_seq" %in% names(df)) paste0("<br>mean n_seq: ", round(as.numeric(df$mean_n_seq), 1)) else "",
      "<br>mean n_hap: ",   round(as.numeric(df$mean_n_hap), 2),
      "<br>mean FGT: ",     round(as.numeric(df$mean_fgt_excess), 2))
  }

  map_bin_vcoords <- function(assoc, bs) {
    half_bs <- bs * 0.5
    all_contigs <- c(assoc$contig1, assoc$contig2)
    all_bins    <- c(assoc$bin1, assoc$bin2)
    all_mids    <- all_bins + half_bs
    bin_key     <- paste(all_contigs, all_bins, sep = "\t")

    uniq_idx <- !duplicated(bin_key)
    uniq_ctg <- all_contigs[uniq_idx]
    uniq_mid <- all_mids[uniq_idx]

    idf <- data.frame(contig = uniq_ctg,
                      start  = as.integer(uniq_mid),
                      end    = as.integer(uniq_mid + 1L),
                      stringsAsFactors = FALSE)
    filtered <- cxt_contig2view_interval(idf)
    if (is.null(filtered) || nrow(filtered) == 0 || !"vstart" %in% names(filtered)) return(NULL)
    # use midpoint of the mapped interval as the virtual coordinate
    vc_mid <- (filtered$vstart + filtered$vend) * 0.5
    vc_lookup <- setNames(vc_mid, paste(filtered$contig, floor((filtered$start - half_bs) / bs) * bs, sep = "\t"))

    key1 <- paste(assoc$contig1, assoc$bin1, sep = "\t")
    key2 <- paste(assoc$contig2, assoc$bin2, sep = "\t")
    list(vc1 = vc_lookup[key1], vc2 = vc_lookup[key2])
  }

  tri_ylim_for <- function(xlim, profile) {
    aspect <- as.numeric(profile$triangle_aspect)
    if (is.na(aspect) || aspect <= 0) aspect <- 0.15
    c(0, diff(xlim) * aspect)
  }

  plot_binned <- function(profile, gg, assoc, style, orientation, xlim, use_hover = FALSE, binsize = profile$binsize) {
    bs <- as.integer(binsize)
    half_bs <- bs * 0.5
    is_triangle <- orientation == "triangle"

    vc <- map_bin_vcoords(assoc, bs)
    if (is.null(vc))
      return(list(plot = gg, legends = list()))

    vc1 <- vc$vc1
    vc2 <- vc$vc2

    keep <- !is.na(vc1) & !is.na(vc2)
    if (is_triangle) {
      diamond_cx <- (vc1 + vc2) * 0.5
      keep <- keep & diamond_cx >= (xlim[1] - half_bs) & diamond_cx <= (xlim[2] + half_bs)
    }
    if (!any(keep))
      return(list(plot = gg, legends = list()))

    assoc <- assoc[keep, ]
    vc1   <- unname(vc1[keep])
    vc2   <- unname(vc2[keep])
    n     <- nrow(assoc)

    vals  <- style_col(assoc, style)
    fill  <- get_fill_color(vals, style)

    cat(sprintf("[allele_matrix] binned: %d cells (binsize=%d, hover=%s)\n", n, bs, use_hover))

    if (!is_triangle) {
      plot_df <- data.frame(
        xmin = c(vc1 - half_bs, vc2 - half_bs),
        xmax = c(vc1 + half_bs, vc2 + half_bs),
        ymin = c(vc2 - half_bs, vc1 - half_bs),
        ymax = c(vc2 + half_bs, vc1 + half_bs),
        fill = rep.int(fill, 2L),
        stringsAsFactors = FALSE)
      if (use_hover) plot_df$hover <- rep.int(hover_binned(assoc), 2L)

      gg <- gg +
        ggplot2::geom_rect(
          data = plot_df,
          if (use_hover)
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, text = hover)
          else
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
          color = NA) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_cartesian(xlim = xlim, ylim = xlim)
    } else {
      is_sorted <- vc1 <= vc2
      xs <- ifelse(is_sorted, vc1 - half_bs, vc2 - half_bs)
      xe <- ifelse(is_sorted, vc1 + half_bs, vc2 + half_bs)
      ys <- ifelse(is_sorted, vc2 - half_bs, vc1 - half_bs)
      ye <- ifelse(is_sorted, vc2 + half_bs, vc1 + half_bs)

      poly_df <- data.frame(
        x     = c((xs+ys)*0.5, (xe+ys)*0.5, (xe+ye)*0.5, (xs+ye)*0.5),
        y     = pmax(c((ys-xs)*0.5, (ys-xe)*0.5, (ye-xe)*0.5, (ye-xs)*0.5), 0),
        group = rep.int(seq_len(n), 4L),
        fill  = rep.int(fill, 4L),
        stringsAsFactors = FALSE)
      if (use_hover) poly_df$hover <- rep.int(hover_binned(assoc), 4L)

      gg <- gg +
        ggplot2::geom_polygon(
          data = poly_df,
          if (use_hover)
            ggplot2::aes(x = x, y = y, group = group, fill = fill, text = hover)
          else
            ggplot2::aes(x = x, y = y, group = group, fill = fill),
          color = NA) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_cartesian(xlim = xlim, ylim = tri_ylim_for(xlim, profile))
    }

    return(list(plot = gg, legends = list()))
  }

  plot_full <- function(profile, gg, assoc, sites, style, orientation, xlim, use_hover = FALSE) {
    is_triangle <- orientation == "triangle"

    site_df <- data.frame(
      contig  = sites$ref_contig,
      start   = sites$ref_start,
      end     = pmax(sites$ref_end, sites$ref_start + 1L),
      site_id = sites$site_id,
      stringsAsFactors = FALSE)

    if (is_triangle) {
      filtered_sites <- cxt_contig2view_interval(site_df)
    } else {
      filtered_sites <- cxt_filter_intervals(site_df, merge_adjacent = FALSE)
    }
    if (is.null(filtered_sites) || nrow(filtered_sites) == 0 || !"vstart" %in% names(filtered_sites))
      return(list(plot = gg, legends = list()))

    site_ids <- filtered_sites$site_id
    vstart   <- filtered_sites$vstart
    vend     <- filtered_sites$vend
    vmid     <- (vstart + vend) * 0.5

    # expand each site symmetrically to site_width, capped at midpoint to neighbor
    site_width <- as.numeric(if (!is.null(profile$site_width)) profile$site_width else 10)
    half_w     <- site_width / 2
    ord        <- order(vmid)
    ms         <- vmid[ord]
    n_s        <- length(ms)
    if (n_s > 1L) {
      between   <- (ms[-n_s] + ms[-1L]) / 2
      lim_left  <- c(-Inf, between)
      lim_right <- c(between, Inf)
    } else {
      lim_left  <- -Inf
      lim_right <- Inf
    }
    vs_exp        <- numeric(n_s)
    ve_exp        <- numeric(n_s)
    vs_exp[ord]   <- pmax(ms - half_w, lim_left)
    ve_exp[ord]   <- pmin(ms + half_w, lim_right)
    vstart <- vs_exp
    vend   <- ve_exp

    # pre-filter sites and assoc to the visible range
    view_margin <- diff(xlim) * 0.5
    in_view     <- vmid >= (xlim[1] - view_margin) & vmid <= (xlim[2] + view_margin)
    if (!any(in_view))
      return(list(plot = gg, legends = list()))
    visible_ids <- site_ids[in_view]
    assoc <- assoc[assoc$site1 %in% visible_ids & assoc$site2 %in% visible_ids, ]
    if (nrow(assoc) == 0)
      return(list(plot = gg, legends = list()))

    ix1  <- match(assoc$site1, site_ids)
    ix2  <- match(assoc$site2, site_ids)
    keep <- !is.na(ix1) & !is.na(ix2)
    if (!any(keep))
      return(list(plot = gg, legends = list()))

    ix1   <- ix1[keep]
    ix2   <- ix2[keep]
    assoc <- assoc[keep, ]
    n     <- nrow(assoc)

    vals  <- style_col(assoc, style)
    fill  <- get_fill_color(vals, style)

    cat(sprintf("[allele_matrix] full: %d pairs (hover=%s)\n", n, use_hover))

    if (!is_triangle) {
      plot_df <- data.frame(
        xmin = c(vstart[ix1], vstart[ix2]),
        xmax = c(vend[ix1],   vend[ix2]),
        ymin = c(vstart[ix2], vstart[ix1]),
        ymax = c(vend[ix2],   vend[ix1]),
        fill = rep.int(fill, 2L),
        stringsAsFactors = FALSE)
      if (use_hover) plot_df$hover <- rep.int(hover_full(assoc), 2L)

      gg <- gg +
        ggplot2::geom_rect(
          data = plot_df,
          if (use_hover)
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, text = hover)
          else
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
          color = NA) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_cartesian(xlim = xlim, ylim = xlim)
    } else {
      is_sorted <- vmid[ix1] <= vmid[ix2]
      xs <- ifelse(is_sorted, vstart[ix1], vstart[ix2])
      xe <- ifelse(is_sorted, vend[ix1],   vend[ix2])
      ys <- ifelse(is_sorted, vstart[ix2], vstart[ix1])
      ye <- ifelse(is_sorted, vend[ix2],   vend[ix1])

      poly_df <- data.frame(
        x     = c((xs+ys)*0.5, (xe+ys)*0.5, (xe+ye)*0.5, (xs+ye)*0.5),
        y     = pmax(c((ys-xs)*0.5, (ys-xe)*0.5, (ye-xe)*0.5, (ye-xs)*0.5), 0),
        group = rep.int(seq_len(n), 4L),
        fill  = rep.int(fill, 4L),
        stringsAsFactors = FALSE)
      if (use_hover) poly_df$hover <- rep.int(hover_full(assoc), 4L)

      ylim_tri <- tri_ylim_for(xlim, profile)
      ymax     <- ylim_tri[2]
      site_edges <- unique(c(vstart, vend))
      grid_df <- data.frame(
        x    = c(site_edges, site_edges),
        xend = c(site_edges - ymax, site_edges + ymax),
        y    = 0,
        yend = ymax)

      gg <- gg +
        # ggplot2::geom_segment(
        #   data = grid_df,
        #   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        #   color = "gray85", linewidth = 0.3) +
        ggplot2::geom_polygon(
          data = poly_df,
          if (use_hover)
            ggplot2::aes(x = x, y = y, group = group, fill = fill, text = hover)
          else
            ggplot2::aes(x = x, y = y, group = group, fill = fill),
          color = NA) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim_tri)
    }

    return(list(plot = gg, legends = list()))
  }

  build_matrix_legend <- function(style) {
    if (style == "cramers_v") {
      cols <- grDevices::colorRampPalette(c("white", "#2171B5", "#08306B"))(8)
      vals <- seq(0, 1, length.out = 8)
      legend_data <- data.frame(y = seq_along(vals), x = 1, color = cols,
                                label = sprintf("%.2f", vals), stringsAsFactors = FALSE)
      gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35,
                                        ymin = y - 0.45, ymax = y + 0.45, fill = color),
                           color = "black", size = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.2) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "Cram\u00e9r's V") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                       plot.margin = ggplot2::margin(8, 8, 8, 8)) +
        ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(vals) + 0.5))
      list(list(gg = gg, height = 280, width = 250, title = "Cram\u00e9r's V"))
    } else if (style == "fgt_excess") {
      cols <- grDevices::colorRampPalette(c("#2166AC", "#D1E5F0", "white", "#FDDBC7", "#B2182B"))(7)
      labels <- c("-1.0", "-0.67", "-0.33", "0", "0.33", "0.67", "1.0")
      legend_data <- data.frame(y = seq_along(labels), x = 1, color = cols,
                                label = labels, stringsAsFactors = FALSE)
      gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35,
                                        ymin = y - 0.45, ymax = y + 0.45, fill = color),
                           color = "black", size = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.2) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "FGT excess") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                       plot.margin = ggplot2::margin(8, 8, 8, 8)) +
        ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(labels) + 0.5))
      list(list(gg = gg, height = 260, width = 250, title = "FGT excess"))
    } else if (style == "n_hap") {
      cols   <- nhap_cols
      labels <- c("1", "2", "3", "4")
      legend_data <- data.frame(y = seq_along(labels), x = 1, color = cols,
                                label = labels, stringsAsFactors = FALSE)
      gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35,
                                        ymin = y - 0.45, ymax = y + 0.45, fill = color),
                           color = "black", size = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.2) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "haplotypes") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                       plot.margin = ggplot2::margin(8, 8, 8, 8)) +
        ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(labels) + 0.5))
      list(list(gg = gg, height = 180, width = 250, title = "haplotypes"))
    } else if (style == "fgt_compat") {
      cols   <- fgt_compat_cols
      labels <- c("compatible", "incompatible")
      legend_data <- data.frame(y = seq_along(labels), x = 1, color = cols,
                                label = labels, stringsAsFactors = FALSE)
      gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35,
                                        ymin = y - 0.45, ymax = y + 0.45, fill = color),
                           color = "black", size = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.2) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "tree compatibility") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                       plot.margin = ggplot2::margin(8, 8, 8, 8)) +
        ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(labels) + 0.5))
      list(list(gg = gg, height = 120, width = 250, title = "tree compatibility"))
    } else if (style == "p_adj") {
      cols   <- padj_cols
      labels <- c("< 0.001", "0.001-0.01", "0.01-0.05", "> 0.05")
      legend_data <- data.frame(y = seq_along(labels), x = 1, color = cols,
                                label = labels, stringsAsFactors = FALSE)
      gg <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.35, xmax = x + 0.35,
                                        ymin = y - 0.45, ymax = y + 0.45, fill = color),
                           color = "black", size = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = label), x = 1.7, hjust = 0, size = 3.2) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = "adjusted p-value") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                       plot.margin = ggplot2::margin(8, 8, 8, 8)) +
        ggplot2::coord_cartesian(xlim = c(0.5, 3), ylim = c(0.5, length(labels) + 0.5))
      list(list(gg = gg, height = 180, width = 250, title = "adjusted p-value"))
    } else {
      list()
    }
  }

  plot_f <- function(profile, gg) {
    assembly <- cxt_get_assembly()
    binsize  <- if (!is.null(profile$binsize)) profile$binsize else "500"
    style    <- if (!is.null(profile$style)) profile$style else "cramers_v"
    orientation <- if (!is.null(profile$orientation)) profile$orientation else "triangle"
    use_hover <- isTRUE(profile$show_hover)
    xlim     <- cxt_get_xlim()

    matrix_legends <- build_matrix_legend(style)

    if (binsize == "auto") {
      view_width <- diff(xlim)
      binsize <- if      (view_width <   8000) "10"
                 else if (view_width <  20000) "50"
                 else if (view_width <  60000) "100"
                 else if (view_width < 150000) "200"
                 else if (view_width < 400000) "500"
                 else                          "1000"
    }

    cat(sprintf("[allele_matrix] binsize=%s, style=%s, orientation=%s, hover=%s\n",
                binsize, style, orientation, use_hover))

    if (binsize == "full") {
      assoc <- if (is.function(profile$assoc_f)) profile$assoc_f(assembly, "full") else NULL
      sites <- if (is.function(profile$sites_f)) profile$sites_f(assembly) else NULL
      if (is.null(assoc) || nrow(assoc) == 0 || is.null(sites) || nrow(sites) == 0)
        return(list(plot = gg, legends = matrix_legends))

      cat(sprintf("[allele_matrix] full data: %d pairs, %d sites\n", nrow(assoc), nrow(sites)))
      result <- plot_full(profile, gg, assoc, sites, style, orientation, xlim, use_hover)
      result$legends <- matrix_legends
      return(result)
    } else {
      assoc <- if (is.function(profile$assoc_f)) profile$assoc_f(assembly, binsize) else NULL
      if (is.null(assoc) || nrow(assoc) == 0)
        return(list(plot = gg, legends = matrix_legends))
      cat(sprintf("[allele_matrix] binned data: %d pairs\n", nrow(assoc)))
      result <- plot_binned(profile, gg, assoc, style, orientation, xlim, use_hover, binsize = binsize)
      result$legends <- matrix_legends
      return(result)
    }
  }

  profile_create(
    id = id, name = name, type = "allele_matrix", height = height,
    is_fixed = is_fixed,
    attr = list(hide_y_ticks = TRUE),
    params = params, plot_f = plot_f,
    assoc_f = assoc_f,
    sites_f = sites_f,
    auto_register = auto_register)
}
