########################################################
# RNA gene expression profile
#
# Reads RNA_ASS_GENE_TABLE (per subject) via get_gene_table_f.
# Renders one horizontal segment per gene (start–end), y-axis = metric.
# Multiple timepoints overlaid with distinct colors.
#
# Metrics:
#   rna_count, rna_rpkm, dna_count, dna_rpkm, rscore,
#   rna_over_dna  (normalized ratio with pseudocount)
#
# For rna_over_dna normalization uses total_mapped from
# RNA_ASS_LIB_STATS via get_lib_stats_f.
########################################################

# re-use the same color palette as the coverage profile
if (!exists("RNA_TIMEPOINT_COLORS"))
    RNA_TIMEPOINT_COLORS <- c(
        "D-1"  = "#2171B5",
        "D+4"  = "#D94801",
        "D+7"  = "#238B45",
        "D+18" = "#6A3D9A"
    )

default_rna_gene_params <- function(id, shared_id = id) {
    list(
        # shared across RNA/DNA cov and RNA gene profiles
        show_Dm1    = list(group_id = shared_id, type = "boolean", default = TRUE),
        show_Dp4    = list(group_id = shared_id, type = "boolean", default = TRUE),
        show_Dp7    = list(group_id = shared_id, type = "boolean", default = TRUE),
        show_Dp18   = list(group_id = shared_id, type = "boolean", default = TRUE),
        pseudocount = list(group_id = shared_id, type = "double",  default = 0.1),
        ratio_ymax  = list(group_id = shared_id, type = "double",  default = 3),
        # gene-profile-specific
        metric    = list(group_id = id, type = "select",
                         choices  = c("rna_count", "rna_rpkm", "dna_count", "dna_rpkm",
                                      "rscore", "rna_over_dna"),
                         default  = "rna_rpkm"),
        host_only = list(group_id = id, type = "boolean", default = FALSE),
        log_y     = list(group_id = id, type = "boolean", default = FALSE),
        height    = list(group_id = id, type = "integer", default = 100)
    )
}

rna_gene_profile <- function(id, name,
                              get_gene_table_f,
                              get_lib_stats_f,
                              sample_types  = c("D-1", "D+4", "D+7", "D+18"),
                              height        = 100,
                              shared_id     = id,
                              auto_register = TRUE) {

    params <- default_rna_gene_params(id, shared_id)

    plot_f <- function(profile, gg) {
        assembly <- cxt_get_assembly()

        genes <- get_gene_table_f(assembly)
        if (is.null(genes) || nrow(genes) == 0)
            return(list(plot = gg, legends = list()))

        metric    <- profile$metric
        host_only <- isTRUE(profile$host_only)
        pseudo    <- profile$pseudocount
        log_y     <- isTRUE(profile$log_y)

        # rscore is only defined for host-assigned genes
        if (metric == "rscore") host_only <- TRUE

        if (host_only)
            genes <- genes[!is.na(genes$host_bin), ]
        if (nrow(genes) == 0)
            return(list(plot = gg, legends = list()))

        # map gene intervals to virtual coordinates; carry row_id for join-back
        genes$.row_id <- seq_len(nrow(genes))
        intervals_df <- data.frame(
            contig  = genes$contig,
            start   = genes$start,
            end     = genes$end,
            .row_id = genes$.row_id,
            stringsAsFactors = FALSE
        )
        mapped <- cxt_filter_intervals(intervals_df)
        if (is.null(mapped) || nrow(mapped) == 0)
            return(list(plot = gg, legends = list()))

        genes <- merge(genes, mapped[, c(".row_id", "vstart", "vend")], by = ".row_id")
        if (nrow(genes) == 0)
            return(list(plot = gg, legends = list()))

        show_flags <- c(
            "D-1"  = isTRUE(profile$show_Dm1),
            "D+4"  = isTRUE(profile$show_Dp4),
            "D+7"  = isTRUE(profile$show_Dp7),
            "D+18" = isTRUE(profile$show_Dp18)
        )
        active_types <- intersect(sample_types, names(show_flags)[show_flags])
        if (length(active_types) == 0)
            return(list(plot = gg, legends = list()))

        stats <- if (metric == "rna_over_dna") get_lib_stats_f(assembly) else NULL

        ratio_data_max <- -Inf

        for (st in active_types) {
            color  <- RNA_TIMEPOINT_COLORS[st]
            if (is.na(color)) color <- "#888888"

            y <- .gene_metric(genes, st, metric, pseudo, stats)
            if (is.null(y)) next

            # track ratio max in raw (linear) space for the ymax cap
            if (metric == "rna_over_dna") {
                m <- suppressWarnings(max(y, na.rm = TRUE))
                if (is.finite(m)) ratio_data_max <- max(ratio_data_max, m)
            }

            if (log_y)
                y <- log2(pmax(y, 1e-9))

            valid <- !is.na(y)
            if (!any(valid)) next

            plot_df <- data.frame(
                xstart = genes$vstart[valid],
                xend   = genes$vend[valid],
                y      = y[valid],
                label  = paste0(genes$gene[valid], "\n", st, " ", metric, ": ",
                                round(y[valid], 4)),
                stringsAsFactors = FALSE
            )

            gg <- gg + ggplot2::geom_segment(
                data      = plot_df,
                ggplot2::aes(x = xstart, xend = xend, y = y, yend = y, text = label),
                color     = color,
                linewidth = 1.2,
                alpha     = 0.8
            )
        }

        gg <- gg + ggplot2::ylab(name)
        if (!log_y)
            gg <- gg + ggplot2::expand_limits(y = 0)

        # for ratio metric, cap axis top at ratio_ymax but zoom in when data is smaller:
        # effective ymax = min(observed data max, ratio_ymax), in the current y-scale.
        # ratio_ymax = 0 → no cap; axis auto-fits (0 as bottom in linear via expand_limits above).
        if (metric == "rna_over_dna" && !is.null(profile$ratio_ymax) &&
            is.finite(profile$ratio_ymax) && profile$ratio_ymax > 0 &&
            is.finite(ratio_data_max)) {
            cap    <- profile$ratio_ymax
            eff    <- min(ratio_data_max, cap)
            ymax_p <- if (log_y) log2(eff) else eff
            ymin_p <- if (log_y) NA else 0
            gg <- gg + ggplot2::coord_cartesian(ylim = c(ymin_p, ymax_p))
        }

        list(plot = gg, legends = list())
    }

    profile_create(
        id               = id,
        name             = name,
        type             = "rna_gene",
        height           = height,
        attr             = list(),
        params           = params,
        plot_f           = plot_f,
        get_gene_table_f = get_gene_table_f,
        get_lib_stats_f  = get_lib_stats_f,
        sample_types     = sample_types,
        auto_register    = auto_register
    )
}

# extract numeric metric vector for one sample_type from gene table
.gene_metric <- function(genes, st, metric, pseudo, stats) {
    if (metric == "rna_over_dna") {
        rna_col <- paste0(st, "_rna_count")
        dna_col <- paste0(st, "_dna_count")
        if (!rna_col %in% names(genes) || !dna_col %in% names(genes)) return(NULL)

        total_rna <- .get_total_mapped(stats, st, "rna")
        total_dna <- .get_total_mapped(stats, st, "dna")

        # add pseudo to raw counts (before normalization); no extra pseudo in ratio
        rna_norm <- (genes[[rna_col]] + pseudo) /
                    (if (!is.na(total_rna) && total_rna > 0) total_rna else 1)
        dna_norm <- (genes[[dna_col]] + pseudo) /
                    (if (!is.na(total_dna) && total_dna > 0) total_dna else 1)
        return(rna_norm / dna_norm)
    }

    # rscore columns are named {st}_rna_rscore in the gene table
    col <- if (metric == "rscore") paste0(st, "_rna_rscore") else paste0(st, "_", metric)
    if (!col %in% names(genes)) return(NULL)
    as.numeric(genes[[col]])
}

