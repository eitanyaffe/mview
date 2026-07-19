########################################################
# RNA coverage profile
#
# Queries srt files live via srt_query (short_R.cpp).
# Shows per-nt depth binned to <= max_points bins.
# Modes: rna | dna | rna+dna | rna_over_dna.
# Multiple timepoints overlaid with distinct colors.
#
# Required by init_short() before use (sources short_R.cpp).
########################################################

RNA_TIMEPOINT_COLORS <- c(
    "D-1"  = "#2171B5",
    "D+4"  = "#D94801",
    "D+7"  = "#238B45",
    "D+18" = "#6A3D9A"
)

# blend a hex color toward white by factor (0 = original, 1 = white)
.lighten_color <- function(col, factor = 0.55) {
    v <- col2rgb(col) / 255
    rgb(v[1] + (1 - v[1]) * factor,
        v[2] + (1 - v[2]) * factor,
        v[3] + (1 - v[3]) * factor)
}

default_rna_cov_params <- function(id, shared_id = id, mode_fixed = NULL) {
    params <- list(
        # shared across RNA/DNA cov and RNA gene profiles
        show_Dm1    = list(group_id = shared_id, type = "boolean", default = TRUE),
        show_Dp4    = list(group_id = shared_id, type = "boolean", default = TRUE),
        show_Dp7    = list(group_id = shared_id, type = "boolean", default = TRUE),
        show_Dp18   = list(group_id = shared_id, type = "boolean", default = TRUE),
        pseudocount = list(group_id = shared_id, type = "double",  default = 0.1),
        ratio_ymax  = list(group_id = shared_id, type = "double",  default = 3),
        map_style   = list(group_id = shared_id, type = "select",
                           choices  = c("non_unique_loose", "non_unique_strict", "unique_strict"),
                           default  = "non_unique_loose"),
        # coverage-profile-specific
        normalize   = list(group_id = id, type = "boolean", default = TRUE),
        max_points  = list(group_id = id, type = "integer", default = 400),
        y_scale     = list(group_id = id, type = "select",
                           choices  = c("log", "linear"), default = "log"),
        height      = list(group_id = id, type = "integer", default = 150)
    )
    # mode selector only when not fixed
    if (is.null(mode_fixed))
        params$mode <- list(group_id = id, type = "select",
                            choices  = c("rna", "dna", "rna+dna", "rna_over_dna"),
                            default  = "rna")
    params
}

# query one srt store for the current view, return a data.frame (contig, mid, coverage)
# or NULL on failure. binsize = 0 → per-nt.
.srt_query_view <- function(store, track, intervals, binsize) {
    if (is.null(store) || !inherits(store, "externalptr")) return(NULL)
    result <- tryCatch(
        srt_query(store, intervals, tracks = track, binsize = binsize),
        error = function(e) { cat(sprintf("srt_query error: %s\n", e$message)); NULL }
    )
    if (is.null(result) || nrow(result) == 0) return(NULL)
    result$mid <- as.integer((result$start + result$end) / 2)
    result
}

# look up total_mapped for (sample_type, lib_kind) from lib_stats data.frame
.get_total_mapped <- function(stats, sample_type, lib_kind) {
    if (is.null(stats)) return(NA_real_)
    ix <- stats$SAMPLE_TYPE == sample_type & stats$LIB_KIND == lib_kind
    if (sum(ix) != 1) return(NA_real_)
    as.numeric(stats$total_mapped[ix])
}

rna_cov_profile <- function(id, name,
                             get_rna_srt_f,
                             get_dna_srt_f,
                             get_lib_stats_f,
                             sample_types  = c("D-1", "D+4", "D+7", "D+18"),
                             height        = 150,
                             shared_id     = id,
                             mode_fixed    = NULL,
                             auto_register = TRUE) {

    params <- default_rna_cov_params(id, shared_id, mode_fixed)

    plot_f <- function(profile, gg) {
        assembly <- cxt_get_assembly()
        xlim     <- cxt_get_xlim()

        if (is.null(xlim) || length(xlim) != 2)
            return(list(plot = gg, legends = list()))

        intervals <- cxt_get_zoom_view()
        if (is.null(intervals) || nrow(intervals) == 0)
            return(list(plot = gg, legends = list()))

        # keep only contig/start/end for srt_query
        iv <- intervals[, c("contig", "start", "end")]

        range_bp <- xlim[2] - xlim[1] + 1
        binsize  <- max(1L, as.integer(ceiling(range_bp / profile$max_points)))

        mode       <- if (!is.null(profile$mode_fixed)) profile$mode_fixed else profile$mode
        normalize  <- profile$normalize
        pseudo     <- profile$pseudocount
        y_scale    <- if (is.null(profile$y_scale)) "log" else profile$y_scale
        map_style  <- if (!is.null(profile$map_style)) profile$map_style else "non_unique_loose"

        # determine which sample_types to show
        show_flags <- c(
            "D-1"  = isTRUE(profile$show_Dm1),
            "D+4"  = isTRUE(profile$show_Dp4),
            "D+7"  = isTRUE(profile$show_Dp7),
            "D+18" = isTRUE(profile$show_Dp18)
        )
        active_types <- intersect(sample_types, names(show_flags)[show_flags])
        if (length(active_types) == 0)
            return(list(plot = gg, legends = list()))

        stats <- get_lib_stats_f(assembly)

        ratio_data_max <- -Inf

        for (st in active_types) {
            color <- RNA_TIMEPOINT_COLORS[st]
            if (is.na(color)) color <- "#888888"

            need_rna <- mode %in% c("rna", "rna+dna", "rna_over_dna")
            need_dna <- mode %in% c("dna", "rna+dna", "rna_over_dna")

            rna_df <- NULL
            dna_df <- NULL

            if (need_rna) {
                rna_store <- get_rna_srt_f(assembly, st, map_style)
                rna_df    <- .srt_query_view(rna_store, "rna", iv, binsize)
            }
            if (need_dna) {
                dna_store <- get_dna_srt_f(assembly, st, map_style)
                dna_df    <- .srt_query_view(dna_store, "dna", iv, binsize)
            }

            total_rna <- .get_total_mapped(stats, st, "rna")
            total_dna <- .get_total_mapped(stats, st, "dna")

            # build one plot_df per mode branch
            if (mode == "rna" && !is.null(rna_df)) {
                cov <- rna_df$coverage
                if (normalize && !is.na(total_rna) && total_rna > 0)
                    cov <- cov / total_rna * 1e6
                plot_df <- data.frame(contig = rna_df$contig,
                                      coord  = rna_df$mid,
                                      value  = cov,
                                      stringsAsFactors = FALSE)
                gg <- .add_cov_layer(gg, plot_df, color, st, y_scale, pseudo)

            } else if (mode == "dna" && !is.null(dna_df)) {
                cov <- dna_df$coverage
                if (normalize && !is.na(total_dna) && total_dna > 0)
                    cov <- cov / total_dna * 1e6
                plot_df <- data.frame(contig = dna_df$contig,
                                      coord  = dna_df$mid,
                                      value  = cov,
                                      stringsAsFactors = FALSE)
                gg <- .add_cov_layer(gg, plot_df, color, paste0(st, " dna"), y_scale, pseudo)

            } else if (mode == "rna+dna") {
                if (!is.null(rna_df)) {
                    cov <- rna_df$coverage
                    if (normalize && !is.na(total_rna) && total_rna > 0)
                        cov <- cov / total_rna * 1e6
                    plot_df <- data.frame(contig = rna_df$contig,
                                          coord  = rna_df$mid,
                                          value  = cov,
                                          stringsAsFactors = FALSE)
                    gg <- .add_cov_layer(gg, plot_df, color, paste0(st, " rna"), y_scale, pseudo)
                }
                if (!is.null(dna_df)) {
                    cov <- dna_df$coverage
                    if (normalize && !is.na(total_dna) && total_dna > 0)
                        cov <- cov / total_dna * 1e6
                    plot_df <- data.frame(contig = dna_df$contig,
                                          coord  = dna_df$mid,
                                          value  = cov,
                                          stringsAsFactors = FALSE)
                    gg <- .add_cov_layer(gg, plot_df, .lighten_color(color), paste0(st, " dna"),
                                         y_scale, pseudo)
                }

            } else if (mode == "rna_over_dna" && !is.null(rna_df) && !is.null(dna_df)) {
                # align bins by mid position (inner join)
                merged <- merge(
                    rna_df[, c("contig", "mid", "coverage")],
                    dna_df[, c("contig", "mid", "coverage")],
                    by = c("contig", "mid"), suffixes = c("_rna", "_dna")
                )
                if (nrow(merged) == 0) next

                # pseudo added to raw counts before normalization; ratio has no additional pseudo
                rna_norm <- (merged$coverage_rna + pseudo) / (if (!is.na(total_rna) && total_rna > 0) total_rna else 1)
                dna_norm <- (merged$coverage_dna + pseudo) / (if (!is.na(total_dna) && total_dna > 0) total_dna else 1)
                ratio    <- rna_norm / dna_norm
                ratio_data_max <- max(ratio_data_max, max(ratio, na.rm = TRUE))
                plot_df  <- data.frame(contig = merged$contig,
                                       coord  = merged$mid,
                                       value  = ratio,
                                       stringsAsFactors = FALSE)
                # ratio is always > 0 (pseudo baked into counts); no pseudo needed for log
                gg <- .add_cov_layer(gg, plot_df, color, paste0(st, " ratio"), y_scale, 0)
            }
        }

        ylab <- if (y_scale == "log") paste0("log10(", name, "+", pseudo, ")") else name
        gg <- gg +
            ggplot2::ylab(ylab) +
            ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "gray85", linewidth = 0.3),
                           panel.grid.major.y = ggplot2::element_line(color = "gray85", linewidth = 0.3),
                           panel.grid.minor.x = ggplot2::element_blank(),
                           panel.grid.minor.y = ggplot2::element_blank())
        if (y_scale == "linear")
            gg <- gg + ggplot2::expand_limits(y = 0)

        # for ratio mode, cap axis top at ratio_ymax but zoom in when data is smaller:
        # effective ymax = min(observed data max, ratio_ymax), in the current y-scale.
        # ratio_ymax = 0 → no cap; axis auto-fits (0 as bottom in linear via expand_limits above).
        if (mode == "rna_over_dna" && !is.null(profile$ratio_ymax) &&
            is.finite(profile$ratio_ymax) && profile$ratio_ymax > 0 &&
            is.finite(ratio_data_max)) {
            cap    <- profile$ratio_ymax
            eff    <- min(ratio_data_max, cap)
            ymax_p <- if (y_scale == "log") log10(eff) else eff
            ymin_p <- if (y_scale == "log") NA else 0
            gg <- gg + ggplot2::coord_cartesian(ylim = c(ymin_p, ymax_p))
        }
        list(plot = gg, legends = list())
    }

    profile_create(
        id              = id,
        name            = name,
        type            = "rna_cov",
        height          = height,
        attr            = list(),
        params          = params,
        plot_f          = plot_f,
        get_rna_srt_f   = get_rna_srt_f,
        get_dna_srt_f   = get_dna_srt_f,
        get_lib_stats_f = get_lib_stats_f,
        sample_types    = sample_types,
        mode_fixed      = mode_fixed,
        auto_register   = auto_register
    )
}

# add a geom_line layer for one (sample_type, mode) series
.add_cov_layer <- function(gg, plot_df, color, label, y_scale, pseudo = 0,
                            lwd = 0.6, lty = 1) {
    mapped <- cxt_filter_coords(plot_df)
    if (is.null(mapped) || nrow(mapped) == 0) return(gg)

    mapped <- mapped[order(mapped$vcoord), ]
    mapped$value <- as.numeric(mapped$value)

    if (y_scale == "log")
        mapped$value <- log10(mapped$value + pseudo)

    # assign a group id per contiguous interval so geom_line doesn't bridge gaps
    vc_diff  <- c(0, diff(mapped$vcoord))
    typical  <- median(vc_diff[vc_diff > 0])
    group_id <- cumsum(vc_diff > 2 * typical) + 1L

    plot_out <- data.frame(
        x     = mapped$vcoord,
        y     = mapped$value,
        group = group_id,
        label = paste0(label, ": ", round(mapped$value, 3)),
        stringsAsFactors = FALSE
    )

    gg + ggplot2::geom_line(
        data      = plot_out,
        ggplot2::aes(x = x, y = y, group = group, text = label),
        color     = color,
        linewidth = lwd,
        linetype  = lty
    )
}
