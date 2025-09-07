# Helper for creating an alignment-based profile
# Supports bin/pileup/full modes based on zoom level

# Load shared utilities for mutation coloring and legends
source("profiles/align/align_utils.r")

default_alignment_params <- list(
  # align_general - parameters used across multiple modes
  plot_style = list(
    group_id = "align_general",
    type = "select",
    choices = c("auto", "bin", "full"),
    default = "auto"
  ),
  mutation_color_mode = list(
    group_id = "mutations",
    type = "select", 
    choices = c("detailed", "type"),
    default = "detailed"
  ),
  full_threshold = list(
    group_id = "mutations",
    type = "integer",
    default = 200000
  ),
  clip_mode = list(
    group_id = "align_filter",
    type = "select",
    choices = c("all", "complete", "allow_one_side_clip", "only_one_side_clipped", "only_two_side_clipped", "only_clipped", "local_align"),
    default = "all"
  ),
  clip_margin = list(
    group_id = "align_filter",
    type = "integer",
    default = 10
  ),
  min_mutations_percent = list(
    group_id = "align_filter",
    type = "select",
    choices = c("0%" = 0.0, "0.01%" = 0.01, "0.1%" = 0.1, "1%" = 1.0, "10%" = 10.0),
    default = 0.0
  ),
  max_mutations_percent = list(
    group_id = "align_filter",
    type = "select",
    choices = c("0%" = 0.0, "0.01%" = 0.01, "0.1%" = 0.1, "1%" = 1.0, "10%" = 10.0),
    default = 10.0
  ),
  min_alignment_length = list(
    group_id = "align_filter",
    type = "integer",
    default = 0
  ),
  max_alignment_length = list(
    group_id = "align_filter",
    type = "integer",
    default = 0
  ),
  min_indel_length = list(
    group_id = "align_filter",
    type = "integer",
    default = 3
  ),
  height = list(
    group_id = "align_general",
    type = "integer",
    default = 400
  ),
  force_max_y = list(
    group_id = "align_general",
    type = "integer",
    default = 0
  ),

  # align_full - parameters specific to full profile mode
  height_style = list(
    group_id = "align_full",
    type = "select",
    choices = c("by_mutations", "by_coord_left", "by_coord_right"), 
    default = "by_mutations"
  ),  
  full_style = list(
    group_id = "align_full",
    type = "select",
    choices = c("none", "by_mutations", "by_strand", "show_mutations"), 
    default = "show_mutations"
  ),

  max_reads = list(
    group_id = "align_full",
    type = "integer",
    default = 100000
  ),
  max_mutations = list(
    group_id = "align_full",
    type = "integer",
    default = 10000
  ),

  full_mutation_lwd = list(
    group_id = "align_full", 
    type = "double",
    default = 1
  ),
  bin_style = list(
    group_id = "align_bin",
    type = "select",
    choices = c("by_seg_density", "by_mut_density", "by_median_mutation_density", "by_genomic_distance", "by_nonref_density", "by_seg_clip_density", "by_non_ref_clip_density"),
    default = "by_seg_density"
  ),
  bin_type = list(
    group_id = "align_bin",
    type = "select",
    choices = c("auto", "500", "1000", "2000", "5000", "10000", "20000", "50000", "100000", "200000", "500000"),
    default = "auto"
  ),
  seg_threshold = list(
    group_id = "align_bin",
    type = "double",
    default = 0.2
  ),
  non_ref_threshold = list(
    group_id = "align_bin",
    type = "double",
    default = 0.9
  ),
  num_threads = list(
    group_id = "align_bin",
    type = "integer",
    default = 0
  ),
  target_bins = list(
    group_id = "align_bin",
    type = "integer",
    default = 250
  ),
  normalize_distrib_bins = list(
    group_id = "align_bin",
    type = "boolean",
    default = FALSE
  ),
  pileup_threshold = list(
    group_id = "align_pileup",
    type = "integer",
    default = 200
  ),
  max_col_dist_percent = list(
    group_id = "align_general",
    type = "select",
    choices = c("auto", "0.1", "1", "10"),
    default = "auto"
  ),
  use_pileup = list(
    group_id = "align_general",
    type = "boolean",
    default = FALSE
  ),
  show_hover = list(
    group_id = "align_general",
    type = "boolean",
    default = TRUE
  )
)

align_profile <- function(id, name,
                          aln_f = NULL,
                          bin_type = "auto",
                          plot_style = "auto",
                          mutation_color_mode = "detailed",
                          full_threshold = 100000,
                          pileup_threshold = 200,
                          use_pileup = FALSE,
                          target_bins = 100,
                          height_style = "by_mutations",
                          max_mutations = 1000,
                          max_reads = 1000,
                          full_style = "by_mutations",
                          clip_mode = "all",
                          clip_margin = 10,
                          min_mutations_percent = 0.0,
                          max_mutations_percent = 10.0,
                          min_alignment_length = 0,
                          max_alignment_length = 0,
                          min_indel_length = 3,
                          full_mutation_lwd = 0.5,
                          force_max_y = 0,
                          show_hover = TRUE,
                          binsizes = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000),
                          params = default_alignment_params,
                          auto_register = TRUE) {
  # Check for required alignment functions
  required_funcs <- c("aln_query_full", "aln_query_pileup", "aln_query_bin")
  missing_funcs <- required_funcs[!sapply(required_funcs, exists, mode = "function")]
  if (length(missing_funcs) > 0) {
    stop(paste("Missing functions:", paste(missing_funcs, collapse = ", ")))
  }

  # Determine display mode based on view range
  get_display_mode <- function(xlim, full_threshold, pileup_threshold, plot_style, use_pileup) {
    if (is.null(xlim) || length(xlim) != 2) {
      return("bin")
    }

    range_bp <- (xlim[2] + 1) - xlim[1]

    if (plot_style == "auto") {
      if (use_pileup && range_bp <= pileup_threshold) {
        return("pileup")
      } else if (range_bp <= (full_threshold + 1)) {
        return("full")
      } else {
        return("bin")
      }
    } else {
      # explicit mode: bin, full
      return(plot_style)
    }
  }

  plot_f <- function(profile, cxt, gg) {
    # Check that intervals dataframe has at least one row
    if (is.null(cxt$intervals) || nrow(cxt$intervals) == 0) {
      warning("intervals dataframe is empty")
      return(list(plot = gg, legends = list()))
    }

    aln <- if (is.function(aln_f)) aln_f(cxt) else aln_f

    # !!! we have multiple alignments, we need to get the correct one
    cache_set("aln_obj", aln)
    
  if (any(!is.element(cxt$contigs, aln_get_contigs(aln)$contig_id))) {
    warning("contigs not found in alignment")
    return(list(plot = gg, legends = list()))
  }


    if (is.null(aln) || !inherits(aln, "externalptr")) {
      warning(sprintf("align_profile '%s': Invalid AlignmentStore pointer.", name))
      return(list(plot = gg, legends = list()))
    }

    mode <- get_display_mode(cxt$mapper$xlim, profile$full_threshold, profile$pileup_threshold, profile$plot_style, profile$use_pileup)
    cat(sprintf("mode: %s\n", mode))

    # Call mode-specific plot function
    if (mode == "full") {
      return(align_profile_full(profile, cxt, aln, gg))
    } else if (mode == "pileup") {
      return(align_profile_pileup(profile, cxt, aln, gg))
    } else { # mode == "bin"
      return(align_profile_bin(profile, cxt, aln, gg))
    }
  }

  # Create profile
  profile_create(
    id = id, name = name, type = "align", height = params$height$default,
    params = params, plot_f = plot_f,
    auto_register = auto_register,
    aln_f = aln_f,
    bin_type = bin_type,
    plot_style = plot_style,
    mutation_color_mode = mutation_color_mode,
    full_threshold = full_threshold,
    pileup_threshold = pileup_threshold,
    use_pileup = use_pileup,
    target_bins = target_bins,
    height_style = height_style,
    max_mutations = max_mutations,
    max_reads = max_reads,
    full_style = full_style,
    clip_mode = clip_mode,
    clip_margin = clip_margin,
    min_mutations_percent = min_mutations_percent,
    max_mutations_percent = max_mutations_percent,
    min_alignment_length = min_alignment_length,
    max_alignment_length = max_alignment_length,
    min_indel_length = min_indel_length,
    full_mutation_lwd = full_mutation_lwd,
    force_max_y = force_max_y,
    show_hover = show_hover,
    binsizes = binsizes
  )
}
