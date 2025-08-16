# Helper for creating an alignment-based profile
# Supports bin/pileup/full modes based on zoom level

# define the shared red color palette used across all mutation visualizations
get_red_palette <- function() {
  return(c("#fee0d2", "#fc9272", "#fb6a4a", "#cb181d", "#a50f15", "#67000d"))
}

# shared red color scale function for both continuous and discrete red colors
get_shared_red_scale <- function(continuous = TRUE, num_colors = 6) {
  red_palette <- get_red_palette()
  
  if (continuous) {
    # return the base colors for continuous scale
    return(c(red_palette[1], red_palette[6]))
  } else {
    # return discrete red colors from light to dark
    return(red_palette[seq_len(min(num_colors, length(red_palette)))])
  }
}

# Get mutation density colors (shared between bin and full modes)
get_mutation_colors <- function(mutations_per_bp) {
  # use shared red palette
  red_palette <- get_red_palette()
  
  # apply same thresholds as discrete scale using vectorized ifelse for efficiency
  colors <- ifelse(mutations_per_bp <= 0, red_palette[1],
            ifelse(mutations_per_bp <= 0.0001, red_palette[2],
            ifelse(mutations_per_bp <= 0.001, red_palette[3],
            ifelse(mutations_per_bp <= 0.01, red_palette[4],
            ifelse(mutations_per_bp <= 0.1, red_palette[5],
                   red_palette[6])))))
  
  return(colors)
}


default_alignment_params <- list(
  # align_general - parameters used across multiple modes
  plot_style = list(
    group_id = "align_general",
    type = "select",
    choices = c("auto_full", "auto_pileup", "bin", "full", "pileup"),
    default = "auto_full"
  ),
  alignment_filter = list(
    group_id = "align_general",
    type = "select",
    choices = c("all", "single", "single_complete", "only_multiple"),
    default = "all"
  ),
  height_style = list(
    group_id = "align_general",
    type = "select",
    choices = c("by_mutations", "by_coord_left", "by_coord_right"), 
    default = "by_mutations"
  ),

  # align_full - parameters specific to full profile mode
  full_style = list(
    group_id = "align_full",
    type = "select",
    choices = c("none", "by_mutations", "by_strand", "show_mutations"), 
    default = "by_mutations"
  ),
  full_threshold = list(
    group_id = "align_full",
    type = "integer",
    default = 50000
  ),
  max_reads = list(
    group_id = "align_full",
    type = "integer",
    default = 1000
  ),
  max_mutations = list(
    group_id = "align_full",
    type = "integer",
    default = 1000
  ),
  full_min_mutations_density = list(
    group_id = "align_full",
    type = "double",
    default = 0
  ),
  full_mutation_lwd = list(
    group_id = "align_full", 
    type = "double",
    default = 0.5
  ),

  # align_bin / align_pileup - parameters for bin and pileup modes
  target_bins = list(
    group_id = "align_bin",
    type = "integer",
    default = 250
  ),
  bin_style = list(
    group_id = "align_bin",
    type = "select",
    choices = c("by_mut_density", "by_genomic_distance", "by_seg_density", "by_nonref_density"),
    default = "by_mut_density"
  ),
  bin_type = list(
    group_id = "align_bin",
    type = "select",
    choices = c("auto", "10", "100", "1000", "5000", "10000"),
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
  normalize_distrib_bins = list(
    group_id = "align_bin",
    type = "boolean",
    default = FALSE
  ),
  pileup_threshold = list(
    group_id = "align_pileup",
    type = "integer",
    default = 1000
  )
)

align_profile <- function(id, name, height = 400,
                          aln_f = NULL,
                          bin_type = "auto",
                          plot_style = "auto_full",
                          full_threshold = 50000,
                          pileup_threshold = 1000,
                          target_bins = 100,
                          height_style = "by_mutations",
                          max_mutations = 1000,
                          max_reads = 1000,
                          full_style = "by_mutations",
                          alignment_filter = "all",
                          full_min_mutations_density = 0,
                          full_mutation_lwd = 0.5,
                          params = default_alignment_params,
                          auto_register = TRUE) {
  # Check for required alignment functions
  required_funcs <- c("aln_query_full", "aln_query_pileup", "aln_query_bin")
  missing_funcs <- required_funcs[!sapply(required_funcs, exists, mode = "function")]
  if (length(missing_funcs) > 0) {
    stop(paste("Missing functions:", paste(missing_funcs, collapse = ", ")))
  }

  # Determine display mode based on view range
  get_display_mode <- function(xlim, full_threshold, pileup_threshold, plot_style) {
    if (is.null(xlim) || length(xlim) != 2) {
      return("bin")
    }

    range_bp <- (xlim[2] + 1) - xlim[1]

    if (plot_style == "auto_full") {
      return(if (range_bp <= full_threshold) "full" else "bin")
    } else if (plot_style == "auto_pileup") {
      return(if (range_bp <= pileup_threshold) "pileup" else "bin")
    } else {
      # explicit mode: bin, full, or pileup
      return(plot_style)
    }
  }

  plot_f <- function(profile, cxt, gg) {
    # Check that intervals dataframe has at least one row
    if (is.null(cxt$intervals) || nrow(cxt$intervals) == 0) {
      warning("intervals dataframe is empty")
      return(gg)
    }

    aln <- if (is.function(aln_f)) aln_f(cxt) else aln_f

    # !!! we have multiple alignments, we need to get the correct one
    cache_set("aln_obj", aln)
    
  if (any(!is.element(cxt$contigs, aln_get_contigs(aln)$contig_id))) {
    warning("contigs not found in alignment")
    return(gg)
  }


    if (is.null(aln) || !inherits(aln, "externalptr")) {
      warning(sprintf("align_profile '%s': Invalid AlignmentStore pointer.", name))
      return(gg)
    }

    mode <- get_display_mode(cxt$mapper$xlim, profile$full_threshold, profile$pileup_threshold, profile$plot_style)
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
    id = id, name = name, type = "align", height = height,
    params = params, plot_f = plot_f,
    auto_register = auto_register,
    aln_f = aln_f,
    bin_type = bin_type,
    plot_style = plot_style,
    full_threshold = full_threshold,
    pileup_threshold = pileup_threshold,
    target_bins = target_bins,
    height_style = height_style,
    max_mutations = max_mutations,
    max_reads = max_reads,
    full_style = full_style,
    alignment_filter = alignment_filter,
    full_min_mutations_density = full_min_mutations_density,
    full_mutation_lwd = full_mutation_lwd
  )
}
