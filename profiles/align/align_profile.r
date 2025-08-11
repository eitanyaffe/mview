# Helper for creating an alignment-based profile
# Supports bin/pileup/full modes based on zoom level

# Get mutation density colors (shared between bin and full modes)
get_mutation_colors <- function(mutations_per_100bp) {
  # use consistent scale: 0 to 1 mutation per 100bp
  colors <- get_color_scale(
    values = mutations_per_100bp,
    colors = c("lightgray", "red"),
    min_val = 0,
    max_val = 1.0,
    num_steps = 20
  )
  return(colors)
}
default_alignment_params <- list(
  plot_style = list(
    group_id = "alignment",
    type = "select",
    choices = c("auto", "full", "bin"),
    default = "auto"
  ),
  alignment_filter = list(
    group_id = "alignment",
    type = "select",
    choices = c("all", "single", "single_complete", "only_multiple"),
    default = "all"
  ),
  full_style = list(
    group_id = "alignment",
    type = "select",
    choices = c("none", "by_mutations", "by_strand", "show_mutations"), 
    default = "by_mutations"
  ),
  target_bins = list(
    group_id = "alignment",
    type = "integer",
    default = 250
  ),
  threshold = list(
    group_id = "alignment",
    type = "integer",
    default = 50000
  ),
  max_reads = list(
    group_id = "alignment",
    type = "integer",
    default = 1000
  ),
  max_mutations = list(
    group_id = "alignment",
    type = "integer",
    default = 1000
  ),
  height_style = list(
    group_id = "alignment",
    type = "select",
    choices = c("by_mutations", "by_coord_left", "by_coord_right"), 
    default = "by_mutations"
  ),
  bin_type = list(
    group_id = "alignment",
    type = "select",
    choices = c("auto", "10", "100", "1000", "5000", "10000"),
    default = "auto"
  )
)

align_profile <- function(id, name, height = 1,
                          aln_f = NULL,
                          bin_type = "auto",
                          plot_style = "auto",
                          threshold = 50000,
                          target_bins = 250,
                          height_style = "by_mutations",
                          max_mutations = 1000,
                          max_reads = 1000,
                          full_style = "by_mutations",
                          alignment_filter = "all",
                          params = default_alignment_params,
                          auto_register = TRUE) {
  # Check for required alignment functions
  required_funcs <- c("aln_query_full", "aln_query_pileup", "aln_query_bin")
  missing_funcs <- required_funcs[!sapply(required_funcs, exists, mode = "function")]
  if (length(missing_funcs) > 0) {
    stop(paste("Missing functions:", paste(missing_funcs, collapse = ", ")))
  }

  # Determine display mode based on view range
  get_display_mode <- function(xlim, threshold, plot_style) {
    if (is.null(xlim) || length(xlim) != 2) {
      return("bin")
    }

    range_bp <- (xlim[2] + 1) - xlim[1]

    if (plot_style != "auto") {
      return(plot_style)
    } else {
      if (range_bp <= threshold) {
        return("full")
      } else {
        return("bin")
      }
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

    mode <- get_display_mode(cxt$mapper$xlim, profile$threshold, profile$plot_style)
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
    threshold = threshold,
    target_bins = target_bins,
    height_style = height_style,
    max_mutations = max_mutations,
    max_reads = max_reads,
    full_style = full_style,
    alignment_filter = alignment_filter
  )
}
