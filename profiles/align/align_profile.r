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

# Get variant colors (shared between full and pileup modes)
get_variant_colors <- function(aln, cxt) {
  # query both full and pileup to get all possible variants
  all_variants <- character(0)
  
  # get variants from pileup data
  tryCatch({
    pileup_df <- aln_query_pileup(aln, cxt$intervals, report_mode_str = "all")
    if (!is.null(pileup_df) && nrow(pileup_df) > 0) {
      all_variants <- c(all_variants, unique(pileup_df$variant))
    }
  }, error = function(e) {
    # continue if pileup query fails
  })
  
  # get variants from full data mutations (if any)
  tryCatch({
    full_df <- aln_query_full(aln, cxt$intervals, "by_mutations", 1000, "all")
    if (!is.null(full_df$mutations) && nrow(full_df$mutations) > 0) {
      all_variants <- c(all_variants, unique(full_df$mutations$desc))
    }
  }, error = function(e) {
    # continue if full query fails
  })
  
  # get unique sorted variants
  unique_variants <- sort(unique(all_variants))
  
  if (length(unique_variants) == 0) {
    return(list())
  }
  
  # create color palette
  if (length(unique_variants) <= 8) {
    # use predefined soft colors for small numbers of variants
    soft_colors <- c("#E8F4FD", "#FFE6E6", "#E6F7E6", "#FFF2E6", 
                     "#F0E6FF", "#E6F7FF", "#FFE6F7", "#F7FFE6")
    variant_colors <- soft_colors[seq_along(unique_variants)]
  } else {
    # use soft pastel colors for larger numbers
    variant_colors <- rainbow(length(unique_variants), s = 0.3, v = 0.9)
  }
  names(variant_colors) <- unique_variants
  
  # standardized gray for REF
  if ("REF" %in% unique_variants) {
    variant_colors["REF"] <- "#D3D3D3"  # lighter than dark gray, darker than light gray
  }
  
  return(variant_colors)
}
default_alignment_params <- list(
  plot_style = list(
    group_id = "alignment",
    type = "select",
    choices = c("auto_full", "auto_pileup", "bin", "full", "pileup"),
    default = "auto_full"
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
  full_threshold = list(
    group_id = "alignment",
    type = "integer",
    default = 50000
  ),
  pileup_threshold = list(
    group_id = "alignment",
    type = "integer",
    default = 1000
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
  ),
  full_min_mutations_density = list(
    group_id = "alignment",
    type = "double",
    default = 0
  ),
  full_mutation_lwd = list(
    group_id = "alignment", 
    type = "double",
    default = 0.5
  )
)

align_profile <- function(id, name, height = 400,
                          aln_f = NULL,
                          bin_type = "auto",
                          plot_style = "auto_full",

                          full_threshold = 50000,
                          pileup_threshold = 1000,
                          target_bins = 250,
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
