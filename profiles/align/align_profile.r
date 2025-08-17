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

# Get fixed colors for variant types (shared between full and pileup modes)
get_variant_type_colors <- function(variants) {
  # define gentle but visible colors for all 12 nucleotide substitutions
  substitution_colors <- c(
    "A:T" = "#e6d19a", "A:C" = "#b3e6f2", "A:G" = "#d9b3ff",
    "T:A" = "#ffe0b3", "T:C" = "#b3e6d1", "T:G" = "#b3d9ff",
    "C:A" = "#ffb3d9", "C:T" = "#d9ffb3", "C:G" = "#ffb3ff",
    "G:A" = "#ffecb3", "G:T" = "#b3e0ff", "G:C" = "#ffccb3"
  )
  
  # colors for gains (3 red shades, darker range, squeezed)
  gain_light <- "#cc6666"   # darker light red for +1
  gain_medium <- "#b30000"  # darker medium red for +2 to +3
  gain_dark <- "#800000"    # darkest red for +4 or more
  
  # colors for losses (3 blue shades, darker range, squeezed)  
  loss_light <- "#6666cc"   # darker light blue for -1
  loss_medium <- "#0000b3"  # darker medium blue for -2 to -3
  loss_dark <- "#000080"    # darkest blue for -4 or more
  
  # ref and unknown colors
  ref_color <- "#e0dede"   # light gray
  unknown_color <- "#000000" # black
  
  # assign colors based on variant pattern
  colors <- character(length(variants))
  
  for (i in seq_along(variants)) {
    variant <- variants[i]
    
    if (variant == "REF") {
      colors[i] <- ref_color
    } else if (variant %in% names(substitution_colors)) {
      # direct substitution match
      colors[i] <- substitution_colors[variant]
    } else if (grepl("^\\+", variant)) {
      # gain pattern: +XXX
      gain_length <- nchar(gsub("^\\+", "", variant))
      colors[i] <- if (gain_length == 1) {
        gain_light
      } else if (gain_length <= 3) {
        gain_medium
      } else {
        gain_dark
      }
    } else if (grepl("^-", variant)) {
      # loss pattern: -XXX
      loss_length <- nchar(gsub("^-", "", variant))
      colors[i] <- if (loss_length == 1) {
        loss_light
      } else if (loss_length <= 3) {
        loss_medium
      } else {
        loss_dark
      }
    } else {
      # unknown variant type
      colors[i] <- unknown_color
    }
  }
  
  return(colors)
}

# get simplified colors for variant types (type mode: sub=orange, deletion=blue, addition=red)
get_mutation_type_colors <- function(variants) {
  # define simplified color scheme
  substitution_color <- "#ff8c00"  # orange for all substitutions
  deletion_color <- "#0066cc"      # blue for all deletions  
  addition_color <- "#cc0000"      # red for all additions
  ref_color <- "#e0dede"           # light gray
  unknown_color <- "#000000"       # black
  
  # assign colors based on variant pattern
  colors <- character(length(variants))
  
  for (i in seq_along(variants)) {
    variant <- variants[i]
    
    if (variant == "REF") {
      colors[i] <- ref_color
    } else if (grepl(":", variant)) {
      # substitution pattern: X:Y
      colors[i] <- substitution_color
    } else if (grepl("^\\+", variant)) {
      # addition pattern: +XXX
      colors[i] <- addition_color
    } else if (grepl("^-", variant)) {
      # deletion pattern: -XXX
      colors[i] <- deletion_color
    } else {
      # unknown variant type
      colors[i] <- unknown_color
    }
  }
  
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
  mutation_color_mode = list(
    group_id = "align_general",
    type = "select", 
    choices = c("detailed", "type"),
    default = "detailed"
  ),
  clip_mode = list(
    group_id = "align_filter",
    type = "select",
    choices = c("all", "complete", "allow_one_side_clip", "only_one_side_clipped", "only_two_side_clipped"),
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
  height_style = list(
    group_id = "align_general",
    type = "select",
    choices = c("by_mutations", "by_coord_left", "by_coord_right"), 
    default = "by_mutations"
  ),
  height = list(
    group_id = "align_general",
    type = "integer",
    default = 400
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

  full_mutation_lwd = list(
    group_id = "align_full", 
    type = "double",
    default = 1.5
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

align_profile <- function(id, name,
                          aln_f = NULL,
                          bin_type = "auto",
                          plot_style = "auto_full",
                          mutation_color_mode = "detailed",
                          full_threshold = 50000,
                          pileup_threshold = 1000,
                          target_bins = 100,
                          height_style = "by_mutations",
                          max_mutations = 1000,
                          max_reads = 1000,
                          full_style = "by_mutations",
                          clip_mode = "all",
                          clip_margin = 10,
                          min_mutations_percent = 0.0,
                          max_mutations_percent = 10.0,
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
    id = id, name = name, type = "align", height = params$height$default,
    params = params, plot_f = plot_f,
    auto_register = auto_register,
    aln_f = aln_f,
    bin_type = bin_type,
    plot_style = plot_style,
    mutation_color_mode = mutation_color_mode,
    full_threshold = full_threshold,
    pileup_threshold = pileup_threshold,
    target_bins = target_bins,
    height_style = height_style,
    max_mutations = max_mutations,
    max_reads = max_reads,
    full_style = full_style,
    clip_mode = clip_mode,
    clip_margin = clip_margin,
    min_mutations_percent = min_mutations_percent,
    max_mutations_percent = max_mutations_percent,
    full_mutation_lwd = full_mutation_lwd
  )
}
