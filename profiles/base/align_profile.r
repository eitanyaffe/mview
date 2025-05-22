# Helper for creating an alignment-based profile
# Supports bin/pileup/full modes based on zoom level

align_profile <- function(id, name, height = 1,
                          aln_f = NULL,
                          bin_type = "auto",
                          thresholds = list(full = 20000, pileup = 20000),
                          params = list(),
                          auto_register = TRUE) {
  # Check for required alignment functions
  required_funcs <- c("aln_query_full", "aln_query_pileup", "aln_query_bin")
  missing_funcs <- required_funcs[!sapply(required_funcs, exists, mode = "function")]
  if (length(missing_funcs) > 0) {
    stop(paste("Missing functions:", paste(missing_funcs, collapse = ", ")))
  }

  # Determine display mode based on view range
  get_display_mode <- function(xlim, thresholds) {
    if (is.null(xlim) || length(xlim) != 2) {
      return("bin")
    }

    range_bp <- (xlim[2] + 1) - xlim[1]

    if (range_bp <= thresholds$full) {
      return("full")
    } else if (range_bp <= thresholds$pileup) {
      return("pileup")
    } else {
      return("bin")
    }
  }

  plot_f <- function(profile, cxt, gg) {
    aln <- if (is.function(aln_f)) aln_f(cxt) else aln_f
    if (is.null(aln) || !inherits(aln, "externalptr")) {
      warning(sprintf("align_profile '%s': Invalid AlignmentStore pointer.", name))
      return(NULL)
    }

    mode <- get_display_mode(cxt$mapper$xlim, profile$thresholds)
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
    aln_f = aln_f,
    bin_type = bin_type,
    thresholds = thresholds,
    auto_register = auto_register
  )
}
