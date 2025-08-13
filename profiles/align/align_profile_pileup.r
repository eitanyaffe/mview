# Pileup mode for alignment profile

align_query_pileup_mode <- function(aln, cxt) {
  # Query pileup data
  df <- aln_query_pileup(aln, cxt$intervals, report_mode_str = "all")

  if (!is.null(df) && nrow(df) > 0) {
    df$coord <- df$position + 1
    return(filter_coords(df, cxt, cxt$mapper$xlim))
  }

  return(NULL)
}

align_profile_pileup <- function(profile, cxt, aln, gg) {
  # Pileup plotting implementation
  # (Currently a placeholder as noted in original code)
  df <- align_query_pileup_mode(aln, cxt)
  if (is.null(df) || nrow(df) == 0) {
    return(gg)
  }
  # baseline at y=0
  gg <- gg + ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
  return(gg)
}
