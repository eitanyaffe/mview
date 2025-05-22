# Full mode for alignment profile

align_query_full_mode <- function(aln, cxt) {
  # Query full alignment data
  df <- aln_query_full(aln, cxt$intervals, height_style_str = "by_mutations")

  # Process alignments
  alignments <- NULL
  if (!is.null(df$alignments) && nrow(df$alignments) > 0) {
    alns <- df$alignments
    alns$contig <- alns$contig_id
    alns$start <- alns$contig_start + 1
    alns$end <- alns$contig_end
    alignments <- filter_segments(alns, cxt, cxt$mapper$xlim)
  }

  # Process mutations
  mutations <- NULL
  if (!is.null(df$mutations) && nrow(df$mutations) > 0) {
    muts <- df$mutations
    muts$contig <- muts$contig_id
    muts$coord <- muts$position + 1
    mutations <- filter_coords(muts, cxt, cxt$mapper$xlim)
  }

  return(list(alignments = alignments, mutations = mutations))
}

align_profile_full <- function(profile, cxt, aln, gg) {
  df <- align_query_full_mode(aln, cxt)

  # Plot alignments
  if (!is.null(df$alignments) && nrow(df$alignments) > 0) {
    cat(sprintf("plotting %d alignments\n", nrow(df$alignments)))
    alignments <- df$alignments
    alignments$rect_ymin <- alignments$height
    alignments$rect_ymax <- alignments$height + 1
    # clip gstart and gend to xlim
    alignments$gstart <- pmax(alignments$gstart, cxt$mapper$xlim[1])
    alignments$gend <- pmin(alignments$gend, cxt$mapper$xlim[2])

    alignments$hover_text <- paste0(
      "Read: ", alignments$read_id, "\n",
      alignments$start, "-", alignments$end, "\n",
      "Mutations: ", alignments$mutation_count
    )

    gg <- gg + ggplot2::geom_rect(
      data = alignments,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = rect_ymin, ymax = rect_ymax,
        text = hover_text
      ),
      fill = "lightgray", color = NA
    )
  }

  # Plot mutations
  if (!is.null(df$mutations) && nrow(df$mutations) > 0) {
    cat(sprintf("plotting %d mutations\n", nrow(df$mutations)))
    mutations <- df$mutations

    # Sample if too many mutations
    if (nrow(mutations) > 10000) {
      mutations <- mutations[sample(nrow(mutations), 10000), ]
      cat("sampled to 10000 mutations\n")
    }

    # Add hover text for mutations
    mutations$hover_text <- paste0(
      "Read: ", mutations$read_id, "\n",
      "Position: ", mutations$coord, "\n",
      "Type: ", mutations$desc
    )

    # Plot vertical segments of length 1 instead of points
    gg <- gg + ggplot2::geom_segment(
      data = mutations,
      ggplot2::aes(
        x = gcoord, xend = gcoord,
        y = height, yend = height + 1,
        color = desc,
        text = hover_text
      ),
      size = 0.5
    ) +
      ggplot2::scale_color_discrete(name = "Mutation") +
      ggplot2::theme(legend.position = "top")
  }

  return(gg)
}
