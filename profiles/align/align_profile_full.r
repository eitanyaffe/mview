# Full mode for alignment profile

# Get colors for alignments based on color mode
get_alignment_colors <- function(alignments, reads, style) {
  if (style == "none") {
    return(rep("lightgray", nrow(alignments)))
  } else if (style == "by_strand") {
    # get corresponding read info
    read_lookup <- setNames(reads$is_reverse, reads$read_id)
    read_is_reverse <- read_lookup[alignments$read_id]
    
    # xor of read and alignment strand
    alignments$along_contig = !xor(alignments$is_reverse, read_is_reverse)
    colors <- ifelse(alignments$along_contig, "lightgreen", "#f5d34d")
    return(colors)
  } else if (style == "by_mutations") {
    # use shared mutation color function
    alignment_lengths <- alignments$end - alignments$start + 1
    mutations_per_100bp <- (alignments$mutation_count / alignment_lengths) * 100
    colors <- get_mutation_colors(mutations_per_100bp)
    return(colors)
  } else if (style == "show_mutations") {
    return(rep("lightgray", nrow(alignments)))
  } else {
    return(rep("lightgray", nrow(alignments)))
  }
}

align_query_full_mode <- function(aln, cxt, 
height_style_str = "by_mutations", 
max_reads = 1000,
alignment_filter = "all") {
  # Query full alignment data
  df <- aln_query_full(aln, cxt$intervals, height_style_str, max_reads, alignment_filter)
  cat(sprintf("done with aln_query_full\n"))
  
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

  # Process reads
  reads <- NULL
  if (!is.null(df$reads) && nrow(df$reads) > 0) {
    reads <- df$reads
    reads$contig <- reads$contig_id
    reads$start <- reads$span_start + 1
    reads$end <- reads$span_end
    reads <- filter_segments(reads, cxt, cxt$mapper$xlim)
  }

  return(list(alignments = alignments, mutations = mutations, reads = reads))
}

align_profile_full <- function(profile, cxt, aln, gg) {
  height_style <- profile$height_style
  if (!is.element(height_style, c("by_mutations", "by_coord_left", "by_coord_right"))) {
    height_style <- "by_mutations"
  }

  df <- align_query_full_mode(
    aln = aln, 
    cxt = cxt, 
    height_style = height_style, 
    max_reads = profile$max_reads, 
    alignment_filter = profile$alignment_filter)
  
  # Plot reads first (under alignments)
  if (!is.null(df$reads) && nrow(df$reads) > 0) {

    cat(sprintf("plotting %d reads\n", nrow(df$reads)))
    reads <- df$reads
    # clip gstart and gend to xlim
    reads$gstart <- pmax(reads$gstart, cxt$mapper$xlim[1])
    reads$gend <- pmin(reads$gend, cxt$mapper$xlim[2])

    reads$hover_text <- paste0(
      "Read: ", reads$read_id, "\n",
      reads$start, "-", reads$end
    )

    # draw read lines at middle height
    gg <- gg + ggplot2::geom_segment(
      data = reads,
      ggplot2::aes(
        x = gstart, xend = gend,
        y = height + 0.5, yend = height + 0.5,
        text = hover_text,
        key = read_id
      ),
      color = "black", size = 0.5
    )
  }

  # Plot alignments
  if (!is.null(df$alignments) && nrow(df$alignments) > 0) {
    cat(sprintf("plotting %d alignments\n", nrow(df$alignments)))

    alignments <- df$alignments
    alignments = alignments[is.element(alignments$read_id, reads$read_id), ]
    alignments$rect_ymin <- alignments$height
    alignments$rect_ymax <- alignments$height + 1

    # determine clipped status
    left_in_plot = alignments$gstart > cxt$mapper$xlim[1]
    right_in_plot = alignments$gend < cxt$mapper$xlim[2]

    clip_read_start = alignments$read_start > 0
    clip_read_end = alignments$read_end < alignments$read_length

    alignments$clipped_left <- ifelse(!alignments$is_reverse, clip_read_start, clip_read_end) & left_in_plot
    alignments$clipped_right <- ifelse(!alignments$is_reverse, clip_read_end, clip_read_start) & right_in_plot
    
    # clip gstart and gend to xlim
    alignments$gstart <- pmax(alignments$gstart, cxt$mapper$xlim[1])
    alignments$gend <- pmin(alignments$gend, cxt$mapper$xlim[2])

    # set colors based on color mode
    alignments$color <- get_alignment_colors(alignments, reads, profile$full_style)

    alignments$align_length <- alignments$end - alignments$start + 1
    alignments$mut_density <- alignments$mutation_count / (alignments$align_length / 1000)
    alignments$hover_text <- paste0(
      "Read: ", alignments$read_id, "\n",
      alignments$start, "-", alignments$end, "\n",
      "Read coords: ", alignments$read_start, "-", alignments$read_end, " / ", alignments$read_length, "\n",
      "Mutations: ", alignments$mutation_count, " (", sprintf("%.2f", alignments$mut_density), " muts/1kb)"
    )

    # draw alignment rectangles
    gg <- gg + ggplot2::geom_rect(
      data = alignments,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = rect_ymin, ymax = rect_ymax,
        text = hover_text,
        key = read_id
      ),
      fill = alignments$color, color = NA
    )

    # draw left clipping indicators
    if (any(alignments$clipped_left)) {
      left_clipped <- alignments[alignments$clipped_left, ]
      gg <- gg + ggplot2::geom_segment(
        data = left_clipped,
        ggplot2::aes(
          x = gstart, xend = gstart,
          y = rect_ymin, yend = rect_ymax,
          text = hover_text,
          key = read_id
        ),
        color = "black", size = 0.5
      )
    }

    # draw right clipping indicators
    if (any(alignments$clipped_right)) {
      right_clipped <- alignments[alignments$clipped_right, ]
      
      gg <- gg + ggplot2::geom_segment(
        data = right_clipped,
        ggplot2::aes(
          x = gend, xend = gend,
          y = rect_ymin, yend = rect_ymax,
          text = hover_text,
          key = read_id
        ),
        color = "black", size = 0.5
      )
    }
  }

  # baseline at y=0
  gg <- gg + ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)

  # Plot mutations
  if (!is.null(df$mutations) && nrow(df$mutations) > 0 && profile$full_style == "show_mutations") {
    cat(sprintf("plotting %d mutations\n", nrow(df$mutations)))
    mutations <- df$mutations
    mutations = mutations[is.element(mutations$read_id, df$reads$read_id), ]

    # sample if too many mutations
    if (nrow(mutations) > profile$max_mutations) {
      mutations <- mutations[sample(nrow(mutations), profile$max_mutations), ]
      cat(sprintf("sampled to %d mutations\n", profile$max_mutations))
    }

    mutations$hover_text <- paste0(
      "Read: ", mutations$read_id, "\n",
      "Position: ", mutations$coord, "\n",
      "Type: ", mutations$desc
    )

    # show rectangles when zoomed in under 100bp, otherwise segments
    xlim_range <- cxt$mapper$xlim[2] - cxt$mapper$xlim[1]
    if (xlim_range < 500) {
      # plot mutation rectangles
      gg <- gg + ggplot2::geom_rect(
        data = mutations,
        ggplot2::aes(
          xmin = gcoord - 0.4, xmax = gcoord + 0.4,
          ymin = height, ymax = height + 1,
          fill = desc,
          text = hover_text
        ),
        color = NA
      ) +
        ggplot2::scale_fill_discrete(name = "Mutation") +
        ggplot2::theme(legend.position = "top")
    } else {
      # plot mutation segments
      gg <- gg + ggplot2::geom_segment(
        data = mutations,
        ggplot2::aes(
          x = gcoord, xend = gcoord,
          y = height, yend = height + 1,
          color = desc,
          text = hover_text
        ),
        size = 0.25
      ) +
        ggplot2::scale_color_discrete(name = "Mutation") +
        ggplot2::theme(legend.position = "top")
    }
  }
  # !!!
  # cache_set("alns", df$reads)
  return(gg)
}
