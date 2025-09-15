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
    colors <- ifelse(alignments$along_contig, "#dcfbdc", "#efcb39")
    return(colors)
  } else if (style == "by_mutations") {
    # use shared mutation color function
    alignment_lengths <- alignments$end - alignments$start + 1
    mutations_per_bp <- (alignments$mutation_count / alignment_lengths)
    colors <- get_mutation_colors(mutations_per_bp)
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
clip_mode = "all",
clip_margin = 10,
min_mutations_percent = 0.0,
max_mutations_percent = 10.0,
min_alignment_length = 0,
max_alignment_length = 0,
min_indel_length = 3) {
  
  # Create cache key based on all relevant parameters
  # Use address of external pointer as unique identifier for alignment
  aln_id <- if (is(aln, "externalptr")) {
    paste0("ptr_", format(aln))
  } else {
    digest::digest(aln, algo = "md5")
  }
  
  cache_key <- paste0("full_query_",
                     digest::digest(list(
                       aln_id = aln_id,
                       xlim = cxt$mapper$xlim,
                       height_style_str = height_style_str,
                       max_reads = max_reads,
                       clip_mode = clip_mode,
                       clip_margin = clip_margin,
                       min_mutations_percent = min_mutations_percent,
                       max_mutations_percent = max_mutations_percent,
                       min_alignment_length = min_alignment_length,
                       max_alignment_length = max_alignment_length,
                       min_indel_length = min_indel_length
                     ), algo = "md5"))

  # Use cache for the full query
  df <- cache(cache_key, {
    aln_query_full(aln, cxt$intervals, height_style_str, max_reads, clip_mode, clip_margin, as.numeric(min_mutations_percent), as.numeric(max_mutations_percent), as.integer(min_alignment_length), as.integer(max_alignment_length), as.integer(min_indel_length))
  })
  cat(sprintf("done with aln_query_full\n"))
  
  # Process alignments
  alignments <- NULL
  if (!is.null(df$alignments) && nrow(df$alignments) > 0) {
    alns <- df$alignments
    alns$contig <- alns$contig_id
    # alntools full mode outputs 1-based coordinates  
    alns$start <- alns$contig_start
    alns$end <- alns$contig_end
    alignments <- filter_segments(alns, cxt, cxt$mapper$xlim)
  }

  # Process mutations
  mutations <- NULL
  if (!is.null(df$mutations) && nrow(df$mutations) > 0) {
    muts <- df$mutations
    muts$contig <- muts$contig_id
    # alntools full mode outputs 1-based coordinates
    muts$coord <- muts$position
    mutations <- filter_coords(muts, cxt, cxt$mapper$xlim)
  }

  # Process reads
  reads <- NULL
  if (!is.null(df$reads) && nrow(df$reads) > 0) {
    reads <- df$reads
    reads$contig <- reads$contig_id
    # alntools full mode outputs 1-based coordinates
    reads$start <- reads$span_start
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
    clip_mode = profile$clip_mode,
    clip_margin = profile$clip_margin,
    min_mutations_percent = as.numeric(profile$min_mutations_percent),
    max_mutations_percent = as.numeric(profile$max_mutations_percent),
    min_alignment_length = as.integer(profile$min_alignment_length),
    max_alignment_length = as.integer(profile$max_alignment_length),
    min_indel_length = as.integer(profile$min_indel_length))
  
  # Plot reads first (under alignments)
  if (!is.null(df$reads) && nrow(df$reads) > 0) {

    cat(sprintf("plotting %d reads\n", nrow(df$reads)))
    reads <- df$reads

    # create hover text only if enabled
    if (profile$show_hover) {
      reads$hover_text <- paste0(
        "Read: ", reads$read_id, "\n",
        reads$start, "-", reads$end
      )
    } else {
      reads$hover_text <- ""
    }

    # draw read lines at middle height
    gg <- gg + ggplot2::geom_segment(
      data = reads,
      ggplot2::aes(
        x = gstart, xend = gend,
        y = height + 0.5, yend = height + 0.5,
        text = hover_text,
        key = read_id
      ),
      color = "gray50", size = 0.3
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
    

    # set colors based on color mode
    alignments$color <- get_alignment_colors(alignments, reads, profile$full_style)

    alignments$align_length <- alignments$end - alignments$start + 1
    alignments$mut_density_pct <- (alignments$mutation_count / alignments$align_length) * 100
    
    # create hover text only if enabled
    if (profile$show_hover) {
      alignments$hover_text <- paste0(
        "Read: ", alignments$read_id, "\n",
        alignments$start, "-", alignments$end, "\n",
        "Read coords: ", alignments$read_start, "-", alignments$read_end, " / ", alignments$read_length, "\n",
        "Mutations: ", alignments$mutation_count, " (", sprintf("%.3f", alignments$mut_density_pct), "%)"
      )
    } else {
      alignments$hover_text <- ""
    }

    # determine if we should show gray borders based on view range
    range_bp <- (cxt$mapper$xlim[2] + 1) - cxt$mapper$xlim[1]
    show_borders <- range_bp <= 200000  # show borders only under 1Mb
    border_color <- if (show_borders) "gray50" else NA
    
    # draw alignment rectangles
    gg <- gg + ggplot2::geom_rect(
      data = alignments,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = rect_ymin, ymax = rect_ymax,
        text = hover_text,
        key = read_id
      ),
      fill = alignments$color, color = border_color, size = 0.2
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

  # Plot mutations using unified function
  if (!is.null(df$mutations) && nrow(df$mutations) > 0 && profile$full_style == "show_mutations") {
    mutations <- df$mutations
    mutations = mutations[is.element(mutations$read_id, df$reads$read_id), ]

    # build hover text for align-specific format
    if (profile$show_hover) {
      mutations$hover_text <- paste0(
        "Read: ", mutations$read_id, "\n",
        "Position: ", mutations$coord, "\n",
        "Type: ", mutations$desc
      )
    } else {
      mutations$hover_text <- ""
    }

    # add ybottom and ytop for align profile (height and height + 1)
    mutations$ybottom <- mutations$height
    mutations$ytop <- mutations$height + 1

    # use unified mutation plotting function
    gg <- plot_mutations_unified(gg, mutations, profile, cxt)
  }
  
  # apply force_max_y if set
  if (!is.null(profile$force_max_y) && profile$force_max_y > 0) {
    gg <- gg + ggplot2::coord_cartesian(ylim = c(NA, profile$force_max_y))
  }
  
  # !!!
  # cache_set("alns", df$reads)
  # build legends based on mutation_color_mode
  legends <- list()
  mutation_color_mode <- if (!is.null(profile$mutation_color_mode)) profile$mutation_color_mode else "detailed"
  if (mutation_color_mode == "detailed") {
    legend_gg <- create_detailed_variant_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 420, width = 350)))
    }
  } else if (mutation_color_mode == "type") {
    legend_gg <- create_simplified_variant_legend()
    if (!is.null(legend_gg)) {
      legends <- c(legends, list(list(gg = legend_gg, height = 140, width = 300)))
    }
  }
  
  return(list(plot = gg, legends = legends))
}
