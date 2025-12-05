# Full mode for alignment profile

# Get colors for alignments based on color mode
get_alignment_colors <- function(alignments, chunks, style) {
  if (style == "none") {
    return(rep("lightgray", nrow(alignments)))
  } else if (style == "by_strand") {
    # strand_flipped is already computed in C++: true if same contig as first alignment
    # in chunk AND different strand. Just use it directly.
    colors <- ifelse(alignments$strand_flipped, "#efcb39", "#dcfbdc")
    
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

align_query_full_mode <- function(aln, 
height_style_str = "by_mutations", 
max_reads = 100000,
clip_mode = "all",
clip_margin = 10,
min_mutations_percent = 0.0,
max_mutations_percent = 10.0,
min_alignment_length = 0,
max_alignment_length = 0,
max_margin = 10,
chunk_type = "break_on_overlap",
min_indel_length = 3) {
  
  intervals <- cxt_get_zoom_view()
  plotted_segs <- cxt_get_plotted_segments()
  xlim <- cxt_get_xlim()
  
  # Create cache key based on all relevant parameters
  # Use address of external pointer as unique identifier for alignment
  aln_id <- if (is(aln, "externalptr")) {
    paste0("ptr_", format(aln))
  } else {
    digest::digest(aln, algo = "md5")
  }
  
  # include plotted segments with strands for cache invalidation when segments flip
  seg_key <- if (nrow(plotted_segs) > 0) {
    paste(plotted_segs$segment, plotted_segs$strand, collapse = ",")
  } else {
    ""
  }
  
  cache_key <- paste0("full_query_",
                     digest::digest(list(
                       aln_id = aln_id,
                       seg_key = seg_key,
                       xlim = xlim,
                       height_style_str = height_style_str,
                       max_reads = max_reads,
                       clip_mode = clip_mode,
                       clip_margin = clip_margin,
                       min_mutations_percent = min_mutations_percent,
                       max_mutations_percent = max_mutations_percent,
                       min_alignment_length = min_alignment_length,
                       max_alignment_length = max_alignment_length,
                       max_margin = max_margin,
                       chunk_type = chunk_type,
                       min_indel_length = min_indel_length
                     ), algo = "md5"))
  
  # Use cache for the full query
  df <- cache(cache_key, {
    aln_query_full(aln, intervals, height_style_str, max_reads, clip_mode, 
    clip_margin, as.numeric(min_mutations_percent), as.numeric(max_mutations_percent), 
    as.integer(min_alignment_length), as.integer(max_alignment_length), as.integer(max_margin), 
    chunk_type, as.integer(min_indel_length))
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
    # use vstart/vend from C++ (already clipped to interval)
    alns$gstart <- alns$vstart
    alns$gend <- alns$vend
    alignments <- alns
  }

  # Process mutations
  mutations <- NULL
  if (!is.null(df$mutations) && nrow(df$mutations) > 0) {
    muts <- df$mutations
    muts$contig <- muts$contig_id
    # alntools full mode outputs 1-based coordinates
    muts$coord <- muts$position
    mutations <- cxt_filter_coords(muts)
  }

  # Process chunks - they already have vcoords from alntools
  chunks <- NULL
  if (!is.null(df$chunks) && nrow(df$chunks) > 0) {
    chunks <- df$chunks
    # chunks now have vstart/vend which are view coordinates
    chunks$gstart <- chunks$vstart
    chunks$gend <- chunks$vend
  }
  
  return(list(alignments = alignments, mutations = mutations, chunks = chunks))
}

align_profile_full <- function(profile, aln, gg) {
  intervals <- cxt_get_zoom_view()
  height_style <- profile$height_style
  if (!is.element(height_style, c("by_mutations", "by_coord_left", "by_coord_right"))) {
    height_style <- "by_mutations"
  }

  df <- align_query_full_mode(
    aln = aln, 
    height_style = height_style, 
    max_reads = profile$max_reads, 
    clip_mode = profile$clip_mode,
    clip_margin = profile$clip_margin,
    min_mutations_percent = as.numeric(profile$min_mutations_percent),
    max_mutations_percent = as.numeric(profile$max_mutations_percent),
    min_alignment_length = as.integer(profile$min_alignment_length),
    max_alignment_length = as.integer(profile$max_alignment_length),
    max_margin = as.integer(profile$max_margin),
    chunk_type = profile$chunk_type,
    min_indel_length = as.integer(profile$min_indel_length))
  
  # Plot chunks first (under alignments)
  if (!is.null(df$chunks) && nrow(df$chunks) > 0) {

    cat(sprintf("plotting %d chunks\n", nrow(df$chunks)))
    chunks <- df$chunks

    # create composite key with profile id for alignment tab lookup
    chunks$click_key <- paste0(profile$id, ":", chunks$read_id)

    # create hover text only if enabled
    if (profile$show_hover) {
      chunks$hover_text <- paste0(
        "Chunk: ", chunks$chunk_id, "\n",
        "Read: ", chunks$read_id, "\n",
        "Alignments: ", chunks$num_alignments
      )
    } else {
      chunks$hover_text <- ""
    }

    # draw chunk lines at middle height
    gg <- gg + ggplot2::geom_segment(
      data = chunks,
      ggplot2::aes(
        x = gstart, xend = gend,
        y = height + 0.5, yend = height + 0.5,
        text = hover_text,
        key = click_key
      ),
      color = "gray50", size = 0.3
    )
  }
  
  # Plot alignments
  chunks <- df$chunks
  if (!is.null(df$alignments) && nrow(df$alignments) > 0) {
    cat(sprintf("plotting %d alignments\n", nrow(df$alignments)))

    alignments <- df$alignments
    # filter to alignments with chunks if chunks exist
    if (!is.null(chunks) && nrow(chunks) > 0) {
      alignments <- alignments[is.element(alignments$chunk_id, chunks$chunk_id), ]
    }
    
    # exit early if no alignments remain after filtering
    if (is.null(alignments) || nrow(alignments) == 0) {
      gg <- gg + ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
      return(list(plot = gg, legends = list()))
    }
    
    alignments$rect_ymin <- alignments$height
    alignments$rect_ymax <- alignments$height + 1

    # create composite key with profile id for alignment tab lookup
    alignments$click_key <- paste0(profile$id, ":", alignments$read_id)

    # determine clipped status
    xlim <- cxt_get_xlim()
    left_in_plot <- alignments$gstart > xlim[1]
    right_in_plot <- alignments$gend < xlim[2]

    # check if read is soft-clipped (doesn't start at 0 or doesn't end at read_length)
    clip_read_start <- alignments$read_start > profile$max_margin
    clip_read_end <- alignments$read_end < alignments$read_length - profile$max_margin

    # account for segment strand when determining clipping side in view coordinates
    # when segment is minus strand, gstart corresponds to higher contig coord (right side)
    plotted_segs <- cxt_get_plotted_segments()
    seg_strand_lookup <- setNames(plotted_segs$strand, plotted_segs$contig)
    aln_seg_strand <- seg_strand_lookup[alignments$contig_id]
    aln_seg_strand[is.na(aln_seg_strand)] <- "+"
    seg_is_minus <- aln_seg_strand == "-"
    
    # effective reverse in view coords: XOR of alignment reverse and segment minus
    effective_reverse <- xor(alignments$is_reverse, seg_is_minus)
    
    # determine if read is actually clipped (soft-clipped or interval-clipped)
    alignments$clipped_left <- ifelse(!effective_reverse, clip_read_start, clip_read_end)
    alignments$clipped_right <- ifelse(!effective_reverse, clip_read_end, clip_read_start)
    
    # only show clipping indicators if read is actually clipped AND edge is visible in view
    # black vertical lines indicate alignment ends in middle of read, not just outside view
    # alignments$clipped_left <- (read_clip_left | interval_clip_left) & left_in_plot
    # alignments$clipped_right <- (read_clip_right | interval_clip_right) & right_in_plot

    # set colors based on color mode
    alignments$color <- get_alignment_colors(alignments, chunks, profile$full_style)

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
    range_bp <- (xlim[2] + 1) - xlim[1]
    show_borders <- range_bp <= 200000  # show borders only under 1Mb
    border_color <- if (show_borders) "gray50" else NA
    
    # draw alignment rectangles
    gg <- gg + ggplot2::geom_rect(
      data = alignments,
      ggplot2::aes(
        xmin = gstart, xmax = gend,
        ymin = rect_ymin, ymax = rect_ymax,
        text = hover_text,
        key = click_key
      ),
      fill = alignments$color, color = border_color, size = 0.2
    )
    
    # draw red rectangles for selected read
    selected_read_id <- cache_get_if_exists("selected_read_id", NULL)
    if (!is.null(selected_read_id) && length(selected_read_id) > 0) {
      selected_alignments <- alignments[alignments$read_id == selected_read_id, ]
      if (nrow(selected_alignments) > 0) {
        gg <- gg + ggplot2::geom_rect(
          data = selected_alignments,
          ggplot2::aes(
            xmin = gstart, xmax = gend,
            ymin = rect_ymin, ymax = rect_ymax,
            text = hover_text,
            key = click_key
          ),
          fill = "red", color = NA
        )
      }
    }

    # draw left clipping indicators
    if (any(alignments$clipped_left)) {
      left_clipped <- alignments[alignments$clipped_left, ]
      gg <- gg + ggplot2::geom_segment(
        data = left_clipped,
        ggplot2::aes(
          x = gstart, xend = gstart,
          y = rect_ymin, yend = rect_ymax,
          text = hover_text,
          key = click_key
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
          key = click_key
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
    # filter mutations to those belonging to alignments in visible chunks
    visible_chunk_alignments <- alignments$alignment_index
    mutations = mutations[is.element(mutations$alignment_index, visible_chunk_alignments), ]

    # only proceed if mutations remain after filtering
    if (!is.null(mutations) && nrow(mutations) > 0) {
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
      gg <- plot_mutations_unified(gg, mutations, profile)
    }
  }
  
  # apply force_max_y if set
  if (!is.null(profile$force_max_y) && profile$force_max_y > 0) {
    gg <- gg + ggplot2::coord_cartesian(ylim = c(NA, profile$force_max_y))
  }
  
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
