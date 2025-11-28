# default parameters
default_axis_params <- list(
  show_nt = list(
    group_id = "axis",
    type = "boolean",
    default = TRUE
  ),
  nt_threshold = list(
    group_id = "axis",
    type = "integer",
    default = 200
  ),
  scale_factor = list(
    group_id = "axis",
    type = "double",
    default = 2
  )
)

# compute nice tick interval for a given range
compute_nice_interval <- function(visible_range, target_ticks = 10) {
  if (visible_range <= 0) return(NULL)
  raw_interval <- visible_range / target_ticks
  magnitude <- 10^floor(log10(raw_interval))
  normalized <- raw_interval / magnitude
  nice <- if (normalized <= 1) 1 else if (normalized <= 2) 2 else if (normalized <= 5) 5 else 10
  max(1, nice * magnitude)
}

# compute tick positions in local coordinates
compute_local_ticks <- function(local_start, local_end, nice_interval) {
  start_tick <- ceiling(local_start / nice_interval) * nice_interval
  end_tick <- floor(local_end / nice_interval) * nice_interval
  if (start_tick > end_tick) return(numeric(0))
  ticks <- seq(start_tick, end_tick, by = nice_interval)
  ticks <- unique(round(ticks))
  ticks[ticks != 0]
}

#' Create a simple axis profile
axis_profile <- function(id = "simple_axis",
                         name = "simple axis",
                         height = 40, is_fixed = TRUE,
                         params = default_axis_params,
                         auto_register = TRUE) {

  data_f <- function() NULL

  # extract sequence for a contig range from cached fasta
  get_sequence <- function(assembly, contig, start_pos, end_pos) {
    tryCatch({
      sequences <- get_fasta(assembly)
      if (is.null(sequences)) return(NULL)
      
      # find contig in fasta
      contig_seq <- NULL
      for (seq_name in names(sequences)) {
        if (seq_name == contig || grepl(paste0("^", contig, "($|\\s)"), seq_name)) {
          contig_seq <- sequences[[seq_name]]
          break
        }
      }
      if (is.null(contig_seq)) return(NULL)
      
      # extract range
      full_seq <- as.character(contig_seq)
      start_pos <- max(1, start_pos)
      end_pos <- min(nchar(full_seq), end_pos)
      if (start_pos > end_pos) return(NULL)
      substr(full_seq, start_pos, end_pos)
    }, error = function(e) NULL)
  }

  # compute tick annotations for current view (used by hover/click handlers)
  get_annotations_f <- function(profile) {
    contig_df <- cxt_get_entire_view()
    if (is.null(contig_df) || nrow(contig_df) == 0) return(NULL)

    zoom_xlim <- cxt_get_xlim()
    if (is.null(zoom_xlim)) zoom_xlim <- c(min(contig_df$vstart), max(contig_df$vend))
    
    # filter to visible contigs using virtual coords
    visible <- contig_df[contig_df$vstart < zoom_xlim[2] & contig_df$vend > zoom_xlim[1], , drop = FALSE]
    if (nrow(visible) != 1) return(NULL)

    # compute local visible range
    vstart <- visible$vstart[1]
    local_start <- visible$start[1]
    local_end <- visible$end[1]
    local_visible_start <- max(local_start, zoom_xlim[1] - vstart + local_start)
    local_visible_end <- min(local_end, zoom_xlim[2] - vstart + local_start)
    
    nice_interval <- compute_nice_interval(local_visible_end - local_visible_start)
    if (is.null(nice_interval)) return(NULL)
    
    local_ticks <- compute_local_ticks(local_visible_start, local_visible_end, nice_interval)
    if (length(local_ticks) == 0) return(NULL)

    # convert to virtual coords for plotting
    global_ticks <- vstart + local_ticks - local_start
    data.frame(
      x = global_ticks,
      y = 0,
      label = format(local_ticks, big.mark = ",", scientific = FALSE)
    )
  }

  plot_f <- function(profile, gg) {
    contig_df <- cxt_get_entire_view()
    if (is.null(contig_df) || nrow(contig_df) == 0) return(list(plot = gg, legends = list()))

    # get zoom range in virtual coords
    xlim <- cxt_get_xlim()
    zoom_xlim <- if (!is.null(xlim) && length(xlim) == 2) {
      range(xlim)
    } else {
      c(min(contig_df$vstart), max(contig_df$vend))
    }
    
    # filter to visible contigs using virtual coords
    visible_df <- contig_df[contig_df$vstart < zoom_xlim[2] & contig_df$vend > zoom_xlim[1], , drop = FALSE]
    if (nrow(visible_df) == 0) return(list(plot = gg, legends = list()))

    window_size <- zoom_xlim[2] - zoom_xlim[1] + 1

    # add contig name labels (clip to zoom range so labels stay visible)
    if (nrow(visible_df) == 1) {
      name_data <- data.frame(x = (zoom_xlim[1] + zoom_xlim[2]) / 2, y = 0.45, label = visible_df$contig)
    } else {
      clipped_start <- pmax(visible_df$vstart, zoom_xlim[1])
      clipped_end <- pmin(visible_df$vend, zoom_xlim[2])
      name_data <- data.frame(x = (clipped_start + clipped_end) / 2, y = 0.55, label = visible_df$contig)
    }
    gg <- gg + ggplot2::geom_text(data = name_data, ggplot2::aes(x = x, y = y, label = label), 
                                   color = "black", size = 3, vjust = 0.5)

    # single contig: draw axis line and ticks
    if (nrow(visible_df) == 1) {
      vstart <- visible_df$vstart[1]
      vend <- visible_df$vend[1]
      local_start <- visible_df$start[1]
      local_end <- visible_df$end[1]
      
      # axis line clipped to zoom range
      axis_start <- max(vstart, zoom_xlim[1])
      axis_end <- min(vend, zoom_xlim[2])
      axis_data <- data.frame(x = axis_start, y = 0.8, xend = axis_end, yend = 0.8)
      gg <- gg + ggplot2::geom_segment(data = axis_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), 
                                        color = "black", size = 0.4)

      # compute visible range in local coords
      local_visible_start <- max(local_start, zoom_xlim[1] - vstart + local_start)
      local_visible_end <- min(local_end, zoom_xlim[2] - vstart + local_start)
      local_range <- local_visible_end - local_visible_start
      
      if (local_range > 0) {
        # add tick marks
        nice_interval <- compute_nice_interval(local_range)
        if (!is.null(nice_interval)) {
          local_ticks <- compute_local_ticks(local_visible_start, local_visible_end, nice_interval)
          if (length(local_ticks) > 0) {
            global_ticks <- vstart + local_ticks - local_start
            tick_data <- data.frame(x = global_ticks, y = 0.8, xend = global_ticks, yend = 0.75)
            gg <- gg + ggplot2::geom_segment(data = tick_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), 
                                              color = "black", size = 0.3)
          }
        }
        
        # show nucleotides if zoomed in enough
        show_nt <- get_param("axis", "show_nt")()
        nt_threshold <- get_param("axis", "nt_threshold")()
        if (show_nt && window_size <= nt_threshold) {
          contig <- visible_df$contig[1]
          assembly <- cxt_get_assembly()
          sequence <- get_sequence(assembly, contig, local_visible_start, local_visible_end)
          
          if (!is.null(sequence) && nchar(sequence) > 0) {
            nt_local_positions <- seq(local_visible_start, local_visible_start + nchar(sequence) - 1)
            nt_global_positions <- vstart + nt_local_positions - local_start
            nt_chars <- strsplit(sequence, "")[[1]]
            
            # font size scales with zoom
            scale_factor <- get_param("axis", "scale_factor")()
            font_size <- scale_factor * max(0.1, min(2.0, 100.0 / window_size))
            
            nt_data <- data.frame(x = nt_global_positions, y = 0.85, label = toupper(nt_chars))
            gg <- gg + ggplot2::geom_text(data = nt_data, ggplot2::aes(x = x, y = y, label = label), 
                                           color = "darkblue", size = font_size, family = "mono", vjust = 0.5)
          }
        }
      }
    }

    list(plot = gg + ggplot2::coord_cartesian(ylim = c(0.35, 0.9), clip = "off"), legends = list())
  }

  profile_create(
    id = id,
    name = name,
    type = "axis",
    height = height, is_fixed = is_fixed,
    attr = list(hide_y_label = TRUE, hide_y_ticks = TRUE),
    params = params,
    data_f = data_f,
    plot_f = plot_f,
    get_annotations = get_annotations_f,
    auto_register = auto_register
  )
}