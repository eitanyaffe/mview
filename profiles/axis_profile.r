  
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
  
#' Create a simple axis profile
#' @param id Unique profile id
#' @param name Display name
#' @param height Height in pixels
#' @param auto_register Whether to register immediately
#' @return Profile object
axis_profile <- function(id = "simple_axis",
                         name = "simple axis",
                         height = 100,
                         params = default_axis_params,
                         auto_register = TRUE) {

  data_f <- function(cxt) NULL

  # helper function to extract sequence for a contig and range from cached fasta data
  get_sequence <- function(assembly, contig, start_pos, end_pos) {
    tryCatch({
      # get cached fasta sequences
      sequences <- get_fasta(assembly)
      
      if (is.null(sequences)) {
        return(NULL)
      }
      
      # find the target contig
      contig_seq <- NULL
      for (seq_name in names(sequences)) {
        if (seq_name == contig || grepl(paste0("^", contig, "($|\\s)"), seq_name)) {
          contig_seq <- sequences[[seq_name]]
          break
        }
      }
      
      if (is.null(contig_seq)) {
        return(NULL)
      }
      
      # get the full sequence string
      full_sequence <- as.character(contig_seq)
      
      # extract the requested range (1-based coordinates)
      if (start_pos < 1) start_pos <- 1
      if (end_pos > nchar(full_sequence)) end_pos <- nchar(full_sequence)
      if (start_pos > end_pos) return(NULL)
      
      substr(full_sequence, start_pos, end_pos)
    }, error = function(e) {
      cat(sprintf("error reading sequence: %s\n", e$message))
      NULL
    })
  }

  # compute coordinate annotations for the current view
  get_annotations_f <- function(profile, cxt) {
    if (is.null(cxt$mapper) || is.null(cxt$mapper$cdf) || nrow(cxt$mapper$cdf) == 0) return(NULL)

    contig_df_all <- cxt$mapper$cdf
    zoom_xlim <- if (!is.null(cxt$mapper$xlim)) range(cxt$mapper$xlim) else c(min(contig_df_all$start), max(contig_df_all$end))
    contig_df <- contig_df_all[contig_df_all$start < zoom_xlim[2] & contig_df_all$end > zoom_xlim[1], , drop = FALSE]

    if (nrow(contig_df) != 1) return(NULL)

    contig_start <- contig_df$start[1]
    contig_end <- contig_df$end[1]

    # compute tick marks in contig-local coordinates
    local_visible_start <- max(0, zoom_xlim[1] - contig_start)
    local_visible_end <- min(contig_end - contig_start, zoom_xlim[2] - contig_start)
    local_visible_range <- local_visible_end - local_visible_start
    if (local_visible_range <= 0) return(NULL)

    target_ticks <- 10
    raw_interval <- local_visible_range / target_ticks
    magnitude <- 10^floor(log10(raw_interval))
    normalized <- raw_interval / magnitude
    nice_interval <- if (normalized <= 1) {
        1 * magnitude
    } else if (normalized <= 2) {
        2 * magnitude
    } else if (normalized <= 5) {
        5 * magnitude
    } else {
        10 * magnitude
    }
    nice_interval <- max(1, nice_interval)

    start_tick_local <- ceiling(local_visible_start / nice_interval) * nice_interval
    end_tick_local <- floor(local_visible_end / nice_interval) * nice_interval
    local_ticks <- if (start_tick_local <= end_tick_local) {
        seq(start_tick_local, end_tick_local, by = nice_interval)
    } else {
        numeric(0)
    }
    if (length(local_ticks) > 0) {
        local_ticks <- round(local_ticks)
        local_ticks <- unique(local_ticks)
        local_ticks <- local_ticks[local_ticks != 0]
    }
    if (length(local_ticks) == 0) return(NULL)

    global_ticks <- contig_start + local_ticks
    cat("axis_profile: returning ", length(global_ticks), " tick positions\n")
    return(data.frame(
        x = global_ticks,
        y = 0,
        label = format(local_ticks, big.mark = ",", scientific = FALSE)
    ))
  }

  plot_f <- function(profile, cxt, gg) {
    if (is.null(cxt$mapper) || is.null(cxt$mapper$cdf) || nrow(cxt$mapper$cdf) == 0) return(list(plot = gg, legends = list()))

    contig_df_all <- cxt$mapper$cdf
    zoom_xlim <- if (!is.null(cxt$mapper$xlim)) range(cxt$mapper$xlim) else c(min(contig_df_all$start), max(contig_df_all$end))
    contig_df <- contig_df_all[contig_df_all$start < zoom_xlim[2] & contig_df_all$end > zoom_xlim[1], , drop = FALSE]

    if (nrow(contig_df) == 0) return(list(plot = gg, legends = list()))
    
    # get threshold parameter
    nt_threshold <- get_param("axis", "nt_threshold")()
    window_size <- zoom_xlim[2] - zoom_xlim[1] + 1

    if (nrow(contig_df) == 1) {
      name_data <- data.frame(x = (zoom_xlim[1] + zoom_xlim[2]) / 2, y = 0.45, label = contig_df$contig)
    } else {
      name_data <- data.frame(x = (contig_df$start + contig_df$end) / 2, y = 0.55, label = contig_df$contig)
    }
    gg <- gg + ggplot2::geom_text(data = name_data, ggplot2::aes(x = x, y = y, label = label), color = "black", size = 3, vjust = 0.5)

    if (nrow(contig_df) == 1) {
      axis_start <- max(min(contig_df$start), zoom_xlim[1])
      axis_end <- min(max(contig_df$end), zoom_xlim[2])
      axis_line_data <- data.frame(x = axis_start, y = 0.8, xend = axis_end, yend = 0.8)
      gg <- gg + ggplot2::geom_segment(data = axis_line_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 0.4)
    }

    if (nrow(contig_df) == 1) {
      contig_start <- contig_df$start[1]
      contig_end <- contig_df$end[1]

      # compute tick marks in contig-local coordinates
      local_visible_start <- floor(max(0, zoom_xlim[1] - contig_start))
      local_visible_end <- ceiling(min(contig_end - contig_start, zoom_xlim[2] - contig_start))
      local_visible_range <- local_visible_end - local_visible_start
      if (local_visible_range > 0) {
        target_ticks <- 10
        raw_interval <- local_visible_range / target_ticks
        magnitude <- 10^floor(log10(raw_interval))
        normalized <- raw_interval / magnitude
        nice_interval <- if (normalized <= 1) 1 * magnitude else if (normalized <= 2) 2 * magnitude else if (normalized <= 5) 5 * magnitude else 10 * magnitude
        nice_interval <- max(1, nice_interval)

        start_tick_local <- ceiling(local_visible_start / nice_interval) * nice_interval
        end_tick_local <- floor(local_visible_end / nice_interval) * nice_interval
        local_ticks <- if (start_tick_local <= end_tick_local) seq(start_tick_local, end_tick_local, by = nice_interval) else numeric(0)
        if (length(local_ticks) > 0) {
          local_ticks <- round(local_ticks)
          local_ticks <- unique(local_ticks)
          local_ticks <- local_ticks[local_ticks != 0]
        }
        if (length(local_ticks) > 0) {
          tick_height <- 0.05
          global_ticks <- contig_start + local_ticks
          tick_data <- data.frame(x = global_ticks, y = 0.8, xend = global_ticks, yend = 0.8 - tick_height)
          gg <- gg + ggplot2::geom_segment(data = tick_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 0.3)
        }
        
        # show nucleotides if enabled and window is small enough
        show_nt <- get_param("axis", "show_nt")()
        if (show_nt && window_size <= nt_threshold) {
          contig <- contig_df$contig[1]
          assembly <- cxt$assembly
          cat("axis_profile: showing nucleotides for contig ", contig, " in assembly ", assembly, "\n")
          sequence <- get_sequence(assembly, contig, 
                                 local_visible_start, local_visible_end)
          
          if (!is.null(sequence) && nchar(sequence) > 0) {
            # create nucleotide positions and labels
            nt_positions <- seq(local_visible_start + 1, local_visible_end)
            nt_chars <- strsplit(sequence, "")[[1]]
            
            # limit to visible range
            visible_indices <- nt_positions >= (zoom_xlim[1] - contig_start + 1) & 
                             nt_positions <= (zoom_xlim[2] - contig_start + 1)
            nt_positions <- nt_positions[visible_indices]
            nt_chars <- nt_chars[visible_indices]
            
            if (length(nt_positions) > 0) {
              global_nt_positions <- contig_start + nt_positions - 1
              
              # calculate font size based on zoom level to avoid collisions
              scale_factor <- get_param("axis", "scale_factor")()
              size_factor <- max(0.1, min(2.0, 100.0 / window_size))
              font_size <- scale_factor * size_factor
              
              # debug output
              cat(sprintf("axis nucleotides: window_size=%.1f, nt_threshold=%d, scale_factor=%.2f, size_factor=%.2f, font_size=%.2f\n", 
                         window_size, nt_threshold, scale_factor, size_factor, font_size))
              
              nt_data <- data.frame(x = global_nt_positions, y = 0.85, label = toupper(nt_chars))
              gg <- gg + ggplot2::geom_text(data = nt_data, 
                                          ggplot2::aes(x = x, y = y, label = label), 
                                          color = "darkblue", size = font_size, 
                                          family = "mono", vjust = 0.5)
            }
          }
        }
      }
    }

    return(list(plot = gg + ggplot2::coord_cartesian(ylim = c(0.35, 0.9), clip = "off"), legends = list()))
  }

  profile_create(
    id = id,
    name = name,
    type = "axis",
    height = height,
    attr = list(hide_y_label = TRUE, hide_y_ticks = TRUE),
    params = params,
    data_f = data_f,
    plot_f = plot_f,
    get_annotations = get_annotations_f,
    auto_register = auto_register
  )
}