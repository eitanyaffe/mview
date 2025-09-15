segments_profile <- function(id, name, height = 30, is_fixed = TRUE,
                            segments_f = NULL,
                            color = "#2E86AB",
                            auto_register = TRUE) {

  plot_f <- function(profile, cxt, gg) {
    # get segments data
    segments <- NULL
    if (is.function(profile$segments_f)) {
      # function provided
      segments <- profile$segments_f(cxt$assembly)
    } else if (is.character(profile$segments_f)) {
      # cache key provided
      if (cache_exists(profile$segments_f)) {
        segments <- cache_get(profile$segments_f)
      }
    } else if (is.data.frame(profile$segments_f)) {
      # direct data frame provided
      segments <- profile$segments_f
    }
    
    if (is.null(segments) || nrow(segments) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # filter segments by current assembly first (if assembly column exists)
    if ("assembly" %in% colnames(segments)) {
      assembly_segments <- segments[segments$assembly == cxt$assembly, ]
      if (is.null(assembly_segments) || nrow(assembly_segments) == 0) {
        return(list(plot = gg, legends = list()))
      }
    } else {
      assembly_segments <- segments
    }
    
    # filter segments using context
    filtered_segments <- filter_segments(assembly_segments, cxt, cxt$mapper$xlim)
    if (is.null(filtered_segments) || nrow(filtered_segments) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    
    # use profile settings directly (minimal approach)
    current_color <- profile$color
    
    # create hover text
    hover_text <- paste0(
      filtered_segments$contig, ": ", 
      filtered_segments$start, "-", filtered_segments$end,
      "\n", filtered_segments$desc
    )
    
    # add segments layer as rectangles
    gg <- gg + 
      ggplot2::geom_rect(
        data = filtered_segments,
        ggplot2::aes(
          xmin = gstart, xmax = gend,
          ymin = -0.15, ymax = 0.15,
          text = hover_text
        ),
        fill = current_color,
        color = "black",
        size = 0.5
      ) +
      ggplot2::geom_text(
        data = filtered_segments,
        ggplot2::aes(
          x = gend + (cxt$mapper$xlim[2] - cxt$mapper$xlim[1]) * 0.017,
          y = -0.05,
          label = id,
          text = hover_text
        ),
        color = "black",
        size = 3,
        hjust = 0,
        vjust = 0.5
      ) +
      ggplot2::ylim(-0.5, 0.5)
    
    return(list(plot = gg, legends = list()))
  }

  # create profile
  profile_create(
    id = id, name = name, type = "segments", height = height, is_fixed = is_fixed,
    attr = list(hide_y_ticks = TRUE),
    params = NULL, plot_f = plot_f,
    segments_f = segments_f,
    color = color,
    auto_register = auto_register
  )
}
