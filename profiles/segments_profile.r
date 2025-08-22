# default parameters for segments profile
default_segments_params <- list(
  height = list(
    group_id = "segments",
    type = "integer",
    default = 40
  )
)

segments_profile <- function(id, name, height = 40,
                            segments_data = NULL,
                            color = "#2E86AB",
                            params = default_segments_params,
                            auto_register = TRUE) {

  plot_f <- function(profile, cxt, gg) {
    # get segments data
    segments <- NULL
    if (is.character(profile$segments_data)) {
      # cache key provided
      if (cache_exists(profile$segments_data)) {
        segments <- cache_get(profile$segments_data)
      }
    } else if (is.data.frame(profile$segments_data)) {
      # direct data frame provided
      segments <- profile$segments_data
    }
    
    if (is.null(segments) || nrow(segments) == 0) {
      return(list(plot = gg, legends = list()))
    }
    
    # filter segments by current assembly first
    assembly_segments <- segments[segments$assembly == cxt$assembly, ]
    if (is.null(assembly_segments) || nrow(assembly_segments) == 0) {
      return(list(plot = gg, legends = list()))
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
          ymin = -0.15, ymax = 0.15
        ),
        fill = current_color,
        color = "black",
        size = 0.5
      ) +
      ggplot2::geom_text(
        data = filtered_segments,
        ggplot2::aes(
          x = gend + (cxt$mapper$xlim[2] - cxt$mapper$xlim[1]) * 0.01,
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
    id = id, name = name, type = "segments", height = height,
    attr = list(hide_y_ticks = TRUE),
    params = params, plot_f = plot_f,
    segments_data = segments_data,
    color = color,
    auto_register = auto_register
  )
}
