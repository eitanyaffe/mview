# profiles/points_profile.r

# Helper for creating a points-based profile
# Depends on profile_create from core/profile_manager.r

#' Create a points profile
#' @param id Unique profile id
#' @param name Display name
#' @param height Relative height weight
#' @param data Data.frame or function(cxt) returning data.frame with 'contig' and 'coord' columns
#' @param y_col Name of y column
#' @param desc_col Name of description column (optional)
#' @param point_aes List of static aesthetics for geom_point
#' @param attr List of additional attributes
#' @param auto_register Whether to register immediately
#' @return Profile object
points_profile <- function(id, name, height = 1,
                           data,
                           y_col = "value", desc_col = "desc",
                           color = "black", size = 1,
                           point_aes = list(),
                           params = list(),
                           auto_register = TRUE) {
  # data_f: get data, map to gcoord using mapper, filter by xlim if present
  data_f <- function(cxt) {
    df <- if (is.function(data)) data(cxt) else data
    if (!all(c("contig", "coord", y_col) %in% names(df))) stop("Missing required columns in points data.")

    # Use the new filter_coords function to handle coordinate mapping and filtering
    df <- filter_coords(df, cxt, cxt$mapper$xlim)

    cat(sprintf("data_f: %s, nrow: %d\n", name, nrow(df)))
    # Return empty data frame if filtering resulted in no rows
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame())
    }

    return(df)
  }
  # plot_f: add points layer, return info_df for hover
  plot_f <- function(profile, cxt, gg) {
    df <- data_f(cxt)
    if (nrow(df) == 0) {
      # Return just the plot object, no info_df
      return(gg)
    }

    # Prepare description for hover text
    if (!is.null(desc_col) && desc_col %in% names(df)) {
      df$hover_text <- df[[desc_col]]
    } else {
      df$hover_text <- paste0(name, ": ", round(df[[y_col]], 2))
    }

    aes_map <- ggplot2::aes(x = gcoord, y = .data[[y_col]], text = hover_text)
    geom_args <- c(list(mapping = aes_map, data = df), point_aes)

    geom_args$color <- verify_color(profile$color %||% color)

    # Suppress warnings specifically for the geom_point call regarding unknown aesthetics like 'text'
    gg2 <- gg + suppressWarnings(do.call(ggplot2::geom_point, geom_args))

    # Return just the plot object
    return(gg2)
  }

  profile_create(
    id = id,
    name = name,
    type = "points",
    height = height,
    params = params,
    data_f = data_f,
    plot_f = plot_f,
    auto_register = auto_register
  )
}
