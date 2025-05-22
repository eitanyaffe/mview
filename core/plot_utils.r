# Utility functions for plotting

#' Get interpolated colors from a palette based on values
#'
#' @param values Numeric vector of values to map to colors
#' @param colors Vector of colors to interpolate between
#' @param min_val Minimum value in the range (maps to first color)
#' @param max_val Maximum value in the range (maps to last color)
#' @param num_steps Number of discrete colors to generate (default: 20)
#' @return Vector of colors corresponding to each value
get_color_scale <- function(values, colors, min_val = 0, max_val = 1, num_steps = 20) {
  # Validate inputs
  if (!is.numeric(values)) {
    stop("values must be numeric")
  }
  if (length(colors) < 2) {
    stop("need at least 2 colors")
  }

  # Create color palette with specified number of steps
  color_palette <- colorRampPalette(colors)(num_steps)

  # Scale and constrain values to range [0,1]
  normalized <- (values - min_val) / (max_val - min_val)
  normalized <- pmin(pmax(normalized, 0), 1)

  # Handle NA or non-finite values
  normalized[is.na(normalized) | !is.finite(normalized)] <- 0

  # Convert to color indices (1 to num_steps)
  color_index <- ceiling(normalized * num_steps)
  color_index[color_index < 1] <- 1

  # Return colors from palette
  return(color_palette[color_index])
}
