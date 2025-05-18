# verify this a valid color, other give black
verify_color <- function(color_val, default_color = "black") {
  # Validate color value and return it if valid, otherwise return "black"
  if (is.null(color_val) || !is.character(color_val) ||
    nchar(color_val) == 0 ||
    !is.character(default_color) ||
    nchar(default_color) == 0) {
    return(default_color)
  }

  # Try to convert to RGB to check if it's a valid color
  valid <- tryCatch(
    {
      grDevices::col2rgb(color_val)
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  return(if (valid) color_val else default_color)
}

print_state <- function(state, title = "State") {
  # Format contig display based on count
  contig_display <- if (length(state$contigs) == 0) {
    "None"
  } else if (length(state$contigs) == 1) {
    as.character(state$contigs[1])
  } else {
    paste(length(state$contigs), "contigs")
  }

  # Format zoom display based on available values
  zoom_display <- if (is.null(state$zoom) || length(state$zoom) != 2) {
    "Full range"
  } else {
    paste(state$zoom[1], "–", state$zoom[2])
  }
  cat(title, ":", contig_display, zoom_display, "\n")
}
