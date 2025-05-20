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

format_bp <- function(bp) {
  rr <- function(x) round(x, 1)
  result <- ifelse(bp < 1000,
    paste0(bp, " bp"),
    ifelse(bp < 1000000,
      paste0(rr(bp / 1000), " kb"),
      ifelse(bp < 1000000000,
        paste0(rr(bp / 1000000), " mb"),
        paste0(rr(bp / 1000000000), " gb")
      )
    )
  )
  return(result)
}

print_state <- function(state, title = "State", show_title = TRUE) {
  # Format contig display based on count
  contig_display <- if (length(state$contigs) == 0) {
    "None"
  } else if (length(state$contigs) == 1) {
    as.character(state$contigs[1])
  } else {
    paste(length(state$contigs), "contigs")
  }

  assembly_display <- if (is.null(state$assembly) || state$assembly == "") {
    "None"
  } else {
    state$assembly
  }

  # Format zoom display based on available values
  zoom_display <- if (is.null(state$zoom) || length(state$zoom) != 2) {
    "Full range"
  } else {
    dd <- round(state$zoom[2]) - round(state$zoom[1])
    format_bp(dd)
  }
  basic_info <- paste0(
    "Subject: ", assembly_display, ", Contigs: ",
    contig_display, ", Zoom: ", zoom_display
  )
  if (show_title) {
    cat(title, ":", basic_info, "\n")
  } else {
    cat(basic_info, "\n")
  }
}
