get_intervals <- function(cdf, zoom) {
  if (!is.null(zoom)) {
    cdf <- cdf[cdf$start <= zoom[2] & cdf$end >= zoom[1], ]
  }
  data.frame(contig = cdf$contig, start = cdf$start + 1, end = cdf$end + 1)
}

build_context <- function(state_contigs, contig_table, zoom, assembly) {
  req(state_contigs)
  req(contig_table)
  valid_contigs <- state_contigs[state_contigs %in% contig_table$contig]
  if (length(valid_contigs) == 0) {
    return(NULL)
  }

  cdf <- contig_table[match(valid_contigs, contig_table$contig), ]

  if (nrow(cdf) == 0) {
    return(NULL)
  }
  cdf$start <- cumsum(c(0, head(cdf$length, -1)))
  cdf$end <- cdf$start + cdf$length
  mapper <- list(
    g2l = function(gcoords) {
      ii <- findInterval(gcoords, c(cdf$start, Inf), left.open = TRUE, all.inside = TRUE)
      if (any(ii <= 0 | ii > nrow(cdf))) stop("Some coordinate indices are out of range")
      data.frame(contig = cdf$contig[ii], coord = gcoords - cdf$start[ii], gcoord = gcoords)
    },
    l2g = function(contigs, coords) {
      idx <- match(contigs, cdf$contig)
      if (any(is.na(idx))) {
        stop(sprintf(
          "Some contigs are not found in the contig table: %s",
          paste(contigs[is.na(idx)], collapse = ", ")
        ))
      }
      cdf$start[idx] + coords
    },
    cdf = cdf
  )

  if (!is.null(zoom)) {
    mapper$xlim <- range(zoom)
    cdf <- cdf[cdf$start <= zoom[2] & cdf$end >= zoom[1], ]
  } else {
    mapper$xlim <- range(cdf$start, cdf$end)
  }
  intervals <- get_intervals(cdf, zoom)
  list(mapper = mapper, zoom = zoom, contigs = state_contigs, assembly = assembly, intervals = intervals)
}

# Safely filters point data and adds global coordinates
filter_coords <- function(df, cxt, xlim = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  # Check required columns
  if (!all(c("contig", "coord") %in% names(df))) {
    warning("missing required columns in filter_coords: contig, coord")
    return(NULL)
  }

  # Map contig names to those in context
  df <- df[df$contig %in% cxt$mapper$cdf$contig, ]
  if (nrow(df) == 0) {
    return(NULL)
  }

  # Add global coordinates
  df$gcoord <- cxt$mapper$l2g(df$contig, df$coord)

  # Filter by xlim if provided
  if (!is.null(xlim)) {
    in_range <- df$gcoord >= xlim[1] & df$gcoord <= xlim[2]
    if (!any(in_range)) {
      return(NULL)
    }
    df <- df[in_range, ]
  }

  return(df)
}

# Safely filters segment data and adds global coordinates
filter_segments <- function(df, cxt, xlim = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  # Check required columns
  if (!all(c("contig", "start", "end") %in% names(df))) {
    warning("missing required columns in filter_segments: contig, start, end")
    return(NULL)
  }

  # Map contig names to those in context
  df <- df[df$contig %in% cxt$mapper$cdf$contig, ]
  if (nrow(df) == 0) {
    return(NULL)
  }

  # Add global coordinates
  df$gstart <- cxt$mapper$l2g(df$contig, df$start)
  df$gend <- cxt$mapper$l2g(df$contig, df$end)

  # Filter by xlim if provided
  if (!is.null(xlim)) {
    # Keep segments with any overlap with xlim (partial or complete)
    has_overlap <- df$gstart <= xlim[2] & df$gend >= xlim[1]
    if (!any(has_overlap)) {
      return(NULL)
    }
    df <- df[has_overlap, ]
  }

  return(df)
}
