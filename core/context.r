build_context <- function(state_contigs, contig_table, zoom) {
    req(state_contigs)
    req(contig_table)
  valid_cids <- intersect(state_contigs, contig_table$cid)
  cdf <- contig_table[match(valid_cids, contig_table$cid), ]
  if (nrow(cdf) == 0) {
    return(NULL)
  }
  cdf$start <- cumsum(c(0, head(cdf$length, -1)))
  cdf$end <- cdf$start + cdf$length
  mapper <- list(
    g2l = function(gcoords) {
      ii <- findInterval(gcoords, c(cdf$start, Inf), left.open = TRUE, all.inside = TRUE)
      if (any(ii <= 0 | ii > nrow(cdf))) stop("Some coordinate indices are out of range")
      data.frame(cid = cdf$cid[ii], coord = gcoords - cdf$start[ii], gcoord = gcoords)
    },
    l2g = function(df) {
      idx <- match(df$cid, cdf$cid)
      if (any(is.na(idx))) {
        stop(sprintf(
          "Some cids are not found in the contig table: %s",
          paste(df$cid[is.na(idx)], collapse = ", ")
        ))
      }
      cdf$start[idx] + df$coord
    },
    cdf = cdf
  )
  if (!is.null(zoom)) {
    mapper$xlim <- range(zoom)
  } else {
    mapper$xlim <- range(cdf$start, cdf$end)
  }
  list(mapper = mapper, zoom = zoom, contigs=state_contigs)
}
