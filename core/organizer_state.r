# organizer state management
# reactive values for segment list editing

# sequence_segments_rv stores a data frame with columns: segment, strand, removed
# this is a local copy for editing before Apply
# removed column: TRUE if segment is marked for removal (not applied until Apply is pressed)
sequence_segments_rv <- reactiveVal(data.frame(segment = character(), strand = character(), removed = logical(), stringsAsFactors = FALSE))

# get the current sequence segments as data frame
get_sequence_segments <- function() {
  return(sequence_segments_rv())
}

# set sequence segments
# accepts either character vector (will add default strand '+') or data frame with segment/strand
# if data frame doesn't have 'removed' column, adds it with FALSE values
set_sequence_segments <- function(value, caller = "unknown") {
  if (is.character(value)) {
    if (length(value) == 0) {
      value <- data.frame(segment = character(), strand = character(), removed = logical(), stringsAsFactors = FALSE)
    } else {
      value <- data.frame(segment = value, strand = rep("+", length(value)), removed = rep(FALSE, length(value)), stringsAsFactors = FALSE)
    }
  } else if (is.data.frame(value)) {
    if (!"removed" %in% names(value)) {
      if (nrow(value) == 0) {
        value$removed <- logical()
      } else {
        value$removed <- rep(FALSE, nrow(value))
      }
    }
  }
  sequence_segments_rv(value)
}

# get just the segment IDs as character vector
get_sequence_segment_ids <- function() {
  df <- sequence_segments_rv()
  if (is.null(df) || nrow(df) == 0) return(character())
  return(df$segment)
}

# flip strand for specified segments
flip_segment_strands <- function(segment_ids) {
  df <- sequence_segments_rv()
  if (is.null(df) || nrow(df) == 0) return()
  
  idx <- df$segment %in% segment_ids
  df$strand[idx] <- ifelse(df$strand[idx] == "+", "-", "+")
  sequence_segments_rv(df)
}

