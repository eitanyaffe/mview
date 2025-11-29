# selector tab state management
# reactive values and debug wrappers for sequence_segments

# reactive values
selected_graph_segment <- reactiveVal(NULL)
selected_graph_edge <- reactiveVal(NULL)
selected_sequence_index <- reactiveVal(NULL)
graph_segments <- reactiveVal(character())
sequence_segments_rv <- reactiveVal(character())
min_support <- reactiveVal(cache_get_if_exists("selector.min_support", 0))
min_percent <- reactiveVal(cache_get_if_exists("selector.min_percent", 0))

# wrapper functions for sequence_segments with debug logging
get_sequence_segments <- function() {
  val <- sequence_segments_rv()
  cat(sprintf("[GET] sequence_segments: length=%d, has_NA=%s\n", 
              length(val), any(is.na(val))))
  if (any(is.na(val))) {
    cat(sprintf("[GET] NA values at indices: %s\n", paste(which(is.na(val)), collapse=", ")))
  }
  return(val)
}

set_sequence_segments <- function(value, caller = "unknown") {
  cat(sprintf("[SET:%s] sequence_segments: length=%d, has_NA=%s\n", 
              caller, length(value), any(is.na(value))))
  if (any(is.na(value))) {
    cat(sprintf("[SET:%s] WARNING: setting NA values at indices: %s\n", 
                caller, paste(which(is.na(value)), collapse=", ")))
    cat(sprintf("[SET:%s] values: %s\n", caller, paste(head(value, 10), collapse=", ")))
    cat(sprintf("[SET:%s] traceback:\n", caller))
    traceback()
  }
  sequence_segments_rv(value)
}

