# selector tab state management
# reactive values and debug wrappers for sequence_segments

# reactive values
selected_graph_segments <- reactiveVal(character())
selected_graph_edge <- reactiveVal(NULL)
selected_sequence_index <- reactiveVal(NULL)
hovered_node <- reactiveVal(NULL)
hovered_edge <- reactiveVal(NULL)
sequence_segments_rv <- reactiveVal(character())
min_support <- reactiveVal(cache_get_if_exists("selector.min_support", 0))
min_percent <- reactiveVal(cache_get_if_exists("selector.min_percent", 0))

# wrapper functions for sequence_segments
get_sequence_segments <- function() {
  return(sequence_segments_rv())
}

set_sequence_segments <- function(value, caller = "unknown") {
  sequence_segments_rv(value)
}

