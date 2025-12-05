# graph state management
# reactive values for the association graph

# graph selection state
selected_graph_segments <- reactiveVal(character())
selected_graph_edge <- reactiveVal(NULL)

# trigger to force graph re-render when selection changes
graph_selection_trigger <- reactiveVal(0)

# hover state
hovered_node <- reactiveVal(NULL)
hovered_edge <- reactiveVal(NULL)

# filter parameters
min_support <- reactiveVal(cache_get_if_exists("selector.min_support", 0))
min_percent <- reactiveVal(cache_get_if_exists("selector.min_percent", 0))

# base segments for graph (before neighbor expansion)
graph_base_segments <- reactiveVal(character())

