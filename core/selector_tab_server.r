# selector tab server
# action functions, observers and renderers

# reset sequence from current plotted segments
reset_sequence <- function() {
  current_segments <- cxt_get_segments()
  if (length(current_segments) > 0) {
    set_sequence_segments(current_segments, "reset_sequence")
  } else {
    set_sequence_segments(character(), "reset_sequence_empty")
  }
}

# update main view with sequence segments
update_main_view <- function() {
  cat("[update_main_view] starting\n")
  seg_ids <- get_sequence_segments()
  if (length(seg_ids) == 0) {
    cat("[update_main_view] no segments, clearing view\n")
    state$segments <- data.frame(
      segment = character(), contig = character(),
      start = integer(), end = integer(), stringsAsFactors = FALSE)
    return()
  }
  
  if (is.null(state$assembly)) {
    cat("[update_main_view] assembly is NULL, returning\n")
    return()
  }
  
  seg_table <- load_bin_segment_table(state$assembly)
  if (is.null(seg_table)) {
    cat("[update_main_view] seg_table is NULL, returning\n")
    return()
  }
  
  if (!all(c("segment", "contig", "start", "end") %in% names(seg_table))) {
    cat("[update_main_view] seg_table missing required columns\n")
    return()
  }
  
  # filter to segments in sequence
  selected_segments <- seg_table[seg_table$segment %in% seg_ids, ]
  cat(sprintf("[update_main_view] found %d matching segments out of %d\n", 
              nrow(selected_segments), length(seg_ids)))
  
  if (nrow(selected_segments) == 0) {
    cat("[update_main_view] no matching segments found\n")
    return()
  }
  
  # reorder to match sequence order
  seg_order <- match(selected_segments$segment, seg_ids)
  selected_segments <- selected_segments[order(seg_order), ]
  
  # get required columns
  result <- data.frame(
    segment = selected_segments$segment,
    contig = selected_segments$contig,
    start = selected_segments$start,
    end = selected_segments$end,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("[update_main_view] setting state$segments with %d segments\n", nrow(result)))
  state$segments <- result
}

# observer for load bin button
observeEvent(input$selectorLoadBinBtn, {
  bin_id <- trimws(input$selectorBinInput)
  if (bin_id == "")
    return()
  
  segments <- load_bin_segments(bin_id)
  if (is.null(segments) || length(segments) == 0)
    return()
  
  graph_segments(segments)
})

# observer for add neighbors button
observeEvent(input$selectorAddNeighborsBtn, {
  seg_id <- selected_graph_segment()
  if (is.null(seg_id) || seg_id == "")
    return()
  
  neighbors <- add_neighbors(seg_id)
  if (is.null(neighbors) || length(neighbors) == 0)
    return()
  
  current_segments <- graph_segments()
  new_segments <- unique(c(current_segments, neighbors))
  graph_segments(new_segments)
})

# observer for min support filter
observeEvent(input$selectorMinSupport, {
  min_support(input$selectorMinSupport)
  cache_set("selector.min_support", input$selectorMinSupport)
})

# observer for min percent filter
observeEvent(input$selectorMinPercent, {
  min_percent(input$selectorMinPercent)
  cache_set("selector.min_percent", input$selectorMinPercent)
})

# observer for reset button
observeEvent(input$selectorResetBtn, {
  reset_sequence()
})

# observer for update button
observeEvent(input$selectorUpdateBtn, {
  update_main_view()
})

# observer for add after button
observeEvent(input$selectorAddAfterBtn, {
  seg_id <- selected_graph_segment()
  seq_idx <- selected_sequence_index()
  
  if (is.null(seg_id) || seg_id == "")
    return()
  
  current_seq <- get_sequence_segments()
  if (length(current_seq) == 0) {
    set_sequence_segments(seg_id, "add_after_empty")
    return()
  }
  
  if (is.null(seq_idx) || seq_idx < 1 || seq_idx > length(current_seq)) {
    set_sequence_segments(c(current_seq, seg_id), "add_after_end")
  } else {
    new_seq <- c(current_seq[1:seq_idx], seg_id, current_seq[(seq_idx+1):length(current_seq)])
    set_sequence_segments(new_seq, "add_after_idx")
  }
})

# observer for add before button
observeEvent(input$selectorAddBeforeBtn, {
  seg_id <- selected_graph_segment()
  seq_idx <- selected_sequence_index()
  
  if (is.null(seg_id) || seg_id == "")
    return()
  
  current_seq <- get_sequence_segments()
  if (length(current_seq) == 0) {
    set_sequence_segments(seg_id, "add_before_empty")
    return()
  }
  
  if (is.null(seq_idx) || seq_idx < 1 || seq_idx > length(current_seq)) {
    set_sequence_segments(c(seg_id, current_seq), "add_before_start")
  } else {
    new_seq <- c(current_seq[1:(seq_idx-1)], seg_id, current_seq[seq_idx:length(current_seq)])
    set_sequence_segments(new_seq, "add_before_idx")
  }
})

# observer for add first button
observeEvent(input$selectorAddFirstBtn, {
  seg_id <- selected_graph_segment()
  if (is.null(seg_id) || seg_id == "")
    return()
  
  current_seq <- get_sequence_segments()
  set_sequence_segments(c(seg_id, current_seq), "add_first")
})

# observer for add last button
observeEvent(input$selectorAddLastBtn, {
  seg_id <- selected_graph_segment()
  if (is.null(seg_id) || seg_id == "")
    return()
  
  current_seq <- get_sequence_segments()
  set_sequence_segments(c(current_seq, seg_id), "add_last")
})

# observer for remove button
observeEvent(input$selectorRemoveBtn, {
  seq_idx <- selected_sequence_index()
  if (is.null(seq_idx) || seq_idx < 1)
    return()
  
  current_seq <- get_sequence_segments()
  if (seq_idx > length(current_seq)) {
    return()
  }
  
  new_seq <- current_seq[-seq_idx]
  set_sequence_segments(new_seq, "remove_btn")
  selected_sequence_index(NULL)
})

# observer for sequence reorder
observeEvent(input$selectorSequenceOrder, {
  new_order <- input$selectorSequenceOrder
  current_seq <- get_sequence_segments()
  # rank_list returns segment IDs in the new order
  if (length(new_order) > 0 && length(new_order) == length(current_seq)) {
    set_sequence_segments(new_order, "reorder")
  }
})

# observer for remove from sequence
observeEvent(input$selectorRemoveFromSequence, {
  seg_id <- input$selectorRemoveFromSequence
  if (is.null(seg_id) || seg_id == "") {
    return()
  }
  
  current_seq <- get_sequence_segments()
  new_seq <- current_seq[current_seq != seg_id]
  set_sequence_segments(new_seq, "remove_from_seq")
})

# render graph
output$selectorGraph <- visNetwork::renderVisNetwork({
  seg_ids <- graph_segments()
  if (length(seg_ids) == 0) {
    return(visNetwork::visNetwork(
      data.frame(id = character(), label = character()),
      data.frame(from = character(), to = character())
    ) %>% visNetwork::visNodes(color = list(background = "#DDE5FF")))
  }
  
  # load data
  count_mat <- load_seg_adj_count(state$assembly)
  total_mat <- load_seg_adj_total(state$assembly)
  associated_mat <- load_seg_adj_associated(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  
  # build nodes
  nodes <- build_graph_nodes(seg_ids)
  
  # build edges
  edges <- build_graph_edges(seg_ids, count_mat, total_mat, associated_mat, bin_segment_table,
                            min_support(), min_percent())
  
  # highlight selected node
  selected_node <- selected_graph_segment()
  if (!is.null(selected_node) && selected_node %in% nodes$id) {
    nodes$color.background <- ifelse(nodes$id == selected_node, "#FF6B6B", "#DDE5FF")
    nodes$color.border <- ifelse(nodes$id == selected_node, "#C92A2A", "#4A90E2")
  } else {
    nodes$color.background <- "#DDE5FF"
    nodes$color.border <- "#4A90E2"
  }
  
  # highlight selected edge
  selected_edge <- selected_graph_edge()
  if (!is.null(selected_edge) && "from" %in% names(selected_edge) && "to" %in% names(selected_edge)) {
    edges$color.color <- ifelse(edges$from == selected_edge$from & edges$to == selected_edge$to,
                                "#FF6B6B", "#848484")
    edges$color.highlight <- ifelse(edges$from == selected_edge$from & edges$to == selected_edge$to,
                                   "#FF6B6B", "#848484")
  } else {
    edges$color.color <- "#848484"
    edges$color.highlight <- "#848484"
  }
  
  # create network
  network <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visEdges(arrows = "to") %>%
    visNetwork::visInteraction(
      zoomView = TRUE,
      dragView = TRUE,
      multiselect = FALSE,
      keyboard = FALSE
    ) %>%
    visNetwork::visOptions(
      manipulation = FALSE
    ) %>%
    visNetwork::visEvents(
      select = paste0("
        function(params) {
          if (params.nodes.length > 0) {
            Shiny.setInputValue('selectorGraphNodeClick', params.nodes[0], {priority: 'event'});
          }
          if (params.edges.length > 0) {
            var edge = params.edges[0];
            Shiny.setInputValue('selectorGraphEdgeClick', {from: edge.from, to: edge.to}, {priority: 'event'});
          }
        }
      ")
    )
  
  # disable pinch zoom by setting touch-related options
  network$x$options$interaction$zoomSpeed <- 0
  network$x$options$interaction$zoomView <- FALSE
  
  return(network)
})

# observer for graph node click
observeEvent(input$selectorGraphNodeClick, {
  selected_graph_segment(input$selectorGraphNodeClick)
})

# observer for graph edge click
observeEvent(input$selectorGraphEdgeClick, {
  selected_graph_edge(input$selectorGraphEdgeClick)
})

# render sequence UI
output$selectorSequenceUI <- renderUI({
  seg_ids <- get_sequence_segments()
  if (length(seg_ids) == 0) {
    return(helpText("Sequence is empty. Click graph nodes to add segments."))
  }
  
  # build UI boxes for each segment
  selected_idx <- selected_sequence_index()
  list_items <- lapply(seq_along(seg_ids), function(i) {
    seg_id <- seg_ids[i]
    is_selected <- !is.null(selected_idx) && selected_idx == i
    bg_color <- if (is_selected) "#FFE5E5" else "#EEF2FF"
    tags$div(
      `data-index` = i,
      style = sprintf("padding:6px; margin:4px; background:%s; border-radius:6px; display:flex; justify-content:space-between; align-items:center; cursor:pointer;", bg_color),
      onclick = sprintf("Shiny.setInputValue('selectorSequenceClick', %d, {priority: 'event'})", i),
      tags$span(seg_id, style = "flex: 1;"),
      tags$button(
        "✕",
        type = "button",
        class = "btn btn-sm",
        style = "padding:0 6px; margin-left:10px;",
        onclick = sprintf(
          "Shiny.setInputValue('selectorRemoveFromSequence', '%s', {priority: 'event'}); event.stopPropagation(); return false;",
          seg_id
        )
      )
    )
  })
  
  names(list_items) <- as.character(seg_ids)
  
  sortable::rank_list(
    text = "Drag to reorder:",
    labels = list_items,
    input_id = "selectorSequenceOrder"
  )
})

# observer for sequence click
observeEvent(input$selectorSequenceClick, {
  selected_sequence_index(input$selectorSequenceClick)
})

# render edge table
output$selectorEdgeTable <- renderDT({
  seg_ids <- graph_segments()
  if (length(seg_ids) == 0) {
    return(datatable(
      data.frame(Message = "No segments loaded. Load a bin to see edges."),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # load data
  count_mat <- load_seg_adj_count(state$assembly)
  total_mat <- load_seg_adj_total(state$assembly)
  associated_mat <- load_seg_adj_associated(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(datatable(
      data.frame(Message = "Association data not available."),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # build edges
  edges <- build_graph_edges(seg_ids, count_mat, total_mat, associated_mat, bin_segment_table,
                            min_support(), min_percent())
  
  if (is.null(edges) || nrow(edges) == 0) {
    return(datatable(
      data.frame(Message = "No edges found (may be filtered out)."),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # build table
  table_data <- build_edge_table(edges)
  
  # highlight selected edge
  selected_edge <- selected_graph_edge()
  if (!is.null(selected_edge) && "from" %in% names(selected_edge) && "to" %in% names(selected_edge)) {
    selected_rows <- which(table_data$seg_src == selected_edge$from & 
                          table_data$seg_tgt == selected_edge$to)
  } else {
    selected_rows <- NULL
  }
  
  datatable(
    table_data,
    rownames = FALSE,
    selection = list(mode = "single", selected = if (length(selected_rows) > 0) selected_rows[1] - 1 else NULL),
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip"
    ),
    filter = "top"
  ) %>% formatRound(columns = c("support", "total", "total_associated"), digits = 0) %>%
    formatRound(columns = c("percent"), digits = 2)
})

# observer for edge table row click
observeEvent(input$selectorEdgeTable_rows_selected, {
  if (is.null(input$selectorEdgeTable_rows_selected)) {
    return()
  }
  
  seg_ids <- graph_segments()
  count_mat <- load_seg_adj_count(state$assembly)
  total_mat <- load_seg_adj_total(state$assembly)
  associated_mat <- load_seg_adj_associated(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  
  edges <- build_graph_edges(seg_ids, count_mat, total_mat, associated_mat, bin_segment_table,
                            min_support(), min_percent())
  
  if (is.null(edges) || nrow(edges) == 0) {
    return()
  }
  
  row_idx <- input$selectorEdgeTable_rows_selected
  if (row_idx > 0 && row_idx <= nrow(edges)) {
    selected_graph_edge(list(from = edges$from[row_idx], to = edges$to[row_idx]))
  }
})

