# selector tab implementation
# Core tab for selecting and ordering segments via interactive graph

# reactive values
selected_graph_segment <- reactiveVal(NULL)
selected_graph_edge <- reactiveVal(NULL)
selected_sequence_index <- reactiveVal(NULL)
graph_segments <- reactiveVal(character())
sequence_segments <- reactiveVal(character())
min_support <- reactiveVal(cache_get_if_exists("selector.min_support", 0))
min_percent <- reactiveVal(cache_get_if_exists("selector.min_percent", 0))

# data loading functions
load_bin_segment_table <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly, null.on.missing = TRUE)
}

load_seg_adj_count <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  get_data("BINNING_SEG_ADJ_count", tag = assembly, null.on.missing = TRUE)
}

load_seg_adj_total <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  get_data("BINNING_SEG_ADJ_total_read_count", tag = assembly, null.on.missing = TRUE)
}

load_seg_adj_associated <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  get_data("BINNING_SEG_ADJ_associated_read_count", tag = assembly, null.on.missing = TRUE)
}

# load segments from bin into graph
load_bin_segments <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  seg_table <- load_bin_segment_table(state$assembly)
  if (is.null(seg_table)) {
    return(NULL)
  }
  
  if (!"bin" %in% names(seg_table) || !"segment" %in% names(seg_table)) {
    return(NULL)
  }
  
  bin_segments <- seg_table[seg_table$bin == bin_id, ]
  if (nrow(bin_segments) == 0) {
    return(NULL)
  }
  
  return(bin_segments$segment)
}

# add neighbors of segment to graph
add_neighbors <- function(segment_id) {
  if (is.null(segment_id) || segment_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  count_mat <- load_seg_adj_count(state$assembly)
  if (is.null(count_mat)) {
    return(NULL)
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(NULL)
  }
  
  # find neighbors where segment is source
  neighbors_src <- count_mat[count_mat$seg_src == segment_id, ]
  neighbors_from_src <- unique(neighbors_src$seg_tgt)
  
  # find neighbors where segment is target
  neighbors_tgt <- count_mat[count_mat$seg_tgt == segment_id, ]
  neighbors_from_tgt <- unique(neighbors_tgt$seg_src)
  
  # combine all neighbors
  all_neighbors <- unique(c(neighbors_from_src, neighbors_from_tgt))
  
  return(all_neighbors)
}

# aggregate edge metrics across libraries
aggregate_edge_metrics <- function(edge_rows, count_mat, total_mat, associated_mat) {
  if (is.null(edge_rows) || nrow(edge_rows) == 0) {
    return(NULL)
  }
  
  # get library columns (exclude metadata columns)
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  
  # find library columns in count_mat
  lib_cols <- names(count_mat)[!names(count_mat) %in% metadata_cols]
  
  # aggregate metrics for each edge
  result <- data.frame(
    seg_src = edge_rows$seg_src,
    seg_tgt = edge_rows$seg_tgt,
    side_src = edge_rows$side_src,
    side_tgt = edge_rows$side_tgt,
    support = numeric(nrow(edge_rows)),
    total = numeric(nrow(edge_rows)),
    total_associated = numeric(nrow(edge_rows)),
    percent = numeric(nrow(edge_rows)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(edge_rows))) {
    seg_src <- edge_rows$seg_src[i]
    seg_tgt <- edge_rows$seg_tgt[i]
    side_src <- edge_rows$side_src[i]
    side_tgt <- edge_rows$side_tgt[i]
    
    # find matching rows in all three matrices
    count_row <- count_mat[count_mat$seg_src == seg_src & 
                          count_mat$seg_tgt == seg_tgt &
                          count_mat$side_src == side_src &
                          count_mat$side_tgt == side_tgt, ]
    total_row <- total_mat[total_mat$seg_src == seg_src & 
                          total_mat$seg_tgt == seg_tgt &
                          total_mat$side_src == side_src &
                          total_mat$side_tgt == side_tgt, ]
    associated_row <- associated_mat[associated_mat$seg_src == seg_src & 
                                    associated_mat$seg_tgt == seg_tgt &
                                    associated_mat$side_src == side_src &
                                    associated_mat$side_tgt == side_tgt, ]
    
    # aggregate across libraries
    support_sum <- 0
    total_sum <- 0
    associated_sum <- 0
    
    if (nrow(count_row) > 0) {
      for (col in lib_cols) {
        if (col %in% names(count_row)) {
          support_sum <- support_sum + as.numeric(count_row[[col]][1])
        }
      }
    }
    
    if (nrow(total_row) > 0) {
      for (col in lib_cols) {
        if (col %in% names(total_row)) {
          total_sum <- total_sum + as.numeric(total_row[[col]][1])
        }
      }
    }
    
    if (nrow(associated_row) > 0) {
      for (col in lib_cols) {
        if (col %in% names(associated_row)) {
          associated_sum <- associated_sum + as.numeric(associated_row[[col]][1])
        }
      }
    }
    
    result$support[i] <- support_sum
    result$total[i] <- total_sum
    result$total_associated[i] <- associated_sum
    result$percent[i] <- if (associated_sum > 0) (support_sum / associated_sum) * 100 else 0
  }
  
  return(result)
}

# build graph nodes
build_graph_nodes <- function(segment_ids) {
  if (length(segment_ids) == 0) {
    return(data.frame(id = character(), label = character(), stringsAsFactors = FALSE))
  }
  
  nodes <- data.frame(
    id = segment_ids,
    label = segment_ids,
    stringsAsFactors = FALSE
  )
  
  return(nodes)
}

# build graph edges with filtering
build_graph_edges <- function(segment_ids, count_mat, total_mat, associated_mat, bin_segment_table, min_support_val, min_percent_val) {
  if (is.null(count_mat) || is.null(total_mat) || is.null(associated_mat)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  if (length(segment_ids) == 0) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  if (!"seg_src" %in% names(count_mat) || !"seg_tgt" %in% names(count_mat)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # filter to edges where both src and tgt are in graph_segments
  edge_rows <- count_mat[count_mat$seg_src %in% segment_ids & 
                        count_mat$seg_tgt %in% segment_ids, ]
  
  if (nrow(edge_rows) == 0) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # aggregate metrics
  edge_metrics <- aggregate_edge_metrics(edge_rows, count_mat, total_mat, associated_mat)
  
  if (is.null(edge_metrics)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # add bin information
  if (!is.null(bin_segment_table) && "segment" %in% names(bin_segment_table) && "bin" %in% names(bin_segment_table)) {
    # create segment -> bin mapping
    seg_to_bin <- setNames(bin_segment_table$bin, bin_segment_table$segment)
    
    # add bin columns
    edge_metrics$bin_src <- seg_to_bin[edge_metrics$seg_src]
    edge_metrics$bin_tgt <- seg_to_bin[edge_metrics$seg_tgt]
    
    # handle missing bins (NA)
    edge_metrics$bin_src[is.na(edge_metrics$bin_src)] <- "N/A"
    edge_metrics$bin_tgt[is.na(edge_metrics$bin_tgt)] <- "N/A"
  } else {
    edge_metrics$bin_src <- "N/A"
    edge_metrics$bin_tgt <- "N/A"
  }
  
  # apply filters
  edge_metrics <- edge_metrics[edge_metrics$support >= min_support_val & 
                              edge_metrics$percent >= min_percent_val, ]
  
  if (nrow(edge_metrics) == 0) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  
  # build edges data frame
  edges <- data.frame(
    from = edge_metrics$seg_src,
    to = edge_metrics$seg_tgt,
    stringsAsFactors = FALSE
  )
  
  # add title for tooltip
  edges$title <- sprintf(
    "%s -> %s\nBin: %s -> %s\nSide: %s -> %s\nSupport: %.0f\nTotal: %.0f\nTotal Associated: %.0f\nPercent: %.2f%%",
    edge_metrics$seg_src,
    edge_metrics$seg_tgt,
    edge_metrics$bin_src,
    edge_metrics$bin_tgt,
    edge_metrics$side_src,
    edge_metrics$side_tgt,
    edge_metrics$support,
    edge_metrics$total,
    edge_metrics$total_associated,
    edge_metrics$percent
  )
  
  # store metrics in edge data for table
  edges$..support.. <- edge_metrics$support
  edges$..total.. <- edge_metrics$total
  edges$..total_associated.. <- edge_metrics$total_associated
  edges$..percent.. <- edge_metrics$percent
  edges$..side_src.. <- edge_metrics$side_src
  edges$..side_tgt.. <- edge_metrics$side_tgt
  edges$..bin_src.. <- edge_metrics$bin_src
  edges$..bin_tgt.. <- edge_metrics$bin_tgt
  
  return(edges)
}

# build edge table
build_edge_table <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0) {
    return(data.frame(
      seg_src = character(),
      seg_tgt = character(),
      bin_src = character(),
      bin_tgt = character(),
      side_src = character(),
      side_tgt = character(),
      support = numeric(),
      total = numeric(),
      total_associated = numeric(),
      percent = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  table_data <- data.frame(
    seg_src = edges$from,
    seg_tgt = edges$to,
    bin_src = edges$..bin_src..,
    bin_tgt = edges$..bin_tgt..,
    side_src = edges$..side_src..,
    side_tgt = edges$..side_tgt..,
    support = edges$..support..,
    total = edges$..total..,
    total_associated = edges$..total_associated..,
    percent = edges$..percent..,
    stringsAsFactors = FALSE
  )
  
  return(table_data)
}

# reset sequence from current plotted segments
reset_sequence <- function() {
  current_segments <- cxt_get_segments()
  if (length(current_segments) > 0) {
    sequence_segments(current_segments)
  } else {
    sequence_segments(character())
  }
}

# update main view with sequence segments
update_main_view <- function() {
  seg_ids <- sequence_segments()
  if (length(seg_ids) == 0) {
    cxt_set_view(NULL)
    return()
  }
  
  if (is.null(state$assembly)) {
    return()
  }
  
  seg_table <- load_bin_segment_table(state$assembly)
  if (is.null(seg_table)) {
    return()
  }
  
  if (!all(c("segment", "contig", "start", "end") %in% names(seg_table))) {
    return()
  }
  
  # filter to segments in sequence
  selected_segments <- seg_table[seg_table$segment %in% seg_ids, ]
  
  if (nrow(selected_segments) == 0) {
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
  
  cxt_set_view(result)
}

# UI function (will be called from server_tabs.r)
selector_tab_ui <- function() {
  tabPanel(
    "Selector",
    fluidRow(
      column(12,
        wellPanel(
          h5("Controls"),
          fluidRow(
            column(3,
              textInput("selectorBinInput", "Bin:", 
                       value = "",
                       placeholder = "Enter bin ID (e.g., b1)")
            ),
            column(3,
              actionButton("selectorLoadBinBtn", "Load Bin", 
                         class = "btn btn-primary")
            ),
            column(3,
              actionButton("selectorAddNeighborsBtn", "Add Neighbors", 
                         class = "btn btn-default")
            ),
            column(3,
              actionButton("selectorResetBtn", "Reset", 
                         class = "btn btn-default"),
              actionButton("selectorUpdateBtn", "Update", 
                         class = "btn btn-primary")
            )
          ),
          fluidRow(
            column(3,
              numericInput("selectorMinSupport", "Min Support:", 
                         value = cache_get_if_exists("selector.min_support", 0),
                         min = 0, step = 1)
            ),
            column(3,
              numericInput("selectorMinPercent", "Min Percent:", 
                         value = cache_get_if_exists("selector.min_percent", 0),
                         min = 0, max = 100, step = 0.1)
            ),
            column(6,
              h6("Sequence Actions:"),
              actionButton("selectorAddAfterBtn", "Add After", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorAddBeforeBtn", "Add Before", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorAddFirstBtn", "Add First", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorAddLastBtn", "Add Last", 
                         class = "btn btn-sm btn-default"),
              actionButton("selectorRemoveBtn", "Remove", 
                         class = "btn btn-sm btn-default")
            )
          )
        )
      )
    ),
    fluidRow(
      column(8,
        h4("Graph (click nodes/edges to select)"),
        visNetwork::visNetworkOutput("selectorGraph", height = "500px")
      ),
      column(4,
        h4("Sequence (drag to reorder, click ✕ to remove)"),
        uiOutput("selectorSequenceUI")
      )
    ),
    fluidRow(
      column(12,
        h4("Edge Data"),
        DTOutput("selectorEdgeTable")
      )
    )
  )
}

# observers and renderers (executed when sourced in server context)
  # observer for load bin button
  observeEvent(input$selectorLoadBinBtn, {
    bin_id <- trimws(input$selectorBinInput)
    if (bin_id == "") {
      showNotification("Please enter a bin ID", type = "warning")
      return()
    }
    
    segments <- load_bin_segments(bin_id)
    if (is.null(segments) || length(segments) == 0) {
      showNotification(sprintf("No segments found for bin: %s", bin_id), type = "warning")
      return()
    }
    
    graph_segments(segments)
    showNotification(sprintf("Loaded %d segments from bin %s", length(segments), bin_id), type = "message")
  })
  
  # observer for add neighbors button
  observeEvent(input$selectorAddNeighborsBtn, {
    seg_id <- selected_graph_segment()
    if (is.null(seg_id) || seg_id == "") {
      showNotification("Please select a segment in the graph first", type = "warning")
      return()
    }
    
    neighbors <- add_neighbors(seg_id)
    if (is.null(neighbors) || length(neighbors) == 0) {
      showNotification(sprintf("No neighbors found for segment: %s", seg_id), type = "warning")
      return()
    }
    
    current_segments <- graph_segments()
    new_segments <- unique(c(current_segments, neighbors))
    graph_segments(new_segments)
    showNotification(sprintf("Added %d neighbors to graph", length(setdiff(new_segments, current_segments))), type = "message")
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
    showNotification("Sequence reset to current plotted segments", type = "message")
  })
  
  # observer for update button
  observeEvent(input$selectorUpdateBtn, {
    update_main_view()
    showNotification("Updated main view with sequence segments", type = "message")
  })
  
  # observer for add after button
  observeEvent(input$selectorAddAfterBtn, {
    seg_id <- selected_graph_segment()
    seq_idx <- selected_sequence_index()
    
    if (is.null(seg_id) || seg_id == "") {
      showNotification("Please select a segment in the graph", type = "warning")
      return()
    }
    
    current_seq <- sequence_segments()
    if (length(current_seq) == 0) {
      sequence_segments(seg_id)
      return()
    }
    
    if (is.null(seq_idx) || seq_idx < 1 || seq_idx > length(current_seq)) {
      # add to end if no valid selection
      sequence_segments(c(current_seq, seg_id))
    } else {
      # insert after selected index
      new_seq <- c(current_seq[1:seq_idx], seg_id, current_seq[(seq_idx+1):length(current_seq)])
      sequence_segments(new_seq)
    }
  })
  
  # observer for add before button
  observeEvent(input$selectorAddBeforeBtn, {
    seg_id <- selected_graph_segment()
    seq_idx <- selected_sequence_index()
    
    if (is.null(seg_id) || seg_id == "") {
      showNotification("Please select a segment in the graph", type = "warning")
      return()
    }
    
    current_seq <- sequence_segments()
    if (length(current_seq) == 0) {
      sequence_segments(seg_id)
      return()
    }
    
    if (is.null(seq_idx) || seq_idx < 1 || seq_idx > length(current_seq)) {
      # add to start if no valid selection
      sequence_segments(c(seg_id, current_seq))
    } else {
      # insert before selected index
      new_seq <- c(current_seq[1:(seq_idx-1)], seg_id, current_seq[seq_idx:length(current_seq)])
      sequence_segments(new_seq)
    }
  })
  
  # observer for add first button
  observeEvent(input$selectorAddFirstBtn, {
    seg_id <- selected_graph_segment()
    if (is.null(seg_id) || seg_id == "") {
      showNotification("Please select a segment in the graph", type = "warning")
      return()
    }
    
    current_seq <- sequence_segments()
    sequence_segments(c(seg_id, current_seq))
  })
  
  # observer for add last button
  observeEvent(input$selectorAddLastBtn, {
    seg_id <- selected_graph_segment()
    if (is.null(seg_id) || seg_id == "") {
      showNotification("Please select a segment in the graph", type = "warning")
      return()
    }
    
    current_seq <- sequence_segments()
    sequence_segments(c(current_seq, seg_id))
  })
  
  # observer for remove button
  observeEvent(input$selectorRemoveBtn, {
    seq_idx <- selected_sequence_index()
    if (is.null(seq_idx) || seq_idx < 1) {
      showNotification("Please select a segment in the sequence", type = "warning")
      return()
    }
    
    current_seq <- sequence_segments()
    if (seq_idx > length(current_seq)) {
      return()
    }
    
    new_seq <- current_seq[-seq_idx]
    sequence_segments(new_seq)
    selected_sequence_index(NULL)
  })
  
  # observer for sequence reorder
  observeEvent(input$selectorSequenceOrder, {
    new_order <- as.numeric(input$selectorSequenceOrder)
    if (length(new_order) > 0 && length(new_order) == length(sequence_segments())) {
      current_seq <- sequence_segments()
      sequence_segments(current_seq[new_order])
    }
  })
  
  # observer for remove from sequence
  observeEvent(input$selectorRemoveFromSequence, {
    seg_id <- input$selectorRemoveFromSequence
    if (is.null(seg_id) || seg_id == "") {
      return()
    }
    
    current_seq <- sequence_segments()
    new_seq <- current_seq[current_seq != seg_id]
    sequence_segments(new_seq)
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
        dragView = TRUE
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
    network$x$options$interaction$touch <- FALSE
    
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
    seg_ids <- sequence_segments()
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

