# graph server
# association graph server logic

#########################################################################
# tracking controls
#########################################################################

# update to view button
observeEvent(input$selectorUpdateToView, {
  state_segs <- state$segments
  if (!is.null(state_segs) && nrow(state_segs) > 0) {
    graph_base_segments(state_segs$segment)
  } else {
    graph_base_segments(character())
  }
})

# update to zoom button
observeEvent(input$selectorUpdateToZoom, {
  zoom_segs <- cxt_get_segments(limit_to_zoom = TRUE)
  if (length(zoom_segs) > 0) {
    graph_base_segments(zoom_segs)
  }
})

# auto-tracking observer
observe({
  auto_mode <- input$selectorAutoUpdate
  if (is.null(auto_mode) || auto_mode == "off") return()
  
  if (auto_mode == "view") {
    state_segs <- state$segments
    if (!is.null(state_segs) && nrow(state_segs) > 0) {
      graph_base_segments(state_segs$segment)
    } else {
      graph_base_segments(character())
    }
  } else if (auto_mode == "zoom") {
    zoom_segs <- cxt_get_segments(limit_to_zoom = TRUE)
    if (length(zoom_segs) > 0) {
      graph_base_segments(zoom_segs)
    }
  }
})

#########################################################################
# computed graph segments from base + neighbors
#########################################################################

computed_graph_segments <- reactive({
  base_segments <- graph_base_segments()
  if (length(base_segments) == 0) {
    return(character())
  }
  
  depth <- if (is.null(input$selectorNeighborDepth)) 0 else input$selectorNeighborDepth
  
  if (depth <= 0) {
    return(base_segments)
  }
  
  # get filter values for neighbor expansion
  min_support_val <- min_support()
  min_percent_val <- min_percent()
  
  # get edge library selection
  edge_lib <- if (!is.null(input$graphEdgeLibs)) input$graphEdgeLibs else "all"
  if (edge_lib == "all") edge_lib <- NULL
  
  # expand neighbors iteratively along filtered edges only
  all_segments <- base_segments
  current_frontier <- base_segments
  
  for (i in seq_len(depth)) {
    new_neighbors <- character()
    for (seg_id in current_frontier) {
      neighbors <- add_filtered_neighbors(seg_id, min_support_val, min_percent_val, lib = edge_lib)
      if (!is.null(neighbors) && length(neighbors) > 0) {
        new_neighbors <- c(new_neighbors, neighbors)
      }
    }
    new_neighbors <- unique(new_neighbors)
    new_neighbors <- new_neighbors[!new_neighbors %in% all_segments]
    if (length(new_neighbors) == 0) break
    all_segments <- c(all_segments, new_neighbors)
    current_frontier <- new_neighbors
  }
  
  return(unique(all_segments))
})

#########################################################################
# filter observers
#########################################################################

observeEvent(input$selectorMinSupport, {
  min_support(input$selectorMinSupport)
  cache_set("selector.min_support", input$selectorMinSupport)
})

observeEvent(input$selectorMinPercent, {
  min_percent(input$selectorMinPercent)
  cache_set("selector.min_percent", input$selectorMinPercent)
})

observeEvent(input$selectorNeighborDepth, {
  cache_set("selector.neighbor_depth", input$selectorNeighborDepth)
})

#########################################################################
# graph action buttons (write directly to state$segments)
#########################################################################

# extracts segment IDs from side-node IDs (strips _L or _R suffix)
extract_segments_from_nodes <- function(node_ids) {
  if (length(node_ids) == 0) return(character())
  segs <- sub("_[LR]$", "", node_ids)
  unique(segs)
}

# goto: navigate to selected segments in main view
observeEvent(input$selectorGotoBtn, {
  node_ids <- selected_graph_segments()
  if (length(node_ids) == 0) return()
  
  seg_ids <- extract_segments_from_nodes(node_ids)
  if (length(seg_ids) == 0) return()
  
  if (is.null(state$assembly)) return()
  seg_table <- get_segments(state$assembly)
  if (is.null(seg_table)) return()
  
  # use match to preserve order of seg_ids
  idx <- match(seg_ids, seg_table$segment)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return()
  selected_segments <- seg_table[idx, ]
  
  regions_module_output$push_undo_state()
  state$segments <- selected_segments
  
  # reset zoom to see full range of new segments
  state$zoom <- NULL
})

# add: append selected graph segments to state$segments
observeEvent(input$selectorAddBtn, {
  node_ids <- selected_graph_segments()
  if (length(node_ids) == 0) return()
  
  seg_ids <- extract_segments_from_nodes(node_ids)
  if (length(seg_ids) == 0) return()
  
  if (is.null(state$assembly)) return()
  seg_table <- get_segments(state$assembly)
  if (is.null(seg_table)) return()
  
  current_segs <- state$segments
  current_ids <- if (is.null(current_segs) || nrow(current_segs) == 0) character() else current_segs$segment
  
  new_ids <- seg_ids[!seg_ids %in% current_ids]
  if (length(new_ids) == 0) return()
  
  new_segments <- seg_table[seg_table$segment %in% new_ids, ]
  if (nrow(new_segments) == 0) return()
  
  regions_module_output$push_undo_state()
  
  if (is.null(current_segs) || nrow(current_segs) == 0) {
    state$segments <- new_segments
  } else {
    # use only columns present in current_segs
    common_cols <- intersect(names(current_segs), names(new_segments))
    state$segments <- rbind(current_segs[, common_cols, drop = FALSE], 
                           new_segments[, common_cols, drop = FALSE])
  }
  
  # reset zoom to see full range of new segments
  state$zoom <- NULL
})

# remove: remove selected graph segments from state$segments
observeEvent(input$selectorRemoveBtn, {
  node_ids <- selected_graph_segments()
  if (length(node_ids) == 0) return()
  
  seg_ids <- extract_segments_from_nodes(node_ids)
  if (length(seg_ids) == 0) return()
  
  current_segs <- state$segments
  if (is.null(current_segs) || nrow(current_segs) == 0) return()
  
  new_segs <- current_segs[!current_segs$segment %in% seg_ids, , drop = FALSE]
  
  regions_module_output$push_undo_state()
  state$segments <- new_segs
  
  # reset zoom to see full range of remaining segments
  state$zoom <- NULL
})

# select all graph nodes
observeEvent(input$selectorSelectAllBtn, {
  all_segs <- computed_graph_segments()
  all_nodes <- c(paste0(all_segs, "_L"), paste0(all_segs, "_R"))
  selected_graph_segments(all_nodes)
})

# clear graph selection
observeEvent(input$selectorClearSelectionBtn, {
  selected_graph_segments(character())
})

# add neighbors: expand graph to include neighbors of selected nodes
observeEvent(input$selectorAddNeighbors, {
  selected_nodes <- selected_graph_segments()
  if (length(selected_nodes) == 0 || is.null(state$assembly)) return()
  
  edge_lib <- if (!is.null(input$graphEdgeLibs)) input$graphEdgeLibs else "all"
  if (edge_lib == "all") edge_lib <- NULL
  neighbor_nodes <- find_neighbor_nodes(selected_nodes, state$assembly, 
                                       min_support(), min_percent(), lib = edge_lib)
  
  if (length(neighbor_nodes) > 0) {
    neighbor_seg_ids <- unique(sub("_[LR]$", "", neighbor_nodes))
    current_base_segs <- graph_base_segments()
    graph_base_segments(unique(c(current_base_segs, neighbor_seg_ids)))
    selected_graph_segments(character())
  }
})

#########################################################################
# graph rendering
#########################################################################

# store base graph data (structure only, no visual properties) for proxy updates
base_graph_nodes <- reactiveVal(NULL)
base_graph_edges <- reactiveVal(NULL)
current_bin_segment_table <- reactiveVal(NULL)

output$selectorGraph <- visNetwork::renderVisNetwork({
  # dependencies: computed_graph_segments (which depends on input$graphEdgeLibs),
  # state$assembly, input$graphDirectedEdges, input$graphEdgeLibs (read directly below),
  # min_support(), min_percent()
  # NOTE: visual inputs (selection, colors, labels, fonts, edge metrics) are NOT dependencies
  seg_ids <- computed_graph_segments()
  
  if (length(seg_ids) == 0) {
    base_graph_nodes(NULL)
    base_graph_edges(NULL)
    return(visNetwork::visNetwork(
      data.frame(id = character(), label = character()),
      data.frame(from = character(), to = character())
    ))
  }
  
  # load data
  count_mat <- load_seg_adj_count(state$assembly)
  total_mat <- load_seg_adj_total(state$assembly)
  associated_mat <- load_seg_adj_associated(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  
  # build nodes (two per segment: left and right sides)
  nodes <- build_graph_nodes(seg_ids)
  
  # build edges (internal + inter-segment)
  # reading input$graphEdgeLibs here establishes reactive dependency
  directed <- if (!is.null(input$graphDirectedEdges)) input$graphDirectedEdges else FALSE
  edge_lib <- if (!is.null(input$graphEdgeLibs)) input$graphEdgeLibs else "all"
  if (edge_lib == "all") edge_lib <- NULL
  edges <- build_graph_edges(seg_ids, count_mat, total_mat, associated_mat, bin_segment_table,
                            min_support(), min_percent(), directed = directed, lib = edge_lib)
  
  # store base nodes/edges (structure only, before visual properties)
  base_graph_nodes(nodes)
  base_graph_edges(edges)
  current_bin_segment_table(bin_segment_table)
  
  # apply visual properties for initial render (using isolate to avoid reactive dependencies)
  # unified observer will handle all subsequent updates via proxy
  selected_nodes <- isolate(selected_graph_segments())
  ordered_segs <- if (!is.null(state$segments) && nrow(state$segments) > 0) state$segments$segment else character()
  seg_colors <- get_segment_colors(nodes$segment, state$assembly, ordered_segs)
  names(seg_colors) <- nodes$segment
  nodes$color.background <- ifelse(nodes$segment %in% names(seg_colors), 
                                   seg_colors[nodes$segment], "#E0E0E0")
  nodes$color.border <- "#666666"
  nodes$borderWidth <- 1
  
  # apply selection highlighting during initial rendering
  if (length(selected_nodes) > 0) {
    selected_mask <- nodes$id %in% selected_nodes
    nodes$color.border[selected_mask] <- "#C92A2A"
    nodes$borderWidth[selected_mask] <- 3
  }
  
  # set node labels based on label type
  label_type <- isolate(if (!is.null(input$graphNodeLabel)) input$graphNodeLabel else "index")
  nodes$label <- sapply(seq_len(nrow(nodes)), function(i) {
    seg <- nodes$segment[i]
    side <- nodes$side[i]
    if (label_type == "id") {
      return(paste0(seg, "_", side))
    }
    if (!seg %in% ordered_segs) return(" ")
    idx <- which(ordered_segs == seg)
    return(paste0(idx, "_", side))
  })
  
  # style edges: internal edges are thick and short, inter-segment edges are thinner
  # apply initial colors and labels based on current metric selection
  if (nrow(edges) > 0) {
    metric <- isolate(if (!is.null(input$graphEdgeMetric)) input$graphEdgeMetric else "none")
    show_labels <- isolate(if (!is.null(input$graphEdgeLabels)) input$graphEdgeLabels else FALSE)
    lib1 <- isolate(input$graphEdgeLib1)
    lib2 <- isolate(input$graphEdgeLib2)
    font_size <- isolate(if (!is.null(input$graphFontSize)) as.numeric(input$graphFontSize) else 38)
    edge_font_size <- font_size
    
    edge_colors <- compute_edge_colors(edges, metric, lib1, lib2)
    if (length(edge_colors) == nrow(edges)) {
      edges$color.color <- edge_colors
    } else {
      edges$color.color <- ifelse(edges$..is_internal.., "#000000", "#848484")
    }
    
    edge_labels <- compute_edge_labels(edges, metric, show_labels, lib1, lib2)
    if (length(edge_labels) == nrow(edges)) {
      edges$label <- edge_labels
      edges$font.size <- edge_font_size
    } else {
      edges$label <- ""
      edges$font.size <- edge_font_size
    }
    
    edges$width <- ifelse(edges$..is_internal.., 18, 6)
    edges$length <- ifelse(edges$..is_internal.., 15, 200)
    
    # add smooth curves for bidirectional edges to avoid overlap
    if (directed && nrow(edges) > 0) {
      # detect bidirectional edges (where A->B and B->A both exist)
      edge_pairs <- paste(edges$from, edges$to, sep = "|")
      reverse_pairs <- paste(edges$to, edges$from, sep = "|")
      is_bidirectional <- edge_pairs %in% reverse_pairs
      
      # apply smooth curves to bidirectional edges (visNetwork expects list format)
      edges$smooth <- I(lapply(seq_len(nrow(edges)), function(i) {
        if (is_bidirectional[i]) {
          list(enabled = TRUE, type = "curvedCW", roundness = 0.1)
        } else {
          FALSE
        }
      }))
      
      # alternate curve direction for reverse edges
      bidirectional_indices <- which(is_bidirectional)
      for (i in bidirectional_indices) {
        reverse_idx <- which(edges$from == edges$to[i] & edges$to == edges$from[i])
        if (length(reverse_idx) > 0) {
          edges$smooth[[reverse_idx[1]]] <- list(enabled = TRUE, type = "curvedCCW", roundness = 0.1)
        }
      }
    } else {
      edges$smooth <- FALSE
    }
    
    # internal edges always light gray and thick
    internal_mask <- !is.na(edges$..is_internal..) & edges$..is_internal..
    edges$color.color[internal_mask] <- "#D3D3D3"
    edges$label[internal_mask] <- ""
  }
  
  # compute layout using igraph's Fruchterman-Reingold algorithm
  if (nrow(edges) > 0 && requireNamespace("igraph", quietly = TRUE)) {
    # create igraph graph from edges (use directed based on input setting)
    g <- igraph::graph_from_data_frame(edges[, c("from", "to"), drop = FALSE], directed = directed, vertices = nodes)
    
    # compute Fruchterman-Reingold layout (force-directed)
    # set seed for reproducible layout
    set.seed(42)
    layout <- igraph::layout_with_fr(g)
    
    # apply layout coordinates to nodes (scale and invert y for vis.js)
    nodes$x <- layout[, 1] * 100
    nodes$y <- -layout[, 2] * 100  # invert for vis.js
  }
  
  nodes$font.color <- "black"
  
  network <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visNodes(
      font = list(size = font_size, color = "black")
    ) %>%
    visNetwork::visPhysics(
      enabled = FALSE
    ) %>%
    visNetwork::visInteraction(
      zoomView = FALSE,
      dragView = TRUE,
      multiselect = TRUE,
      keyboard = FALSE,
      hover = TRUE,
      hoverConnectedEdges = TRUE
    ) %>%
    visNetwork::visOptions(manipulation = FALSE) %>%
    visNetwork::visEvents(
      select = paste0("
        function(params) {
          if (params.nodes.length > 0) {
            Shiny.setInputValue('selectorGraphNodeSelect', params.nodes, {priority: 'event'});
          }
          if (params.edges.length > 0) {
            Shiny.setInputValue('selectorGraphEdgeSelect', params.edges, {priority: 'event'});
          }
        }
      "),
      hoverNode = paste0("
        function(params) {
          Shiny.setInputValue('selectorGraphNodeHover', params.node, {priority: 'event'});
        }
      "),
      blurNode = paste0("
        function(params) {
          Shiny.setInputValue('selectorGraphNodeHover', null, {priority: 'event'});
        }
      "),
      hoverEdge = paste0("
        function(params) {
          Shiny.setInputValue('selectorGraphEdgeHover', params.edge, {priority: 'event'});
        }
      "),
      blurEdge = paste0("
        function(params) {
          Shiny.setInputValue('selectorGraphEdgeHover', null, {priority: 'event'});
        }
      ")
    )
  
  network <- htmlwidgets::onRender(network, "
    function(el, x) {
      var network = this.network;
      network.once('stabilized', function() {
        network.fit();
      });
      Shiny.addCustomMessageHandler('selectorGraphZoom', function(data) {
        var scale = network.getScale() * data.factor;
        scale = Math.max(0.1, Math.min(scale, 10));
        network.moveTo({scale: scale});
      });
    }
  ")
  
  return(network)
})

# unified observer for all visual updates (selection, colors, labels, fonts, edges)
observe({
  # watch all visual inputs
  selected_nodes <- selected_graph_segments()
  scheme_val <- input$segmentColorScheme
  grayscale_val <- input$segmentGrayscale
  label_type <- input$graphNodeLabel
  font_size <- input$graphFontSize
  metric <- input$graphEdgeMetric
  show_labels <- input$graphEdgeLabels
  directed <- input$graphDirectedEdges
  lib1 <- input$graphEdgeLib1
  lib2 <- input$graphEdgeLib2
  
  # get base nodes/edges from reactiveVals
  nodes <- base_graph_nodes()
  edges <- base_graph_edges()
  
  if (is.null(nodes) || nrow(nodes) == 0) return()
  
  # compute node visual properties from base data
  ordered_segs <- if (!is.null(state$segments) && nrow(state$segments) > 0) state$segments$segment else character()
  seg_colors <- get_segment_colors(nodes$segment, state$assembly, ordered_segs)
  names(seg_colors) <- nodes$segment
  nodes$color.background <- ifelse(nodes$segment %in% names(seg_colors), 
                                   seg_colors[nodes$segment], "#E0E0E0")
  nodes$color.border <- "#666666"
  nodes$borderWidth <- 1
  nodes$font.color <- "black"
  
  # apply selection highlighting
  if (length(selected_nodes) > 0) {
    selected_mask <- nodes$id %in% selected_nodes
    nodes$color.border[selected_mask] <- "#C92A2A"
    nodes$borderWidth[selected_mask] <- 3
  }
  
  # set node labels based on label type
  label_type_val <- if (!is.null(label_type)) label_type else "index"
  nodes$label <- sapply(seq_len(nrow(nodes)), function(i) {
    seg <- nodes$segment[i]
    side <- nodes$side[i]
    if (label_type_val == "id") {
      return(paste0(seg, "_", side))
    }
    if (!seg %in% ordered_segs) return(" ")
    idx <- which(ordered_segs == seg)
    return(paste0(idx, "_", side))
  })
  
  # update font size via options
  font_size_val <- if (!is.null(font_size)) as.numeric(font_size) else 38
  visNetwork::visNetworkProxy("selectorGraph") %>%
    visNetwork::visSetOptions(options = list(
      nodes = list(font = list(size = font_size_val, color = "black"))
    ))
  
  # compute edge visual properties from base data
  if (!is.null(edges) && nrow(edges) > 0) {
    metric_val <- if (!is.null(metric)) metric else "none"
    show_labels_val <- if (!is.null(show_labels)) show_labels else FALSE
    lib1_val <- lib1
    lib2_val <- lib2
    edge_font_size <- font_size_val
    
    edge_colors <- compute_edge_colors(edges, metric_val, lib1_val, lib2_val)
    if (length(edge_colors) == nrow(edges)) {
      edges$color.color <- edge_colors
    } else {
      edges$color.color <- ifelse(edges$..is_internal.., "#000000", "#848484")
    }
    
    edge_labels <- compute_edge_labels(edges, metric_val, show_labels_val, lib1_val, lib2_val)
    if (length(edge_labels) == nrow(edges)) {
      edges$label <- edge_labels
    } else {
      edges$label <- ""
    }
    edges$font.size <- edge_font_size
    
    edges$width <- ifelse(edges$..is_internal.., 18, 6)
    edges$length <- ifelse(edges$..is_internal.., 15, 200)
    
    # update smooth curves for bidirectional edges when directed mode changes
    directed_val <- if (!is.null(directed)) directed else FALSE
    if (directed_val && nrow(edges) > 0) {
      # detect bidirectional edges
      edge_pairs <- paste(edges$from, edges$to, sep = "|")
      reverse_pairs <- paste(edges$to, edges$from, sep = "|")
      is_bidirectional <- edge_pairs %in% reverse_pairs
      
      # apply smooth curves to bidirectional edges (visNetwork expects list format)
      edges$smooth <- I(lapply(seq_len(nrow(edges)), function(i) {
        if (is_bidirectional[i]) {
          list(enabled = TRUE, type = "curvedCW", roundness = 0.1)
        } else {
          FALSE
        }
      }))
      
      # alternate curve direction for reverse edges
      bidirectional_indices <- which(is_bidirectional)
      for (i in bidirectional_indices) {
        reverse_idx <- which(edges$from == edges$to[i] & edges$to == edges$from[i])
        if (length(reverse_idx) > 0) {
          edges$smooth[[reverse_idx[1]]] <- list(enabled = TRUE, type = "curvedCCW", roundness = 0.1)
        }
      }
    } else {
      edges$smooth <- FALSE
    }
    
    # internal edges always light gray and thick
    internal_mask <- !is.na(edges$..is_internal..) & edges$..is_internal..
    edges$color.color[internal_mask] <- "#D3D3D3"
    edges$width[internal_mask] <- 18
    edges$width[!internal_mask] <- 6
    edges$label[internal_mask] <- ""
    
    # update edges via proxy
    visNetwork::visNetworkProxy("selectorGraph") %>%
      visNetwork::visSetOptions(options = list(
        edges = list(font = list(size = edge_font_size))
      )) %>%
      visNetwork::visUpdateEdges(edges)
  }
  
  # update nodes via proxy
  visNetwork::visNetworkProxy("selectorGraph") %>%
    visNetwork::visUpdateNodes(nodes)
})

# update library dropdowns when assembly changes
observeEvent(state$assembly, {
  if (is.null(state$assembly)) return()
  
  count_mat <- load_seg_adj_count(state$assembly)
  if (is.null(count_mat)) return()
  
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  lib_cols <- names(count_mat)[!names(count_mat) %in% metadata_cols]
  
  if (length(lib_cols) > 0) {
    lib_choices <- setNames(lib_cols, lib_cols)
    # restore cached library selections if they're still available
    cached_lib1 <- cache_get_if_exists("graph.edge_lib1", NULL)
    cached_lib2 <- cache_get_if_exists("graph.edge_lib2", NULL)
    selected_lib1 <- if (!is.null(cached_lib1) && cached_lib1 %in% lib_cols) cached_lib1 else NULL
    selected_lib2 <- if (!is.null(cached_lib2) && cached_lib2 %in% lib_cols) cached_lib2 else NULL
    updateSelectInput(session, "graphEdgeLib1", choices = lib_choices, selected = selected_lib1)
    updateSelectInput(session, "graphEdgeLib2", choices = lib_choices, selected = selected_lib2)
    
    # populate edge libs dropdown with "all" + each library
    edge_libs_choices <- c("all" = "all", lib_choices)
    cached_edge_libs <- cache_get_if_exists("graph.edge_libs", "all")
    selected_edge_libs <- if (cached_edge_libs %in% names(edge_libs_choices)) cached_edge_libs else "all"
    updateSelectInput(session, "graphEdgeLibs", choices = edge_libs_choices, selected = selected_edge_libs)
  }
}, ignoreNULL = FALSE)

# cache graph controls
observeEvent(input$graphEdgeMetric, {
  cache_set("graph.edge_metric", input$graphEdgeMetric)
})

observeEvent(input$graphEdgeLabels, {
  cache_set("graph.edge_labels", input$graphEdgeLabels)
})

observeEvent(input$graphNodeLabel, {
  cache_set("graph.node_label", input$graphNodeLabel)
})

observeEvent(input$graphFontSize, {
  cache_set("graph.font_size", input$graphFontSize)
})

observeEvent(input$graphDirectedEdges, {
  cache_set("graph.directed_edges", input$graphDirectedEdges)
})

observeEvent(input$graphEdgeLib1, {
  if (!is.null(input$graphEdgeLib1) && input$graphEdgeLib1 != "") {
    cache_set("graph.edge_lib1", input$graphEdgeLib1)
  }
})

observeEvent(input$graphEdgeLib2, {
  if (!is.null(input$graphEdgeLib2) && input$graphEdgeLib2 != "") {
    cache_set("graph.edge_lib2", input$graphEdgeLib2)
  }
})

observeEvent(input$graphEdgeLibs, {
  if (!is.null(input$graphEdgeLibs) && input$graphEdgeLibs != "") {
    cache_set("graph.edge_libs", input$graphEdgeLibs)
  }
})


observeEvent(input$selectorGraphNodeSelect, {
  nodes <- input$selectorGraphNodeSelect
  if (is.null(nodes)) nodes <- character()
  selected_graph_segments(nodes)
})

observeEvent(input$selectorGraphEdgeSelect, {
  edges <- input$selectorGraphEdgeSelect
  if (length(edges) > 0) {
    selected_graph_edge(edges[1])
  }
})

observeEvent(input$selectorGraphNodeHover, {
  hovered_node(input$selectorGraphNodeHover)
  # clear edge hover when hovering a node
  if (!is.null(input$selectorGraphNodeHover)) hovered_edge(NULL)
}, ignoreNULL = FALSE)

observeEvent(input$selectorGraphEdgeHover, {
  hovered_edge(input$selectorGraphEdgeHover)
  # clear node hover when hovering an edge
  if (!is.null(input$selectorGraphEdgeHover)) hovered_node(NULL)
}, ignoreNULL = FALSE)

#########################################################################
# graph zoom controls
#########################################################################

observeEvent(input$selectorGraphReset, {
  visNetwork::visNetworkProxy("selectorGraph") %>%
    visNetwork::visFit()
})

observeEvent(input$selectorGraphZoomIn, {
  session$sendCustomMessage("selectorGraphZoom", list(factor = 1.5))
})

observeEvent(input$selectorGraphZoomOut, {
  session$sendCustomMessage("selectorGraphZoom", list(factor = 0.67))
})

#########################################################################
# UI renderers
#########################################################################

output$selectorSelectedSegmentText <- renderText({
  node_ids <- selected_graph_segments()
  if (length(node_ids) == 0) {
    return("None")
  }
  seg_ids <- extract_segments_from_nodes(node_ids)
  return(paste(seg_ids, collapse = ", "))
})

output$selectorHoverInfo <- renderUI({
  node_hover <- hovered_node()
  edge_hover <- hovered_edge()
  
  # show hovered node details (side-nodes have segment field)
  if (!is.null(node_hover) && node_hover != "") {
    nodes <- base_graph_nodes()
    edges <- base_graph_edges()
    
    # extract segment from side-node
    seg_id <- node_hover
    side <- NULL
    if (!is.null(nodes) && nrow(nodes) > 0 && node_hover %in% nodes$id) {
      node_row <- nodes[nodes$id == node_hover, ]
      seg_id <- node_row$segment[1]
      side <- node_row$side[1]
    }
    
    seg_table <- get_segments(state$assembly)
    bin_segment_table <- load_bin_segment_table(state$assembly)
    
    # compute degree from edges (count edges connected to this node)
    degree <- 0
    if (!is.null(edges) && nrow(edges) > 0) {
      degree <- sum(edges$from == node_hover) + sum(edges$to == node_hover)
    }
    
    if (!is.null(seg_table) && seg_id %in% seg_table$segment) {
      seg_row <- seg_table[seg_table$segment == seg_id, ]
      contig <- seg_row$contig[1]
      start <- seg_row$start[1]
      end <- seg_row$end[1]
      length_kb <- round((end - start + 1) / 1000, 1)
      
      bin <- "unassigned"
      if (!is.null(bin_segment_table) && seg_id %in% bin_segment_table$segment) {
        bin <- bin_segment_table$bin[bin_segment_table$segment == seg_id][1]
      }
      
      side_text <- if (!is.null(side)) paste0(" (", side, ")") else ""
      return(tags$div(
        tags$strong("Segment: "), seg_id, side_text, tags$br(),
        tags$strong("Contig: "), contig, tags$br(),
        tags$strong("Coords: "), format(start, big.mark = ","), " - ", format(end, big.mark = ","), tags$br(),
        tags$strong("Length: "), length_kb, " kb", tags$br(),
        tags$strong("Bin: "), bin, tags$br(),
        tags$strong("Degree: "), degree
      ))
    }
    return(tags$p(tags$strong("Segment: "), seg_id, tags$br(), tags$strong("Degree: "), degree))
  }
  
  # show hovered edge details
  if (!is.null(edge_hover) && edge_hover != "") {
    edges <- base_graph_edges()
    
    if (!is.null(edges) && nrow(edges) > 0 && "id" %in% names(edges)) {
      edge_idx <- which(!is.na(edges$id) & edges$id == edge_hover)
      if (length(edge_idx) > 0) {
        edge_row <- edges[edge_idx[1], ]
        is_internal <- if (!is.null(edge_row$..is_internal..) && !is.na(edge_row$..is_internal..)) edge_row$..is_internal.. else FALSE
        
        # internal edge: show segment info
        if (is_internal) {
          seg_id <- edge_row$..segment..
          return(tags$div(
            tags$strong("Internal edge"), tags$br(),
            tags$strong("Segment: "), seg_id
          ))
        }
        
        # inter-segment edge
        seg_from <- edge_row$..seg_from..
        seg_to <- edge_row$..seg_to..
        side_from <- edge_row$..side_from..
        side_to <- edge_row$..side_to..
        support <- edge_row$..support..
        total <- edge_row$..total..
        total_associated <- edge_row$..total_associated..
        pct <- edge_row$..percent..
        is_directed <- if (!is.null(edge_row$arrows) && !is.na(edge_row$arrows) && edge_row$arrows == "to") TRUE else FALSE
        
        format_val <- function(x) if (!is.na(x)) format(x, big.mark = ",") else "N/A"
        format_pct <- function(x) if (!is.na(x)) paste0(round(x, 1), "%") else "N/A"
        
        edge_label <- if (is_directed) {
          paste0(seg_from, "_", side_from, " → ", seg_to, "_", side_to)
        } else {
          paste0(seg_from, "_", side_from, " - ", seg_to, "_", side_to)
        }
        
        # get library columns from edge data
        lib_cols <- character()
        pct_cols <- names(edge_row)[grepl("^\\.\\.pct_.*\\.\\.$", names(edge_row))]
        for (col in pct_cols) {
          lib_name <- gsub("^\\.\\.pct_(.*)\\.\\.$", "\\1", col)
          if (lib_name != "") lib_cols <- c(lib_cols, lib_name)
        }
        lib_cols <- unique(lib_cols)
        
        # get current metric to determine if we show change
        metric <- if (!is.null(input$graphEdgeMetric)) input$graphEdgeMetric else "none"
        lib1 <- input$graphEdgeLib1
        lib2 <- input$graphEdgeLib2
        
        # build table rows
        support_vals <- c()
        assoc_vals <- c()
        pct_vals <- c()
        
        for (lib in lib_cols) {
          sup_col <- paste0("..sup_", lib, "..")
          assoc_col <- paste0("..assoc_", lib, "..")
          pct_col <- paste0("..pct_", lib, "..")
          sup_val <- if (sup_col %in% names(edge_row)) edge_row[[sup_col]] else 0
          assoc_val <- if (assoc_col %in% names(edge_row)) edge_row[[assoc_col]] else 0
          pct_val <- if (pct_col %in% names(edge_row)) edge_row[[pct_col]] else 0
          sup_val[is.na(sup_val)] <- 0
          assoc_val[is.na(assoc_val)] <- 0
          pct_val[is.na(pct_val)] <- 0
          support_vals <- c(support_vals, round(sup_val))
          assoc_vals <- c(assoc_vals, round(assoc_val))
          pct_vals <- c(pct_vals, round(pct_val))
        }
        
        # add mean column
        support_vals <- c(support_vals, round(support))
        assoc_vals <- c(assoc_vals, round(total_associated))
        pct_vals <- c(pct_vals, round(pct))
        
        # build table
        header_cols <- c(lib_cols, "mean")
        table_rows <- list()
        
        # Support row
        support_cells <- lapply(support_vals, function(x) tags$td(format(x, big.mark = ","), style = "border: 1px solid #ccc; padding: 4px; text-align: right;"))
        table_rows <- c(table_rows, list(tags$tr(
          tags$td(tags$strong("Support"), style = "border: 1px solid #ccc; padding: 4px; text-align: left;"),
          support_cells
        )))
        
        # Total Associated row
        assoc_cells <- lapply(assoc_vals, function(x) tags$td(format(x, big.mark = ","), style = "border: 1px solid #ccc; padding: 4px; text-align: right;"))
        table_rows <- c(table_rows, list(tags$tr(
          tags$td(tags$strong("Total Associated"), style = "border: 1px solid #ccc; padding: 4px; text-align: left;"),
          assoc_cells
        )))
        
        # Percent row
        pct_cells <- lapply(pct_vals, function(x) tags$td(paste0(x, "%"), style = "border: 1px solid #ccc; padding: 4px; text-align: right;"))
        table_rows <- c(table_rows, list(tags$tr(
          tags$td(tags$strong("Percent"), style = "border: 1px solid #ccc; padding: 4px; text-align: left;"),
          pct_cells
        )))
        
        result_div <- tags$div(
          tags$strong("Edge: "), edge_label, tags$br(), tags$br(),
          tags$table(
            style = "border-collapse: collapse; font-size: 12px;",
            tags$thead(
              tags$tr(
                tags$th("", style = "border: 1px solid #ccc; padding: 4px; text-align: left;"),
                lapply(header_cols, function(h) tags$th(h, style = "border: 1px solid #ccc; padding: 4px; text-align: center;"))
              )
            ),
            tags$tbody(table_rows)
          )
        )
        
        return(result_div)
      }
    }
    return(tags$p(tags$strong("Edge: "), edge_hover))
  }
  
  # show selected segments if nothing hovered
  node_ids <- selected_graph_segments()
  if (length(node_ids) > 0) {
    seg_ids <- extract_segments_from_nodes(node_ids)
    return(tags$p(tags$strong("Selected: "), paste(seg_ids, collapse = ", ")))
  }
  
  return(tags$p("Hover over nodes or edges", style = "color: #999;"))
})

