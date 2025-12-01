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
  
  # expand neighbors iteratively
  all_segments <- base_segments
  current_frontier <- base_segments
  
  for (i in seq_len(depth)) {
    new_neighbors <- character()
    for (seg_id in current_frontier) {
      neighbors <- add_neighbors(seg_id)
      if (!is.null(neighbors)) {
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

#########################################################################
# graph action buttons (write directly to state$segments)
#########################################################################

# goto: navigate to selected segments in main view
observeEvent(input$selectorGotoBtn, {
  seg_ids <- selected_graph_segments()
  if (length(seg_ids) == 0) return()
  
  if (is.null(state$assembly)) return()
  seg_table <- get_segments(state$assembly)
  if (is.null(seg_table)) return()
  
  selected_segments <- seg_table[seg_table$segment %in% seg_ids, ]
  if (nrow(selected_segments) == 0) return()
  
  regions_module_output$push_undo_state()
  state$segments <- selected_segments
})

# add: append selected graph segments to state$segments
observeEvent(input$selectorAddBtn, {
  seg_ids <- selected_graph_segments()
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
    state$segments <- rbind(current_segs, new_segments)
  }
})

# remove: remove selected graph segments from state$segments
observeEvent(input$selectorRemoveBtn, {
  seg_ids <- selected_graph_segments()
  if (length(seg_ids) == 0) return()
  
  current_segs <- state$segments
  if (is.null(current_segs) || nrow(current_segs) == 0) return()
  
  new_segs <- current_segs[!current_segs$segment %in% seg_ids, , drop = FALSE]
  
  regions_module_output$push_undo_state()
  state$segments <- new_segs
})

# clear graph selection
observeEvent(input$selectorClearSelectionBtn, {
  selected_graph_segments(character())
})

#########################################################################
# graph rendering
#########################################################################

# store current graph data for proxy updates
current_graph_nodes <- reactiveVal(NULL)
current_graph_edges <- reactiveVal(NULL)
current_bin_segment_table <- reactiveVal(NULL)

output$selectorGraph <- visNetwork::renderVisNetwork({
  seg_ids <- computed_graph_segments()
  if (length(seg_ids) == 0) {
    current_graph_nodes(NULL)
    current_graph_edges(NULL)
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
  
  # color nodes by selected style
  color_by <- if (is.null(input$selectorColorBy)) "bin" else input$selectorColorBy
  if (color_by == "bin" && !is.null(bin_segment_table) && "bin_color" %in% names(bin_segment_table)) {
    seg_to_color <- setNames(bin_segment_table$bin_color, bin_segment_table$segment)
    nodes$color.background <- ifelse(nodes$id %in% names(seg_to_color), 
                                     seg_to_color[nodes$id], "#CCCCCC")
    nodes$color.border <- "#666666"
  } else {
    nodes$color.background <- "#DDE5FF"
    nodes$color.border <- "#4A90E2"
  }
  nodes$borderWidth <- 1
  
  if (nrow(edges) > 0) {
    edges$color.color <- "#848484"
  }
  
  current_graph_nodes(nodes)
  current_graph_edges(edges)
  current_bin_segment_table(bin_segment_table)
  
  network <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visEdges(arrows = "to") %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visPhysics(
      stabilization = list(enabled = TRUE, iterations = 200),
      solver = "forceAtlas2Based"
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
      Shiny.addCustomMessageHandler('selectorGraphZoom', function(data) {
        var scale = network.getScale() * data.factor;
        scale = Math.max(0.1, Math.min(scale, 10));
        network.moveTo({scale: scale});
      });
    }
  ")
  
  return(network)
})

# update node highlighting via proxy
observeEvent(selected_graph_segments(), {
  nodes <- current_graph_nodes()
  if (is.null(nodes) || nrow(nodes) == 0) return()
  
  selected_nodes <- selected_graph_segments()
  bin_segment_table <- current_bin_segment_table()
  
  color_by <- if (is.null(input$selectorColorBy)) "bin" else input$selectorColorBy
  if (color_by == "bin" && !is.null(bin_segment_table) && "bin_color" %in% names(bin_segment_table)) {
    seg_to_color <- setNames(bin_segment_table$bin_color, bin_segment_table$segment)
    nodes$color.background <- ifelse(nodes$id %in% names(seg_to_color), 
                                     seg_to_color[nodes$id], "#CCCCCC")
    nodes$color.border <- "#666666"
  } else {
    nodes$color.background <- "#DDE5FF"
    nodes$color.border <- "#4A90E2"
  }
  nodes$borderWidth <- 1
  
  if (length(selected_nodes) > 0) {
    nodes$color.border <- ifelse(nodes$id %in% selected_nodes, "#C92A2A", nodes$color.border)
    nodes$borderWidth <- ifelse(nodes$id %in% selected_nodes, 3, 1)
  }
  
  visNetwork::visNetworkProxy("selectorGraph") %>%
    visNetwork::visUpdateNodes(nodes)
}, ignoreNULL = FALSE)

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
})

observeEvent(input$selectorGraphEdgeHover, {
  hovered_edge(input$selectorGraphEdgeHover)
})

#########################################################################
# graph zoom controls
#########################################################################

observeEvent(input$selectorGraphReset, {
  visNetwork::visNetworkProxy("selectorGraph") %>%
    visNetwork::visFit()
})

observeEvent(input$selectorGraphZoomIn, {
  session$sendCustomMessage("selectorGraphZoom", list(factor = 2))
})

observeEvent(input$selectorGraphZoomOut, {
  session$sendCustomMessage("selectorGraphZoom", list(factor = 0.5))
})

#########################################################################
# UI renderers
#########################################################################

output$selectorSelectedSegmentText <- renderText({
  seg_ids <- selected_graph_segments()
  if (length(seg_ids) == 0) {
    return("None")
  }
  return(paste(seg_ids, collapse = ", "))
})

output$selectorHoverInfo <- renderUI({
  node_hover <- hovered_node()
  edge_hover <- hovered_edge()
  
  # show hovered node details
  if (!is.null(node_hover) && node_hover != "") {
    seg_id <- node_hover
    seg_table <- get_segments(state$assembly)
    bin_segment_table <- load_bin_segment_table(state$assembly)
    
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
      
      return(tags$div(
        tags$strong("Segment: "), seg_id, tags$br(),
        tags$strong("Contig: "), contig, tags$br(),
        tags$strong("Coords: "), format(start, big.mark = ","), " - ", format(end, big.mark = ","), tags$br(),
        tags$strong("Length: "), length_kb, " kb", tags$br(),
        tags$strong("Bin: "), bin
      ))
    }
    return(tags$p(tags$strong("Segment: "), seg_id))
  }
  
  # show hovered edge details
  if (!is.null(edge_hover) && edge_hover != "") {
    edges <- current_graph_edges()
    from_seg <- NULL
    to_seg <- NULL
    support <- NA
    total_associated <- NA
    total <- NA
    pct <- NA
    
    if (!is.null(edges) && nrow(edges) > 0 && "id" %in% names(edges)) {
      edge_idx <- which(edges$id == edge_hover)
      if (length(edge_idx) > 0) {
        from_seg <- edges$from[edge_idx[1]]
        to_seg <- edges$to[edge_idx[1]]
        if ("..support.." %in% names(edges)) support <- edges$..support..[edge_idx[1]]
        if ("..total_associated.." %in% names(edges)) total_associated <- edges$..total_associated..[edge_idx[1]]
        if ("..total.." %in% names(edges)) total <- edges$..total..[edge_idx[1]]
        if ("..percent.." %in% names(edges)) pct <- edges$..percent..[edge_idx[1]]
      }
    }
    
    if (is.null(from_seg)) {
      parts <- strsplit(edge_hover, "_")[[1]]
      if (length(parts) >= 2) {
        from_seg <- parts[1]
        to_seg <- parts[2]
      }
    }
    
    if (!is.null(from_seg) && !is.null(to_seg)) {
      return(tags$div(
        tags$strong("Edge: "), from_seg, " → ", to_seg, tags$br(),
        tags$strong("Supporting reads: "), if (!is.na(support)) format(support, big.mark = ",") else "N/A", tags$br(),
        tags$strong("Total associated: "), if (!is.na(total_associated)) format(total_associated, big.mark = ",") else "N/A", tags$br(),
        tags$strong("Total reads: "), if (!is.na(total)) format(total, big.mark = ",") else "N/A", tags$br(),
        tags$strong("Percent: "), if (!is.na(pct)) paste0(round(pct, 1), "%") else "N/A"
      ))
    }
    return(tags$p(tags$strong("Edge: "), edge_hover))
  }
  
  # show selected segments if nothing hovered
  seg_ids <- selected_graph_segments()
  if (length(seg_ids) > 0) {
    return(tags$p(tags$strong("Selected: "), paste(seg_ids, collapse = ", ")))
  }
  
  return(tags$p("Hover over nodes or edges", style = "color: #999;"))
})

