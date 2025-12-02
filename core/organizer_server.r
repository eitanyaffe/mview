# organizer server
# segment list server logic

#########################################################################
# organizer JavaScript
#########################################################################

selector_organizer_js <- HTML("
window.selectorSelectedItems = new Set();

window.initSelectorSortable = function() {
  if (typeof Sortable === 'undefined') return;

  var el = document.getElementById('selectorSegmentList');
  if (!el) return;

  // click to toggle selection
  el.addEventListener('click', function(e) {
    var li = e.target.closest('li');
    if (!li) return;
    
    var id = li.getAttribute('data-id');
    if (e.altKey || e.metaKey) {
      // multi-select with Alt/Cmd
      if (window.selectorSelectedItems.has(id)) {
        window.selectorSelectedItems.delete(id);
        li.classList.remove('selected');
      } else {
        window.selectorSelectedItems.add(id);
        li.classList.add('selected');
      }
    } else {
      // single select
      el.querySelectorAll('li.selected').forEach(function(item) {
        item.classList.remove('selected');
      });
      window.selectorSelectedItems.clear();
      window.selectorSelectedItems.add(id);
      li.classList.add('selected');
    }
  });

  window.selectorSortable = new Sortable(el, {
    animation: 150,
    onEnd: function(evt) {
      var items = Array.from(el.children).map(function(x) { return x.getAttribute('data-id'); });
      if (typeof Shiny !== 'undefined') {
        Shiny.setInputValue('selectorSegmentOrder', items, {priority: 'event'});
      }
    }
  });
};

window.selectorFlipSegment = function(segId) {
  if (typeof Shiny !== 'undefined') {
    Shiny.setInputValue('selectorFlipSegment', {id: segId, time: Date.now()}, {priority: 'event'});
  }
};

window.selectorRemoveSegment = function(segId) {
  if (typeof Shiny !== 'undefined') {
    Shiny.setInputValue('selectorRemoveSegment', {id: segId, time: Date.now()}, {priority: 'event'});
  }
};

window.selectorRemoveSelected = function() {
  if (window.selectorSelectedItems.size === 0) return;
  
  var selectedIds = Array.from(window.selectorSelectedItems);
  window.selectorSelectedItems.clear();
  
  if (typeof Shiny !== 'undefined') {
    Shiny.setInputValue('selectorRemoveSegments', selectedIds, {priority: 'event'});
  }
};

window.selectorFlipSelected = function() {
  if (window.selectorSelectedItems.size === 0) return;
  
  var selectedIds = Array.from(window.selectorSelectedItems);
  
  if (typeof Shiny !== 'undefined') {
    Shiny.setInputValue('selectorFlipBtn', selectedIds, {priority: 'event'});
  }
};

Shiny.addCustomMessageHandler('selectorResetList', function(data) {
  var container = document.getElementById('selectorSegmentList');
  if (!container) return;
  
  var parent = container.parentNode;
  container.remove();
  
  var temp = document.createElement('div');
  temp.innerHTML = data.html;
  var newList = temp.firstChild;
  parent.appendChild(newList);
  
  window.initSelectorSortable();
});
")

#########################################################################
# helper functions
#########################################################################

# compute text color based on background luminance
get_text_color <- function(hex_color) {
  hex_color <- gsub("^#", "", hex_color)
  if (nchar(hex_color) != 6) return("white")
  
  r <- strtoi(substr(hex_color, 1, 2), base = 16)
  g <- strtoi(substr(hex_color, 3, 4), base = 16)
  b <- strtoi(substr(hex_color, 5, 6), base = 16)
  
  luminance <- (0.299 * r + 0.587 * g + 0.114 * b) / 255
  if (luminance > 0.5) "black" else "white"
}

# build segment list item
build_selector_item <- function(seg_id, strand, seg_table, bin_segment_table, seg_color = NULL) {
  seg_with_strand <- paste0(seg_id, strand)
  label_text <- seg_with_strand
  bg_color <- if (!is.null(seg_color) && !is.na(seg_color)) seg_color else "#888888"
  
  if (!is.null(seg_table) && seg_id %in% seg_table$segment) {
    seg_row <- seg_table[seg_table$segment == seg_id, ]
    contig <- seg_row$contig[1]
    length_bp <- seg_row$end[1] - seg_row$start[1] + 1
    length_kb <- round(length_bp / 1000, 1)
    
    bin <- NA
    if (!is.null(bin_segment_table) && seg_id %in% bin_segment_table$segment) {
      bin_row <- bin_segment_table[bin_segment_table$segment == seg_id, ]
      bin <- bin_row$bin[1]
    }
    
    label_parts <- c(seg_with_strand, paste0(length_kb, "kb"))
    if (!is.na(bin)) {
      label_parts <- c(label_parts, bin)
    }
    label_text <- paste(label_parts, collapse = " | ")
  }
  
  text_color <- get_text_color(bg_color)
  
  strand_label <- if (strand == "+") "F" else "R"
  tags$li(
    style = sprintf("background:%s; color:%s; display:flex; align-items:center; gap:8px;", bg_color, text_color),
    `data-id` = seg_id,
    `data-strand` = strand,
    tags$button(
      class = "flip-btn",
      onclick = sprintf("window.selectorFlipSegment('%s'); event.stopPropagation();", seg_id),
      style = "background:transparent; border:1px solid; border-radius:3px; padding:0 5px; cursor:pointer; font-size:14px; line-height:1.2; min-width:20px;",
      strand_label
    ),
    tags$span(label_text, style = "flex:1; margin-right:8px;"),
    tags$button(
      class = "remove-btn",
      onclick = sprintf("window.selectorRemoveSegment('%s'); event.stopPropagation();", seg_id),
      style = "background:transparent; border:none; padding:2px 4px; cursor:pointer; font-size:12px; opacity:0.6;",
      "\u2715"
    )
  )
}

#########################################################################
# sync sequence from state$segments
#########################################################################

observe({
  state_segs <- state$segments
  if (is.null(state_segs) || nrow(state_segs) == 0) {
    set_sequence_segments(data.frame(segment = character(), strand = character(), stringsAsFactors = FALSE), "sync_from_state_empty")
  } else {
    strand_col <- if ("strand" %in% names(state_segs)) state_segs$strand else rep("+", nrow(state_segs))
    set_sequence_segments(data.frame(segment = state_segs$segment, strand = strand_col, stringsAsFactors = FALSE), "sync_from_state")
  }
})

#########################################################################
# pending changes detection
#########################################################################

has_pending_changes <- reactive({
  seq_df <- get_sequence_segments()
  seq_ids <- if (nrow(seq_df) == 0) character() else seq_df$segment
  seq_strands <- if (nrow(seq_df) == 0) character() else seq_df$strand
  
  state_segs <- state$segments
  state_ids <- if (is.null(state_segs) || nrow(state_segs) == 0) character() else state_segs$segment
  state_strands <- if (is.null(state_segs) || nrow(state_segs) == 0) character() else {
    if ("strand" %in% names(state_segs)) state_segs$strand else rep("+", nrow(state_segs))
  }
  
  !identical(seq_ids, state_ids) || !identical(seq_strands, state_strands)
})

#########################################################################
# organizer button handlers
#########################################################################

# apply button - update state$segments with current order and strand
observeEvent(input$selectorApplyBtn, {
  seq_df <- get_sequence_segments()
  if (nrow(seq_df) == 0) {
    regions_module_output$push_undo_state()
    state$segments <- state$segments[0, ]
    return()
  }
  
  if (is.null(state$assembly)) return()
  
  seg_table <- get_segments(state$assembly)
  if (is.null(seg_table)) return()
  
  regions_module_output$push_undo_state()
  
  selected_segments <- seg_table[seg_table$segment %in% seq_df$segment, ]
  if (nrow(selected_segments) == 0) return()
  
  seg_order <- match(selected_segments$segment, seq_df$segment)
  selected_segments <- selected_segments[order(seg_order), ]
  
  strand_lookup <- setNames(seq_df$strand, seq_df$segment)
  selected_segments$strand <- strand_lookup[selected_segments$segment]
  
  state$segments <- selected_segments
})

# flip segment strand (immediate, with undo)
observeEvent(input$selectorFlipSegment, {
  flip_data <- input$selectorFlipSegment
  if (is.null(flip_data) || is.null(flip_data$id)) return()
  
  seg_id <- flip_data$id
  current_segs <- state$segments
  if (is.null(current_segs) || nrow(current_segs) == 0) return()
  
  idx <- which(current_segs$segment == seg_id)
  if (length(idx) == 0) return()
  
  if (!"strand" %in% names(current_segs)) {
    current_segs$strand <- "+"
  }
  
  regions_module_output$push_undo_state()
  current_segs$strand[idx] <- if (current_segs$strand[idx] == "+") "-" else "+"
  state$segments <- current_segs
})

# remove single segment (immediate, with undo)
observeEvent(input$selectorRemoveSegment, {
  remove_data <- input$selectorRemoveSegment
  if (is.null(remove_data) || is.null(remove_data$id)) return()
  
  seg_id <- remove_data$id
  current_segs <- state$segments
  if (is.null(current_segs) || nrow(current_segs) == 0) return()
  
  idx <- which(current_segs$segment == seg_id)
  if (length(idx) == 0) return()
  
  regions_module_output$push_undo_state()
  state$segments <- current_segs[-idx, , drop = FALSE]
})

# reset button - revert to state$segments
observeEvent(input$selectorResetBtn, {
  state_segs <- state$segments
  if (is.null(state_segs) || nrow(state_segs) == 0) {
    set_sequence_segments(data.frame(segment = character(), strand = character(), stringsAsFactors = FALSE), "reset_empty")
  } else {
    strand_col <- if ("strand" %in% names(state_segs)) state_segs$strand else rep("+", nrow(state_segs))
    set_sequence_segments(data.frame(segment = state_segs$segment, strand = strand_col, stringsAsFactors = FALSE), "reset")
  }
  
  seg_table <- get_segments(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  seq_df <- get_sequence_segments()
  
  # get colors for all segments
  seg_colors <- get_segment_colors(seq_df$segment, state$assembly, seq_df$segment)
  
  list_items <- lapply(seq_len(nrow(seq_df)), function(i) {
    seg_id <- seq_df$segment[i]
    seg_color <- if (seg_id %in% names(seg_colors)) seg_colors[seg_id] else NULL
    build_selector_item(seg_id, seq_df$strand[i], seg_table, bin_segment_table, seg_color)
  })
  
  session$sendCustomMessage("selectorResetList", list(
    html = as.character(tags$ul(id = "selectorSegmentList", list_items))
  ))
})

# track order changes from drag
observeEvent(input$selectorSegmentOrder, {
  new_order <- input$selectorSegmentOrder
  current_df <- get_sequence_segments()
  
  if (nrow(current_df) > 0 && length(new_order) > 0) {
    strand_lookup <- setNames(current_df$strand, current_df$segment)
    new_strands <- strand_lookup[new_order]
    new_strands[is.na(new_strands)] <- "+"
    set_sequence_segments(data.frame(segment = new_order, strand = new_strands, stringsAsFactors = FALSE), "reorder")
  } else {
    set_sequence_segments(new_order, "reorder")
  }
})

# remove selected segments
observeEvent(input$selectorRemoveSegments, {
  ids_to_remove <- input$selectorRemoveSegments
  if (length(ids_to_remove) == 0) return()
  
  current_df <- get_sequence_segments()
  new_df <- current_df[!current_df$segment %in% ids_to_remove, , drop = FALSE]
  set_sequence_segments(new_df, "remove_selected")
  
  seg_table <- get_segments(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  
  # get colors for all segments
  seg_colors <- get_segment_colors(new_df$segment, state$assembly, new_df$segment)
  
  list_items <- lapply(seq_len(nrow(new_df)), function(i) {
    seg_id <- new_df$segment[i]
    seg_color <- if (seg_id %in% names(seg_colors)) seg_colors[seg_id] else NULL
    build_selector_item(seg_id, new_df$strand[i], seg_table, bin_segment_table, seg_color)
  })
  
  session$sendCustomMessage("selectorResetList", list(
    html = as.character(tags$ul(id = "selectorSegmentList", list_items))
  ))
})

# flip strand for selected segments
observeEvent(input$selectorFlipBtn, {
  flip_segment_strands(input$selectorFlipBtn)
  
  seg_table <- get_segments(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  seq_df <- get_sequence_segments()
  
  # get colors for all segments
  seg_colors <- get_segment_colors(seq_df$segment, state$assembly, seq_df$segment)
  
  list_items <- lapply(seq_len(nrow(seq_df)), function(i) {
    seg_id <- seq_df$segment[i]
    seg_color <- if (seg_id %in% names(seg_colors)) seg_colors[seg_id] else NULL
    build_selector_item(seg_id, seq_df$strand[i], seg_table, bin_segment_table, seg_color)
  })
  
  session$sendCustomMessage("selectorResetList", list(
    html = as.character(tags$ul(id = "selectorSegmentList", list_items))
  ))
})

#########################################################################
# UI renderers
#########################################################################

output$selectorOrganizerButtons <- renderUI({
  seq_df <- get_sequence_segments()
  has_segments <- nrow(seq_df) > 0
  
  remove_class <- if (has_segments) "" else "faded"
  
  tags$div(
    class = "selector-buttons",
    tags$button("Remove", type = "button", 
                class = paste("btn btn-default btn-sm", remove_class),
                onclick = "window.selectorRemoveSelected();"),
    tags$button("Clear", type = "button", 
                class = paste("btn btn-default btn-sm", remove_class),
                onclick = "window.selectorClearSelection();")
  )
})

output$selectorSegmentListUI <- renderUI({
  seq_df <- get_sequence_segments()
  # trigger re-render when color scheme changes
  input$segmentColorScheme
  
  if (nrow(seq_df) == 0) {
    return(tags$div(
      tags$script(selector_organizer_js),
      tags$p("No segments. Add from graph using buttons.", style = "color: #999;")
    ))
  }
  
  seg_table <- get_segments(state$assembly)
  bin_segment_table <- load_bin_segment_table(state$assembly)
  
  # get colors for all segments
  seg_colors <- get_segment_colors(seq_df$segment, state$assembly, seq_df$segment)
  
  list_items <- lapply(seq_len(nrow(seq_df)), function(i) {
    seg_id <- seq_df$segment[i]
    seg_color <- if (seg_id %in% names(seg_colors)) seg_colors[seg_id] else NULL
    build_selector_item(seg_id, seq_df$strand[i], seg_table, bin_segment_table, seg_color)
  })
  
  tags$div(
    tags$script(selector_organizer_js),
    tags$ul(id = "selectorSegmentList", list_items),
    tags$script(HTML("
      setTimeout(function() {
        if (typeof Sortable !== 'undefined') {
          window.initSelectorSortable();
        }
      }, 100);
    "))
  )
})

