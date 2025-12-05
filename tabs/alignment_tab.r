# ---- Alignment Tab Functions ----

# source alignment utilities for mutation coloring
source("profiles/align/align_utils.r")

# set the tab panel UI
set_tab_panel_f(function() {
  tabPanel(
    "Alignments",
    tags$style(HTML("
      #alignment_plot_container {
        border: 1px solid #ddd;
        border-radius: 4px;
        min-height: 80px;
      }
      #alignment_plot_container .ui-resizable-s {
        height: 8px !important;
        background: #e9ecef;
        border-top: 1px solid #d0d7de;
        cursor: ns-resize;
      }
    ")),
    fluidRow(
      column(12,
        div(
          style = "padding: 10px; margin-bottom: 5px;",
          fluidRow(
            column(3,
              selectInput("alignment_profile_select", "Profile:", choices = NULL, width = "100%")
            ),
            column(8,
              textInput("alignment_read_id_input", "Read ID:", value = "", width = "100%", placeholder = "Enter read ID")
            ),
            column(1,
              br(),
              actionButton("alignment_clear_btn", "Clear", style = "margin-top: 5px;")
            )
          )
        )
      )
    ),
    fluidRow(
      column(12,
        div(
          style = "background-color: #f5f5f5; padding: 10px; margin-bottom: 5px; border-radius: 4px; border: 1px solid #ddd;",
          uiOutput("alignmentReadInfo")
        )
      )
    ),
    fluidRow(
      column(12,
        div(
          style = "padding-left: 10px; margin-bottom: 10px;",
          checkboxInput("alignment_show_index", "show alignment index on plot", value = TRUE)
        )
      )
    ),
    fluidRow(
      column(12,
        shinyjqui::jqui_resizable(
          div(
            id = "alignment_plot_container",
            style = "height: 200px; padding: 5px;",
            plotly::plotlyOutput("alignmentPlot", height = "100%")
          ),
          options = list(handles = "s", minHeight = 80)
        )
      )
    ),
    fluidRow(
      column(12,
        style = "margin-top: 10px;",
        DTOutput("alignmentDetailsTable")
      )
    )
  )
})

# ---- Read Info Box Renderer ----

output$alignmentReadInfo <- renderUI({
  result <- state$clicked_read_data
  profile_id <- state$clicked_profile_id
  
  if (is.null(result) || is.null(result$alignments) || nrow(result$alignments) == 0) {
    return(div(
      style = "color: #666; font-style: italic;",
      "Click on a read in the main plot to view its alignments"
    ))
  }
  
  alignments <- result$alignments
  read_id <- alignments$read_id[1]
  read_length <- result$read_length
  num_alignments <- nrow(alignments)
  num_contigs <- length(unique(alignments$contig_id))
  
  div(
    style = "font-family: monospace;",
    tags$span(style = "font-weight: bold; font-size: 14px;", read_id),
    tags$span(style = "margin-left: 20px; color: #555;",
      sprintf("Length: %s bp", format(read_length, big.mark = ","))
    ),
    tags$span(style = "margin-left: 20px; color: #555;",
      sprintf("Alignments: %d", num_alignments)
    ),
    tags$span(style = "margin-left: 20px; color: #555;",
      sprintf("Contigs: %d", num_contigs)
    ),
    tags$span(style = "margin-left: 20px; color: #888;",
      sprintf("Profile: %s", profile_id)
    )
  )
})

# ---- Profile Dropdown Updater ----

# helper to get session object
get_session <- function() {
  if (exists("session", envir = parent.frame())) {
    return(get("session", envir = parent.frame()))
  }
  return(NULL)
}

# update profile dropdown choices based on available alignment objects
observe({
  aln_list <- cache_get_if_exists("aln_obj", list())
  session <- get_session()
  
  if (is.null(aln_list) || !is.list(aln_list) || length(aln_list) == 0) {
    if (!is.null(session)) {
      updateSelectInput(session, "alignment_profile_select", choices = list("(none)" = ""), selected = "")
    }
    return()
  }
  
  profile_ids <- names(aln_list)
  if (length(profile_ids) == 0) {
    if (!is.null(session)) {
      updateSelectInput(session, "alignment_profile_select", choices = list("(none)" = ""), selected = "")
    }
    return()
  }
  
  choices <- setNames(profile_ids, profile_ids)
  current_selection <- input$alignment_profile_select
  
  # if current selection is not in choices, use first available
  if (is.null(current_selection) || !is.element(current_selection, profile_ids)) {
    current_selection <- profile_ids[1]
  }
  
  if (!is.null(session)) {
    updateSelectInput(session, "alignment_profile_select", choices = choices, selected = current_selection)
  }
})

# ---- Helper Function to Load Read ----

# shared function to load read data from profile and read_id
load_read_from_profile <- function(profile_id, read_id) {
  if (is.null(profile_id) || profile_id == "" || is.null(read_id) || read_id == "") {
    state$clicked_read_data <- NULL
    state$clicked_profile_id <- NULL
    cache_set("selected_read_id", NULL)
    current_val <- refresh_trigger()
    refresh_trigger(current_val + 1)
    return()
  }
  
  # get alignment list and look up by profile id
  aln_list <- cache_get_if_exists("aln_obj", list())
  if (is.null(aln_list) || !is.list(aln_list)) {
    return()
  }
  
  aln <- aln_list[[profile_id]]
  if (is.null(aln)) {
    cat(sprintf("alignment not found for profile: %s\n", profile_id))
    return()
  }
  
  cat(sprintf("loading read from profile: %s, read_id: %s\n", profile_id, read_id))
  
  result <- tryCatch({
    aln_alignments_from_read_id(aln, read_id)
  }, error = function(e) {
    cat(sprintf("error getting alignments: %s\n", e$message))
    NULL
  })
  
  if (!is.null(result) && !is.null(result$alignments) && nrow(result$alignments) > 0) {
    state$clicked_read_data <- result
    state$clicked_profile_id <- profile_id
    cache_set("selected_read_id", read_id)
    # trigger profile refresh to show selected read highlighting
    current_val <- refresh_trigger()
    refresh_trigger(current_val + 1)
  } else {
    state$clicked_read_data <- NULL
    state$clicked_profile_id <- NULL
    cache_set("selected_read_id", NULL)
    # trigger profile refresh to clear highlighting
    current_val <- refresh_trigger()
    refresh_trigger(current_val + 1)
  }
}

# ---- Text Input Observer ----

# handle text input changes for read_id
observeEvent(input$alignment_read_id_input, {
  profile_id <- input$alignment_profile_select
  read_id <- input$alignment_read_id_input
  
  # only process if both profile and read_id are provided
  if (is.null(profile_id) || profile_id == "" || is.null(read_id) || read_id == "") {
    return()
  }
  
  load_read_from_profile(profile_id, read_id)
}, ignoreInit = TRUE)

# ---- Profile Selection Observer ----

# handle profile dropdown changes
observeEvent(input$alignment_profile_select, {
  profile_id <- input$alignment_profile_select
  read_id <- input$alignment_read_id_input
  
  # only process if both profile and read_id are provided
  if (is.null(profile_id) || profile_id == "" || is.null(read_id) || read_id == "") {
    return()
  }
  
  load_read_from_profile(profile_id, read_id)
}, ignoreInit = TRUE)

# ---- Clear Button Observer ----

# handle clear button click
observeEvent(input$alignment_clear_btn, {
  session <- get_session()
  if (!is.null(session)) {
    updateTextInput(session, "alignment_read_id_input", value = "")
  }
  state$clicked_read_data <- NULL
  state$clicked_profile_id <- NULL
  cache_set("selected_read_id", NULL)
  # trigger profile refresh to clear highlighting
  current_val <- refresh_trigger()
  refresh_trigger(current_val + 1)
})

# ---- Plotly Click Event Handler ----

# handle plotly click events from the combined plot to display alignment details
observeEvent(eventExpr = plotly::event_data("plotly_click"), {
  event_data <- plotly::event_data("plotly_click")
  req(event_data)
  req(event_data$key)
  
  # parse composite key: "profile_id:read_id"
  click_key <- event_data$key[[1]]
  key_parts <- strsplit(click_key, ":", fixed = TRUE)[[1]]
  if (length(key_parts) < 2) {
    cat(sprintf("invalid click key format: %s\n", click_key))
    return()
  }
  profile_id <- key_parts[1]
  read_id <- paste(key_parts[-1], collapse = ":")  # handle read_ids with colons
  
  # update UI inputs to reflect clicked read
  session <- get_session()
  if (!is.null(session)) {
    updateSelectInput(session, "alignment_profile_select", selected = profile_id)
    updateTextInput(session, "alignment_read_id_input", value = read_id)
  }
  
  # load the read data
  load_read_from_profile(profile_id, read_id)
})

# ---- Alignment Plot Renderer ----

output$alignmentPlot <- plotly::renderPlotly({
  result <- state$clicked_read_data
  show_index <- input$alignment_show_index
  
  if (is.null(result) || is.null(result$alignments) || nrow(result$alignments) == 0) {
    # empty plot with message
    gg <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                        label = "Click on a read in the main plot to view its alignments",
                        size = 5, color = "gray50") +
      ggplot2::theme_void() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
    return(plotly::ggplotly(gg, tooltip = "none"))
  }
  
  alignments <- result$alignments
  mutations <- result$mutations
  read_length <- result$read_length
  
  # assign heights based on contig ordering
  contig_order <- unique(alignments$contig_id[order(
    alignments$read_start,
    alignments$read_end,
    -alignments$num_mutations
  )])
  contig_heights <- setNames(seq_along(contig_order) - 1, contig_order)
  alignments$height <- contig_heights[alignments$contig_id]
  
  # sort alignments and add running index
  alignments <- alignments[order(alignments$height, alignments$read_start), ]
  alignments$index <- seq_len(nrow(alignments))
  
  # prepare alignment rectangles
  alignments$ymin <- alignments$height
  alignments$ymax <- alignments$height + 0.9
  
  # hover text
  alignments$hover_text <- paste0(
    "Contig: ", alignments$contig_id, "\n",
    "Read coords: ", alignments$read_start, "-", alignments$read_end, "\n",
    "Contig coords: ", alignments$contig_start, "-", alignments$contig_end, "\n",
    "Strand: ", ifelse(alignments$is_reverse, "reverse", "forward"), "\n",
    "Mutations: ", alignments$num_mutations
  )
  
  # create base plot
  gg <- ggplot2::ggplot()
  
  # draw strand indicator dots first (bottom layer)
  alignments$strand_x_filled <- ifelse(alignments$is_reverse, alignments$read_end, alignments$read_start)
  alignments$strand_x_open <- ifelse(alignments$is_reverse, alignments$read_start, alignments$read_end)
  alignments$strand_y <- (alignments$ymin + alignments$ymax) / 2
  
  # filled circle on one side
  gg <- gg + ggplot2::geom_point(
    data = alignments,
    ggplot2::aes(x = strand_x_filled, y = strand_y),
    color = "black",
    fill = "black",
    shape = 21,
    size = 2
  )
  
  # open circle on opposite side
  gg <- gg + ggplot2::geom_point(
    data = alignments,
    ggplot2::aes(x = strand_x_open, y = strand_y),
    color = "black",
    fill = "white",
    shape = 21,
    size = 2
  )
  
  # draw alignment rectangles
  gg <- gg + ggplot2::geom_rect(
    data = alignments,
    ggplot2::aes(
      xmin = read_start, xmax = read_end,
      ymin = ymin, ymax = ymax,
      text = hover_text
    ),
    fill = "lightgray",
    color = "gray50",
    size = 0.2
  )
  
  # plot mutations if available
  if (!is.null(mutations) && nrow(mutations) > 0) {
    # map alignment heights to mutations
    height_map <- setNames(alignments$height, alignments$alignment_index)
    mutations$height <- height_map[as.character(mutations$alignment_index)]
    
    mutations$ymin <- mutations$height
    mutations$ymax <- mutations$height + 0.9
    
    # get mutation colors using align_utils function
    mutations$fill_color <- get_variant_type_colors(mutations$desc)
    
    mutations$hover_text <- paste0(
      "Type: ", mutations$type, "\n",
      "Read coord: ", mutations$read_coord, "\n",
      "Contig coord: ", mutations$contig_coord, "\n",
      "Desc: ", mutations$desc
    )
    
    gg <- gg + ggplot2::geom_segment(
      data = mutations,
      ggplot2::aes(
        x = read_coord, xend = read_coord,
        y = ymin, yend = ymax,
        text = hover_text
      ),
      color = mutations$fill_color,
      size = 0.8
    )
  }
  
  # add index labels in the middle of each alignment rectangle (if enabled)
  if (!is.null(show_index) && show_index) {
    alignments$label_x <- (alignments$read_start + alignments$read_end) / 2
    alignments$label_y <- (alignments$ymin + alignments$ymax) / 2
    gg <- gg + ggplot2::geom_text(
      data = alignments,
      ggplot2::aes(x = label_x, y = label_y, label = index),
      color = "white",
      size = 3.5,
      fontface = "bold"
    )
  }
  
  # add read length indicator line at y = -0.2
  gg <- gg + ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = read_length, y = -0.2, yend = -0.2),
    color = "black", size = 0.5
  )
  
  # add baseline
  gg <- gg + ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
  
  # styling
  max_height <- max(alignments$height, na.rm = TRUE) + 1
  gg <- gg + 
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(color = "gray50", size = 0.3)
    ) +
    ggplot2::labs(x = "Read Position") +
    ggplot2::xlim(-read_length * 0.02, read_length * 1.02) +
    ggplot2::ylim(-0.5, max_height + 0.5)
  
  # convert to plotly with click events
  p <- plotly::ggplotly(gg, tooltip = "text", source = "alignmentPlotSource")
  p <- plotly::config(p, displayModeBar = FALSE)
  p
})

# ---- DataTable Renderer ----

output$alignmentDetailsTable <- renderDT({
  result <- state$clicked_read_data
  
  if (is.null(result) || is.null(result$alignments) || nrow(result$alignments) == 0) {
    return(datatable(
      data.frame(Message = "Click on a read in the main plot to see details here."), 
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  align_data <- result$alignments
  
  # assign heights based on contig ordering (same as in plot)
  contig_order <- unique(align_data$contig_id[order(
    align_data$read_start,
    align_data$read_end,
    -align_data$num_mutations
  )])
  contig_heights <- setNames(seq_along(contig_order) - 1, contig_order)
  align_data$height <- contig_heights[align_data$contig_id]
  
  # sort and add index
  align_data <- align_data[order(align_data$height, align_data$read_start), ]
  align_data$index <- seq_len(nrow(align_data))
  align_data$row_id <- align_data$index
  
  # add go to button column
  align_data$`Go to` <- sapply(align_data$row_id, function(id) {
    as.character(actionButton(
      paste0("goto_aln_btn_", id), "Go",
      onclick = sprintf("Shiny.setInputValue('goto_alignment_location_trigger', %d, {priority: 'event'})", id)
    ))
  })
  
  # columns to display
  display_cols <- c("index", "contig_id", "read_start", "read_end", "contig_start", "contig_end", 
                    "is_reverse", "num_mutations", "height", "Go to")
  
  datatable(
    align_data[, display_cols],
    escape = FALSE,
    rownames = FALSE,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      lengthMenu = c(10, 20, 50),
      dom = "lftip"
    ),
    selection = "none"
  )
})

# ---- Go To Button Observer ----

observeEvent(input$goto_alignment_location_trigger, {
  req(input$goto_alignment_location_trigger)
  row_idx <- as.integer(input$goto_alignment_location_trigger)
  
  result <- state$clicked_read_data
  req(result, result$alignments)
  
  # assign heights and sort the same way as the table
  align_data <- result$alignments
  contig_order <- unique(align_data$contig_id[order(
    align_data$read_start,
    align_data$read_end,
    -align_data$num_mutations
  )])
  contig_heights <- setNames(seq_along(contig_order) - 1, contig_order)
  align_data$height <- contig_heights[align_data$contig_id]
  align_data <- align_data[order(align_data$height, align_data$read_start), ]
  
  req(row_idx > 0, row_idx <= nrow(align_data))
  
  selected_alignment <- align_data[row_idx, ]
  req("contig_id" %in% colnames(selected_alignment))
  req("contig_start" %in% colnames(selected_alignment))
  req("contig_end" %in% colnames(selected_alignment))
  
  # push current region to undo before changing
  regions_module_output$push_undo_state()
  
  # get segments for selected contig
  segments <- get_segments(state$assembly)
  selected_segments <- segments[segments$contig == selected_alignment$contig_id, ]
  state$segments <- selected_segments
  gstart <- selected_alignment$contig_start
  gend <- selected_alignment$contig_end
  dd <- (gend - gstart) / 2
  state$zoom <- c(gstart - dd, gend + dd)
})
