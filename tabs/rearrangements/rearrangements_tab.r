# rearrangement tab implementation

# load rearrangements utilities
source("tabs/rearrangements/rearrangements_utils.r", local = TRUE)
source("core/frequency_plots.r")

# validate and extract tab parameters
tab <- get_tab_by_id("rearrangements")
if (is.null(tab)) {
  stop("rearrangements tab not found during loading")
}

required_params <- c("get_aln_f", "library_ids")
missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
if (length(missing_params) > 0) {
  stop(sprintf("rearrangements tab missing required parameters: %s", paste(missing_params, collapse = ", ")))
}

get_aln_f <- tab$get_aln_f
library_ids <- tab$library_ids
max_gap <- tab$max_gap %||% 10
min_element_length <- tab$min_element_length %||% 50
min_anchor_length <- tab$min_anchor_length %||% 200
max_anchor_mutations_percent <- tab$max_anchor_mutations_percent %||% 0.01
max_element_mutation_percent <- tab$max_element_mutation_percent %||% 0.1

if (!is.function(get_aln_f)) {
  stop(sprintf("get_aln_f must be a function, got: %s", class(get_aln_f)))
}
if (!is.character(library_ids) || length(library_ids) == 0) {
  stop("library_ids must be a non-empty character vector")
}

# set the tab panel UI
set_tab_panel_f(function() {
  tabPanel(
    "Rearrangements",
    # top row: parameters and controls (25%) + plot (75%)
    fluidRow(
      column(3,
        wellPanel(
          h5("Rearrangement Parameters"),
          numericInput("maxGapInput", "Max Gap:", 
                      value = cache_get_if_exists("rearrange.max_gap", 10), 
                      min = 0, step = 1, width = "100%"),
          numericInput("minElementLengthInput", "Min Element Length:", 
                      value = cache_get_if_exists("rearrange.min_element_length", 50), 
                      min = 1, step = 1, width = "100%"),
          numericInput("minAnchorLengthInput", "Min Anchor Length:", 
                      value = cache_get_if_exists("rearrange.min_anchor_length", 200), 
                      min = 1, step = 1, width = "100%"),
          numericInput("maxAnchorMutationsPercentInput", "Max Anchor Mutations %:", 
                      value = cache_get_if_exists("rearrange.max_anchor_mutations_percent", 0.1), 
                      min = 0, max = 1, step = 0.001, width = "100%"),
          numericInput("maxElementMutationPercentInput", "Max Element Mutations %:", 
                      value = cache_get_if_exists("rearrange.max_element_mutation_percent", 1), 
                      min = 0, max = 1, step = 0.001, width = "100%"),
          br(),
          actionButton("updateRearrangementsBtn", "Update Rearrangements", class = "btn-primary", width = "100%"),
          br(), br(),
          verbatimTextOutput("rearrangementCountText", placeholder = TRUE),
          br(),
          checkboxInput("autoUpdateProfilesRearrangeChk", "Auto-update profiles", 
                       value = cache_get_if_exists("auto_update_profiles_rearrange", FALSE), width = "100%"),
          br(),
          h5("Filtering"),
          numericInput("rearrangementSpanFilter", "Min Span:", 
                      value = 0, min = 0, max = 1, step = 0.2, width = "100%"),
          numericInput("rearrangementSupportFilter", "Min Total Support:", 
                      value = 1, min = 1, step = 1, width = "100%")
        )
      ),
      column(9,
        create_frequency_plot_ui("rearrangement", library_ids)
      )
    ),
    # bottom row: rearrangements data with internal tabs (full width)
    fluidRow(
      column(12,
        h4("Rearrangements Data"),
        fluidRow(
          column(12,
            actionButton("gotoRearrangementsBtn", "Goto", class = "btn-secondary"),
            actionButton("clearRearrangementsBtn", "Clear Selection", class = "btn-secondary"),
            br(), br()
          )
        ),
        tabsetPanel(
          id = "rearrangementDataTabs",
          tabPanel("Events", 
                   DTOutput("eventsTable"),
                   value = "events_tab")
        )
      )
    )
  )
})

# ---- Rearrangement Data Management ----

# simplified query function that uses the utility function (like variants)
query_rearrangements <- function(assembly, contigs, zoom) {
  # create tab config object for the utility function
  tab_config <- list(
    get_aln_f = get_aln_f,
    library_ids = library_ids,
    max_gap = max_gap,
    min_element_length = min_element_length,
    min_anchor_length = min_anchor_length,
    max_anchor_mutations_percent = max_anchor_mutations_percent,
    max_element_mutation_percent = max_element_mutation_percent
  )
  
  return(query_rearrangements_for_context(assembly, contigs, zoom, tab_config))
}

# filter rearrangements by span (frequency range)
filter_rearrangements_by_span <- function(rearrange_data, min_span) {
  if (is.null(rearrange_data) || is.null(rearrange_data$support) || is.null(rearrange_data$coverage)) {
    return(rearrange_data)
  }
  
  # calculate frequency for each event across libraries
  support_matrix <- rearrange_data$support
  coverage_matrix <- rearrange_data$coverage
  
  # avoid division by zero
  freq_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
  
  # calculate span (max - min frequency) for each event
  event_spans <- apply(freq_matrix, 1, function(row) {
    valid_freqs <- row[!is.na(row) & is.finite(row)]
    if (length(valid_freqs) == 0) return(0)
    max(valid_freqs) - min(valid_freqs)
  })
  
  # filter events meeting span threshold
  keep_events <- event_spans >= min_span
  
  # filter all components
  filtered_data <- list(
    events = rearrange_data$events[keep_events, ],
    support = rearrange_data$support[keep_events, , drop = FALSE],
    coverage = rearrange_data$coverage[keep_events, , drop = FALSE],
    library_ids = rearrange_data$library_ids
  )
  
  return(filtered_data)
}

# filter rearrangements by total support
filter_rearrangements_by_support <- function(rearrange_data, min_support) {
  if (is.null(rearrange_data) || is.null(rearrange_data$events)) {
    return(rearrange_data)
  }
  
  events_df <- rearrange_data$events
  keep_events <- events_df$total_support >= min_support
  
  # filter all components
  filtered_data <- list(
    events = events_df[keep_events, , drop = FALSE],
    read_events = if (!is.null(rearrange_data$read_events)) {
      rearrange_data$read_events[rearrange_data$read_events$event_id %in% events_df$event_id[keep_events], , drop = FALSE]
    } else NULL,
    support = if (!is.null(rearrange_data$support)) rearrange_data$support[keep_events, , drop = FALSE] else NULL,
    coverage = if (!is.null(rearrange_data$coverage)) rearrange_data$coverage[keep_events, , drop = FALSE] else NULL,
    library_ids = rearrange_data$library_ids
  )
  
  return(filtered_data)
}

# helper function to format shim columns for display
format_shim_column <- function(shim_vector, max_display_length = 10) {
  sapply(shim_vector, function(shim) {
    if (is.na(shim) || nchar(shim) == 0) {
      return("")
    }
    if (nchar(shim) <= max_display_length) {
      return(shim)
    }
    # show first 5 + "..." + last 5
    paste0(substr(shim, 1, 5), "...", substr(shim, nchar(shim) - 4, nchar(shim)))
  })
}


# ---- Event Handlers ----

# shared function to process and store filtered rearrangements data
# handles filtering, sorting, coloring, and state updates
process_filtered_rearrangements <- function(raw_data, span_filter, support_filter, clear_selection = TRUE) {
  if (is.null(raw_data)) {
    state$filtered_rearrange_data <- NULL
    state$rearrangements <- NULL
    cache_set("rearrangements.current", NULL)
    return()
  }
  
  # clear selection if requested
  if (clear_selection) {
    selected_rearrangements(NULL)
    cache_set("rearrangements.selected", NULL)
  }
  
  # apply span filter
  filtered_data <- filter_rearrangements_by_span(raw_data, span_filter)
  
  # apply support filter
  filtered_data <- filter_rearrangements_by_support(filtered_data, support_filter)
  
  # sort the events by multiple columns for comprehensive ordering
  if (!is.null(filtered_data$events) && nrow(filtered_data$events) > 0) {
    # sort by: contig, out_clip, in_clip, element_contig, element_start, element_end
    events_df <- filtered_data$events
    sort_order <- order(
      events_df$contig,
      events_df$out_clip,
      events_df$in_clip,
      events_df$element_contig %||% "",  # handle missing element columns
      events_df$element_start %||% 0,
      events_df$element_end %||% 0
    )
    filtered_data$events <- filtered_data$events[sort_order, ]
    # also reorder the support and coverage matrices to match
    if (!is.null(filtered_data$support)) {
      filtered_data$support <- filtered_data$support[sort_order, , drop = FALSE]
    }
    if (!is.null(filtered_data$coverage)) {
      filtered_data$coverage <- filtered_data$coverage[sort_order, , drop = FALSE]
    }
  }
  
  # store filtered data
  state$filtered_rearrange_data <- filtered_data
  
  # add colors and store events dataframe in global state and cache for profiles
  if (!is.null(filtered_data$events)) {
    colored_events <- add_rearrangement_colors(filtered_data$events)
    state$rearrangements <- colored_events
    cache_set("rearrangements.current", colored_events)
  } else {
    state$rearrangements <- NULL
    cache_set("rearrangements.current", NULL)
  }
  
  # refresh profile plots only if auto-update is enabled
  if (input$autoUpdateProfilesRearrangeChk %||% FALSE) {
    if (exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
  } else {
    # mark plots as needing refresh when auto-update is disabled
    if (exists("invalidate_plot") && is.function(invalidate_plot)) {
      invalidate_plot()
    }
  }
}

# function to get selected events from DT table using actual displayed order
# provides a single source of truth for getting selected events from the table
get_selected_rearrangements <- function() {
  selected_rows <- input$eventsTable_rows_selected
  rearrange_data <- state$filtered_rearrange_data
  
  # return NULL if no selection or no data
  if (is.null(selected_rows) || length(selected_rows) == 0 || 
      is.null(rearrange_data) || is.null(rearrange_data$events) || 
      nrow(rearrange_data$events) == 0) {
    return(NULL)
  }
  
  # get the same data that was passed to DT
  display_df <- rearrange_data$events
  
  # validate row indices
  valid_rows <- selected_rows[selected_rows <= nrow(display_df)]
  if (length(valid_rows) == 0) {
    return(NULL)
  }
  
  # return the selected events from the displayed data
  return(display_df[valid_rows, ])
}

# reactive value to track selected rearrangements for highlighting
selected_rearrangements <- reactiveVal(NULL)

# observer for table row selection
observeEvent(input$eventsTable_rows_selected, {
  selected_events <- get_selected_rearrangements()
  
  # update reactive value and cache
  selected_rearrangements(selected_events)
  cache_set("rearrangements.selected", selected_events)
  
  # refresh plots to show/remove highlighting
  if (input$autoUpdateProfilesRearrangeChk %||% FALSE) {
    if (exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
  } else {
    # mark plots as needing refresh when highlighting changes but auto-update is disabled
    if (exists("invalidate_plot") && is.function(invalidate_plot)) {
      invalidate_plot()
    }
  }
}, ignoreNULL = FALSE)

# observer for auto-update profiles checkbox
observeEvent(input$autoUpdateProfilesRearrangeChk, {
  cache_set("auto_update_profiles_rearrange", input$autoUpdateProfilesRearrangeChk)
  
  # if auto-update is enabled and plots were out of sync, refresh now
  if (input$autoUpdateProfilesRearrangeChk && exists("plot_updated") && is.function(plot_updated)) {
    if (!plot_updated() && exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
  }
})

# goto button handler  
observeEvent(input$gotoRearrangementsBtn, {
  selected_events <- get_selected_rearrangements()
  
  if (is.null(selected_events)) {
    showNotification("Please select rearrangements to navigate to", type = "warning")
    return()
  }
  
  # build context to get global coordinates
  contigs_table <- get_contigs(state$assembly)
  if (is.null(contigs_table)) {
    showNotification("Cannot get contig information", type = "error")
    return()
  }
  
  cxt <- build_context(state$contigs, contigs_table, state$zoom, state$assembly)
  if (is.null(cxt)) {
    showNotification("Cannot build context for navigation", type = "error")
    return()
  }
  
  # convert to global coordinates using the context mapper
  # use out_clip as the coordinate for rearrangements
  selected_events$gcoord <- cxt$mapper$l2g(selected_events$contig, selected_events$out_clip)
  
  # calculate spanning range with appropriate margin
  min_coord <- min(selected_events$gcoord)
  max_coord <- max(selected_events$gcoord)
  
  if (nrow(selected_events) == 1) {
    # single event: minimum 10kb window
    window_size <- 10000  # 10kb minimum window
    center <- selected_events$gcoord[1]
    half_window <- window_size / 2
    zoom_start <- center - half_window
    zoom_end <- center + half_window
  } else {
    # multiple events: 10% margin on each side
    span <- max_coord - min_coord
    margin <- span * 0.1
    zoom_start <- min_coord - margin
    zoom_end <- max_coord + margin
  }
  
  # push current region to undo before changing
  regions_module_output$push_undo_state()
  
  # set zoom to calculated range
  state$zoom <- c(zoom_start, zoom_end)
  
  showNotification(sprintf("Navigated to %d selected rearrangements", nrow(selected_events)), type = "message")                                                                                                         
})

# clear selection button handler
observeEvent(input$clearRearrangementsBtn, {
  # clear the table selection
  DT::dataTableProxy("eventsTable") %>% 
    DT::selectRows(NULL)
  
  showNotification("Selection cleared", type = "message")
})

# update button handler
observeEvent(input$updateRearrangementsBtn, {
  # query rearrangements
  rearrange_data <- query_rearrangements(state$assembly, state$contigs, state$zoom)
  
  # store raw data
  state$raw_rearrange_data <- rearrange_data
  
  process_filtered_rearrangements(rearrange_data, input$rearrangementSpanFilter %||% 0.5, input$rearrangementSupportFilter %||% 1, clear_selection = TRUE)
})

# span filter change handler
observeEvent(input$rearrangementSpanFilter, {
  process_filtered_rearrangements(state$raw_rearrange_data, input$rearrangementSpanFilter, input$rearrangementSupportFilter %||% 1, clear_selection = TRUE)
})

# support filter change handler
observeEvent(input$rearrangementSupportFilter, {
  process_filtered_rearrangements(state$raw_rearrange_data, input$rearrangementSpanFilter %||% 0.5, input$rearrangementSupportFilter, clear_selection = TRUE)
})

# ---- Output Renderers ----

# rearrangement count text output
output$rearrangementCountText <- renderText({
  rearrange_data <- state$filtered_rearrange_data
  raw_data <- state$raw_rearrange_data
  
  if (is.null(raw_data) || is.null(raw_data$events)) {
    return("No data loaded")
  }
  
  total_events <- nrow(raw_data$events)
  
  if (is.null(rearrange_data) || is.null(rearrange_data$events)) {
    return(sprintf("Total: %d events", total_events))
  }
  
  filtered_events <- nrow(rearrange_data$events)
  
  return(sprintf("Events: %d/%d", filtered_events, total_events))
})

output$eventsTable <- renderDT({
  rearrange_data <- state$filtered_rearrange_data
  
  if (is.null(rearrange_data) || is.null(rearrange_data$events) || nrow(rearrange_data$events) == 0) {
    return(datatable(
      data.frame(Message = "Click 'Update' to load rearrangements for the current view"),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # format the events table for display
  display_df <- rearrange_data$events
  
  # round frequency to 3 decimal places
  if ("frequency" %in% names(display_df)) {
    display_df$frequency <- round(display_df$frequency, 3)
  }
  
  # format shim columns if they exist
  shim_cols <- c("left_shim", "right_shim", "middle_shim")
  for (col in shim_cols) {
    if (col %in% names(display_df)) {
      display_df[[paste0(col, "_display")]] <- format_shim_column(display_df[[col]])
    }
  }
  
  # complete column mapping for all fields
  col_names <- c(
    "Event ID" = "event_id",
    "Type" = "type",
    "Contig" = "contig",
    "Out Clip" = "out_clip",
    "In Clip" = "in_clip",
    "Element Contig" = "element_contig",
    "Element Strand" = "element_strand",
    "Element Start" = "element_start",
    "Element End" = "element_end",
    "Libraries" = "library_count",
    "Total Support" = "total_support",
    "Total Coverage" = "total_coverage",
    "Frequency" = "frequency",
    "Left Shim" = "left_shim_display",
    "Right Shim" = "right_shim_display",
    "Middle Shim" = "middle_shim_display"
  )
  
  # filter to available columns
  available_cols <- names(display_df)
  col_names <- col_names[col_names %in% available_cols]
  
  # create datatable with all available columns
  dt <- datatable(
    display_df[, col_names, drop = FALSE],
    rownames = FALSE,
    colnames = names(col_names),
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip",
      ordering = FALSE
    ),
    selection = list(mode = "multiple", target = "row"),
    filter = "top"
  )
  
  # add type column styling
  if ("type" %in% names(display_df)) {
    type_colors <- c("large_insert" = "red", "large_delete" = "blue", "large_invert" = "orange", "unknown" = "black")
    dt <- dt %>% formatStyle("type", color = styleEqual(names(type_colors), type_colors))
  }
  
  dt
})


# ---- Plot Renderer ----

# frequency plot output renderer
output$rearrangementFrequencyPlot <- plotly::renderPlotly({
  
  # get current data
  raw_data <- state$filtered_rearrange_data
  
  # get current settings with defaults
  plot_type <- input$rearrangementPlotType %||% "temporal"
  plot_value <- input$rearrangementPlotValue %||% "frequency"
  x_lib <- input$rearrangementXLib %||% library_ids[1]
  y_lib <- input$rearrangementYLib %||% (if(length(library_ids) > 1) library_ids[2] else library_ids[1])
  jitter_enabled <- input$rearrangementJitter %||% FALSE
  
  # get selected items
  selected_items <- get_selected_rearrangements()
  
  # check if we have data
  has_data <- !is.null(raw_data) && !is.null(raw_data$events)
  
  # prepare items dataframe if we have data
  items_df <- NULL
  if (has_data) {
    items_df <- add_rearrangement_colors(raw_data$events)
    items_df$id <- items_df$event_id
    items_df$label <- paste(items_df$type, items_df$contig, items_df$out_clip, sep = " ")
  }
  
  # render the plot using the cleaned internal function
  render_frequency_plot_internal(has_data, items_df, raw_data$support, raw_data$coverage,
                                plot_type, plot_value, x_lib, y_lib, 
                                jitter_enabled, selected_items, library_ids, 
                                "Click 'Update' to load rearrangements")
})

# frequency plot observers for caching
observeEvent(input$rearrangementPlotType, {
  cache_set("rearrangement_frequency_plot_type", input$rearrangementPlotType)
})

observeEvent(input$rearrangementPlotValue, {
  cache_set("rearrangement_frequency_plot_value", input$rearrangementPlotValue)
})

observeEvent(input$rearrangementXLib, {
  cache_set("rearrangement_frequency_plot_x_lib", input$rearrangementXLib)
})

observeEvent(input$rearrangementYLib, {
  cache_set("rearrangement_frequency_plot_y_lib", input$rearrangementYLib)
})

observeEvent(input$rearrangementJitter, {
  cache_set("rearrangement_frequency_plot_jitter", input$rearrangementJitter)
})

# export function for PDF generation
rearrangements_export_pdf <- function(region_info) {
  # query rearrangements for the specified region using tab config (like variants)
  raw_data <- query_rearrangements_for_context(region_info$assembly, region_info$contigs, region_info$context_zoom, get_tab_by_id("rearrangements"))
  
  # apply current filters
  span_filter <- input$rearrangementSpanFilter %||% cache_get_if_exists("rearrange.span_filter", 0.5)
  support_filter <- input$rearrangementSupportFilter %||% cache_get_if_exists("rearrange.support_filter", 1)
  
  # filter the data
  filtered_data <- filter_rearrangements_by_span(raw_data, span_filter)
  filtered_data <- filter_rearrangements_by_support(filtered_data, support_filter)
  
  # get current plot settings
  plot_type <- input$rearrangementPlotType %||% "temporal"
  plot_value <- input$rearrangementPlotValue %||% "frequency"
  x_lib <- input$rearrangementXLib %||% library_ids[1]
  y_lib <- input$rearrangementYLib %||% (if(length(library_ids) > 1) library_ids[2] else library_ids[1])
  jitter_enabled <- input$rearrangementJitter %||% FALSE
  
  # check if we have data
  has_data <- !is.null(filtered_data) && !is.null(filtered_data$events)
  
  # prepare items dataframe if we have data
  items_df <- NULL
  if (has_data) {
    items_df <- add_rearrangement_colors(filtered_data$events)
    items_df$id <- items_df$event_id
    items_df$label <- paste(items_df$type, items_df$contig, items_df$out_clip, sep = " ")
  }
  
  # create the plot using the unified export function
  return(create_frequency_plot_for_export(has_data, items_df, filtered_data$support, filtered_data$coverage,
                                         plot_type, plot_value, x_lib, y_lib, 
                                         jitter_enabled, library_ids, 
                                         title = "Rearrangements"))
}

# export function for table generation
rearrangements_export_table <- function(region_info) {
  # query rearrangements for the specified region using tab config (like variants)
  raw_data <- query_rearrangements_for_context(region_info$assembly, region_info$contigs, region_info$context_zoom, get_tab_by_id("rearrangements"))
  
  # apply current filters
  span_filter <- input$rearrangementSpanFilter %||% cache_get_if_exists("rearrange.span_filter", 0.5)
  support_filter <- input$rearrangementSupportFilter %||% cache_get_if_exists("rearrange.support_filter", 1)
  
  # filter the data
  filtered_data <- filter_rearrangements_by_span(raw_data, span_filter)
  filtered_data <- filter_rearrangements_by_support(filtered_data, support_filter)
  
  # return the filtered events dataframe (or NULL if no data)
  if (!is.null(filtered_data) && !is.null(filtered_data$events)) {
    return(filtered_data$events)
  }
  
  return(NULL)
}

# register the export functions for this tab
register_tab_export_function("rearrangements", "pdf", rearrangements_export_pdf)
register_tab_export_function("rearrangements", "table", rearrangements_export_table)
