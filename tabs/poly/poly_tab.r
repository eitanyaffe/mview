# poly tab implementation

# load poly utilities
source("tabs/poly/poly_utils.r", local = TRUE)

# validate and extract tab parameters
tab <- get_tab_by_id("poly")
if (is.null(tab)) {
  stop("poly tab not found during loading")
}

# extract required parameters
get_unify_table_f <- tab$get_unify_table_f
get_unify_support_f <- tab$get_unify_support_f
get_unify_coverage_f <- tab$get_unify_coverage_f
get_abundance_f <- tab$get_abundance_f
library_ids <- tab$library_ids

# validate functions
required_params <- c("get_unify_table_f", "get_unify_support_f", "get_unify_coverage_f", 
                    "get_abundance_f", "library_ids")
missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
if (length(missing_params) > 0) {
  stop(sprintf("poly tab missing required parameters: %s", paste(missing_params, collapse = ", ")))
}

if (!is.character(library_ids) || length(library_ids) == 0) {
  stop("library_ids must be a non-empty character vector")
}

if (!is.function(get_unify_table_f)) {
  stop(sprintf("get_unify_table_f must be a function, got: %s", class(get_unify_table_f)))
}
if (!is.function(get_unify_support_f)) {
  stop(sprintf("get_unify_support_f must be a function, got: %s", class(get_unify_support_f)))
}
if (!is.function(get_unify_coverage_f)) {
  stop(sprintf("get_unify_coverage_f must be a function, got: %s", class(get_unify_coverage_f)))
}
if (!is.function(get_abundance_f)) {
  stop(sprintf("get_abundance_f must be a function, got: %s", class(get_abundance_f)))
}

# reactive values for state management
poly_auto_update <- reactiveVal(cache_get_if_exists("poly.auto_update", "view"))
poly_host_bin <- reactiveVal(cache_get_if_exists("poly.host_bin", ""))
poly_min_span <- reactiveVal(cache_get_if_exists("poly.min_span", 0))

# flag to track if user has manually changed the host bin
# if TRUE, don't auto-update from last_selected_genome
poly_host_manually_set <- reactiveVal(FALSE)
poly_color_mode <- reactiveVal(cache_get_if_exists("poly.color_mode", "rainbow"))
poly_frozen_set <- reactiveVal(NULL)
poly_plot_type <- reactiveVal(cache_get_if_exists("poly.plot_type", "temporal"))
poly_x_lib <- reactiveVal(cache_get_if_exists("poly.x_lib", library_ids[1]))
poly_y_lib <- reactiveVal(cache_get_if_exists("poly.y_lib", if(length(library_ids) > 1) library_ids[2] else library_ids[1]))
poly_jitter <- reactiveVal(cache_get_if_exists("poly.jitter", FALSE))
poly_hover <- reactiveVal(cache_get_if_exists("poly.hover", TRUE))
poly_selected_uids <- reactiveVal(NULL)
poly_show_locals <- reactiveVal(cache_get_if_exists("poly.show_locals", TRUE))
poly_show_rearrangements <- reactiveVal(cache_get_if_exists("poly.show_rearrangements", TRUE))
poly_show_elements <- reactiveVal(cache_get_if_exists("poly.show_elements", TRUE))

# set the tab panel UI
set_tab_panel_f(function() {
  tabPanel(
    "Poly",
    # top controls panel
    fluidRow(
      column(12,
        wellPanel(
          h5("Controls"),
          fluidRow(
            column(2,
              selectInput("polyAutoUpdateInput", "Auto-update:",
                         choices = list("View" = "view", "Zoom" = "zoom", "Host" = "host"),
                         selected = cache_get_if_exists("poly.auto_update", "view"),
                         width = "100%")
            ),
            column(3,
              textInput("polyHostInput", "Host:",
                       value = cache_get_if_exists("poly.host_bin", ""),
                       placeholder = "Enter host bin ID",
                       width = "100%")
            ),
            column(2,
              numericInput("polyMinSpanInput", "Min Span:",
                          value = cache_get_if_exists("poly.min_span", 0),
                          min = 0, max = 1, step = 0.1, width = "100%")
            ),
            column(2,
              selectInput("polyColorModeInput", "Color Mode:",
                         choices = list("Rainbow" = "rainbow", "Sweep" = "sweep", "Type" = "type"),
                         selected = cache_get_if_exists("poly.color_mode", "rainbow"),
                         width = "100%")
            ),
            column(3,
              h6("Show Variants:"),
              checkboxInput("polyShowLocalsInput", "Locals",
                          value = cache_get_if_exists("poly.show_locals", TRUE)),
              checkboxInput("polyShowRearrangementsInput", "Rearrangements",
                          value = cache_get_if_exists("poly.show_rearrangements", TRUE)),
              checkboxInput("polyShowElementsInput", "Elements",
                          value = cache_get_if_exists("poly.show_elements", TRUE))
            )
          )
        )
      )
    ),
    # plot controls row
    fluidRow(
      column(12,
        wellPanel(
          h5("Plot Controls"),
          fluidRow(
            column(2,
              selectInput("polyPlotTypeInput", "Plot Type:",
                         choices = list("Temporal" = "temporal", "Scatter" = "scatter"),
                         selected = cache_get_if_exists("poly.plot_type", "temporal"),
                         width = "100%")
            ),
            column(2,
              conditionalPanel(
                condition = "input.polyPlotTypeInput == 'scatter'",
                selectInput("polyXLibInput", "X Library:",
                           choices = setNames(library_ids, library_ids),
                           selected = cache_get_if_exists("poly.x_lib", library_ids[1]),
                           width = "100%")
              )
            ),
            column(2,
              conditionalPanel(
                condition = "input.polyPlotTypeInput == 'scatter'",
                selectInput("polyYLibInput", "Y Library:",
                           choices = setNames(library_ids, library_ids),
                           selected = cache_get_if_exists("poly.y_lib", if(length(library_ids) > 1) library_ids[2] else library_ids[1]),
                           width = "100%")
              )
            ),
            column(2,
              checkboxInput("polyJitterInput", "Jitter",
                          value = cache_get_if_exists("poly.jitter", FALSE),
                          width = "100%"),
              checkboxInput("polyHoverInput", "Hover",
                          value = cache_get_if_exists("poly.hover", TRUE),
                          width = "100%")
            )
          )
        )
      )
    ),
    # temporal plots row
    fluidRow(
      column(3,
        conditionalPanel(
          condition = "input.polyAutoUpdateInput == 'host' && input.polyHostInput != ''",
          h4("Host Abundance"),
          plotly::plotlyOutput("polyHostAbundancePlot", height = "300px")
        )
      ),
      column(3,
        h4("Support"),
        plotly::plotlyOutput("polySupportPlot", height = "300px")
      ),
      column(3,
        h4("Coverage"),
        plotly::plotlyOutput("polyCoveragePlot", height = "300px")
      ),
      column(3,
        h4("Percent"),
        plotly::plotlyOutput("polyPercentPlot", height = "300px")
      )
    ),
    # tables row
    fluidRow(
      column(12,
        h4("Variant Tables"),
        fluidRow(
          column(12,
            actionButton("polyGotoSiteBtn", "Goto Site", class = "btn-secondary"),
            actionButton("polyGotoElementBtn", "Goto Element", class = "btn-secondary"),
            actionButton("polyClearSelectionBtn", "Clear Selection", class = "btn-secondary"),
            br(), br()
          )
        ),
        DTOutput("polyVariantsTable")
      )
    )
  )
})

# ---- Data Loading Functions ----

# load unified data for current assembly
load_unified_data_for_assembly <- function(assembly) {
  if (is.null(assembly)) {
    return(NULL)
  }
  
  tab_config <- list(
    get_unify_table_f = get_unify_table_f,
    get_unify_support_f = get_unify_support_f,
    get_unify_coverage_f = get_unify_coverage_f,
    library_ids = library_ids
  )
  
  return(load_unified_data(assembly, tab_config))
}

# load host abundance data
load_host_abundance_data <- function(assembly, host_bin) {
  if (is.null(assembly) || is.null(host_bin) || host_bin == "") {
    return(NULL)
  }
  
  abundance_mat <- get_abundance_f(assembly)
  if (is.null(abundance_mat)) {
    return(NULL)
  }
  
  # filter by bin
  if (!"bin" %in% names(abundance_mat)) {
    return(NULL)
  }
  
  bin_row <- abundance_mat[abundance_mat$bin == host_bin, ]
  if (nrow(bin_row) == 0) {
    return(NULL)
  }
  
  # get library columns (all columns except bin)
  lib_cols <- names(bin_row)[names(bin_row) != "bin"]
  
  # match library columns to library_ids (try both with and without assembly prefix)
  values <- numeric(length(library_ids))
  lib_ids_used <- character(length(library_ids))
  
  for (i in seq_along(library_ids)) {
    lib_id <- library_ids[i]
    lib_ids_used[i] <- lib_id
    
    # try exact match first
    if (lib_id %in% lib_cols) {
      values[i] <- as.numeric(bin_row[[lib_id]])
    } else {
      # try with assembly prefix
      prefixed_lib_id <- paste0(assembly, "_", lib_id)
      if (prefixed_lib_id %in% lib_cols) {
        values[i] <- as.numeric(bin_row[[prefixed_lib_id]])
      } else {
        values[i] <- 0
      }
    }
  }
  
  return(list(libraries = lib_ids_used, values = values))
}

# process and filter unified data based on current settings
process_filtered_poly_data <- function(raw_data, auto_update, host_bin, min_span, contigs, zoom, assembly, clear_selection = TRUE) {
  if (is.null(raw_data)) {
    state$filtered_poly_data <- NULL
    state$poly_split_data <- NULL
    return()
  }
  
  # determine filter mode based on auto_update setting
  if (auto_update == "host" && !is.null(host_bin) && host_bin != "") {
    # host mode with host: filter by host
    filter_mode <- "host"
  } else if (auto_update == "view") {
    # view mode: filter by contigs
    filter_mode <- "view"
  } else if (auto_update == "zoom") {
    # zoom mode: filter by zoom
    filter_mode <- "zoom"
  } else {
    # host mode without host: show all
    filter_mode <- "view"
    contigs <- NULL  # show all contigs
  }
  
  # filter by mode
  filtered_data <- filter_by_mode(raw_data, filter_mode, host_bin, contigs, zoom, assembly, NULL)
  
  if (is.null(filtered_data)) {
    state$filtered_poly_data <- NULL
    state$poly_split_data <- NULL
    cat(sprintf("poly: mode=%s, bin=%s, variants=0\n", auto_update, if(is.null(host_bin) || host_bin == "") "none" else host_bin))
    return()
  }
  
  # filter by span
  filtered_data <- filter_by_span(filtered_data, min_span)
  
  if (is.null(filtered_data)) {
    state$filtered_poly_data <- NULL
    state$poly_split_data <- NULL
    cat(sprintf("poly: mode=%s, bin=%s, variants=0\n", auto_update, if(is.null(host_bin) || host_bin == "") "none" else host_bin))
    return()
  }
  
  # filter by selected types before splitting
  if (!is.null(filtered_data) && !is.null(filtered_data$table)) {
    type_filter <- rep(TRUE, nrow(filtered_data$table))
    
    if (!poly_show_locals()) {
      type_filter <- type_filter & (filtered_data$table$type != "local")
    }
    if (!poly_show_rearrangements()) {
      type_filter <- type_filter & (filtered_data$table$type != "rearrange")
    }
    if (!poly_show_elements()) {
      type_filter <- type_filter & (filtered_data$table$type != "element")
    }
    
    # apply type filter
    if (any(type_filter)) {
      filtered_data$table <- filtered_data$table[type_filter, , drop = FALSE]
      if (!is.null(filtered_data$support)) {
        filtered_data$support <- filtered_data$support[type_filter, , drop = FALSE]
      }
      if (!is.null(filtered_data$coverage)) {
        filtered_data$coverage <- filtered_data$coverage[type_filter, , drop = FALSE]
      }
    } else {
      # no types selected, set to NULL
      filtered_data <- NULL
    }
  }
  
  # store filtered data
  state$filtered_poly_data <- filtered_data
  
  # split by type (for tables)
  if (!is.null(filtered_data)) {
    split_data <- split_by_type(filtered_data)
  } else {
    split_data <- list(locals = NULL, rearrangements = NULL, elements = NULL)
  }
  
  # filter split data by selected types (for tables)
  if (!poly_show_locals()) {
    split_data$locals <- NULL
  }
  if (!poly_show_rearrangements()) {
    split_data$rearrangements <- NULL
  }
  if (!poly_show_elements()) {
    split_data$elements <- NULL
  }
  
  state$poly_split_data <- split_data
  
  # print summary: mode, bin, number of variants
  num_variants <- if (!is.null(filtered_data) && !is.null(filtered_data$table)) nrow(filtered_data$table) else 0
  cat(sprintf("poly: mode=%s, bin=%s, variants=%d\n", auto_update, if(is.null(host_bin) || host_bin == "") "none" else host_bin, num_variants))
}

# ---- Event Handlers ----

# observer for auto-update change
observeEvent(input$polyAutoUpdateInput, {
  poly_auto_update(input$polyAutoUpdateInput)
  cache_set("poly.auto_update", input$polyAutoUpdateInput)
  
  # reload and filter data (only if not in host mode)
  if (input$polyAutoUpdateInput != "host") {
    raw_data <- state$raw_poly_data
    
    # if raw_data is null, try to load it
    if (is.null(raw_data) && !is.null(state$assembly)) {
      raw_data <- load_unified_data_for_assembly(state$assembly)
      state$raw_poly_data <- raw_data
    }
    
    if (!is.null(raw_data)) {
      process_filtered_poly_data(raw_data, input$polyAutoUpdateInput, 
                                poly_host_bin(), poly_min_span(),
                                get_state_contigs(), state$zoom, state$assembly)
    }
  }
})

# observer for host input change (filter automatically when in host mode)
observeEvent(input$polyHostInput, {
  host_bin <- trimws(input$polyHostInput)
  poly_host_bin(host_bin)
  cache_set("poly.host_bin", host_bin)
  
  # check if user typed something different from selected genome
  current_genome <- last_selected_genome()
  if (is.null(current_genome) || host_bin == "") {
    # no genome selected or cleared input - allow auto-update
    poly_host_manually_set(FALSE)
  } else if (host_bin != current_genome) {
    # user typed something different - mark as manually set
    poly_host_manually_set(TRUE)
  } else {
    # matches selected genome - allow auto-update
    poly_host_manually_set(FALSE)
  }
  
  # if in host mode, filter automatically
  if (poly_auto_update() == "host") {
    raw_data <- state$raw_poly_data
    
    # if raw_data is null, try to load it
    if (is.null(raw_data) && !is.null(state$assembly)) {
      raw_data <- load_unified_data_for_assembly(state$assembly)
      state$raw_poly_data <- raw_data
    }
    
    if (!is.null(raw_data)) {
      process_filtered_poly_data(raw_data, "host", host_bin,
                                poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
    }
  }
})

# observer for last selected genome - update host bin if not manually set
observeEvent(last_selected_genome(), {
  # check if poly tab is registered
  tab <- get_tab_by_id("poly")
  if (is.null(tab)) {
    return()
  }
  
  # only update if user hasn't manually changed it
  if (!poly_host_manually_set()) {
    selected_genome <- last_selected_genome()
    if (!is.null(selected_genome) && selected_genome != "") {
      # update the reactive value and input
      poly_host_bin(selected_genome)
      updateTextInput(session, "polyHostInput", value = selected_genome)
      cache_set("poly.host_bin", selected_genome)
      
      # if in host mode, filter automatically
      if (poly_auto_update() == "host") {
        raw_data <- state$raw_poly_data
        
        # if raw_data is null, try to load it
        if (is.null(raw_data) && !is.null(state$assembly)) {
          raw_data <- load_unified_data_for_assembly(state$assembly)
          state$raw_poly_data <- raw_data
        }
        
        if (!is.null(raw_data)) {
          process_filtered_poly_data(raw_data, "host", selected_genome,
                                    poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
        }
      }
    }
  }
})

# observer for min span change
observeEvent(input$polyMinSpanInput, {
  # validate min_span - must be a valid numeric value
  min_span_numeric <- suppressWarnings(as.numeric(input$polyMinSpanInput))
  if (is.na(min_span_numeric) || !is.finite(min_span_numeric)) {
    # reset to 0 if invalid
    min_span_value <- 0
    updateNumericInput(session, "polyMinSpanInput", value = 0)
    showNotification("Invalid min span value, reset to 0", type = "warning", duration = 3)
  } else {
    min_span_value <- min_span_numeric
  }
  
  poly_min_span(min_span_value)
  cache_set("poly.min_span", min_span_value)
  
  # re-filter data
  raw_data <- state$raw_poly_data
  if (!is.null(raw_data)) {
    process_filtered_poly_data(raw_data, poly_auto_update(), poly_host_bin(),
                              min_span_value, get_state_contigs(), state$zoom, state$assembly)
  }
})

# observer for color mode change
observeEvent(input$polyColorModeInput, {
  poly_color_mode(input$polyColorModeInput)
  cache_set("poly.color_mode", input$polyColorModeInput)
  # plots will automatically re-render
})

# observer for plot type
observeEvent(input$polyPlotTypeInput, {
  poly_plot_type(input$polyPlotTypeInput)
  cache_set("poly.plot_type", input$polyPlotTypeInput)
})

# observer for x library
observeEvent(input$polyXLibInput, {
  poly_x_lib(input$polyXLibInput)
  cache_set("poly.x_lib", input$polyXLibInput)
})

# observer for y library
observeEvent(input$polyYLibInput, {
  poly_y_lib(input$polyYLibInput)
  cache_set("poly.y_lib", input$polyYLibInput)
})

# observer for jitter
observeEvent(input$polyJitterInput, {
  poly_jitter(input$polyJitterInput)
  cache_set("poly.jitter", input$polyJitterInput)
})

# observer for hover
observeEvent(input$polyHoverInput, {
  poly_hover(input$polyHoverInput)
  cache_set("poly.hover", input$polyHoverInput)
})

# observer for show locals checkbox
observeEvent(input$polyShowLocalsInput, {
  poly_show_locals(input$polyShowLocalsInput)
  cache_set("poly.show_locals", input$polyShowLocalsInput)
  # re-filter data with new type selection
  raw_data <- state$raw_poly_data
  if (!is.null(raw_data)) {
    process_filtered_poly_data(raw_data, poly_auto_update(), poly_host_bin(),
                              poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
  }
})

# observer for show rearrangements checkbox
observeEvent(input$polyShowRearrangementsInput, {
  poly_show_rearrangements(input$polyShowRearrangementsInput)
  cache_set("poly.show_rearrangements", input$polyShowRearrangementsInput)
  # re-filter data with new type selection
  raw_data <- state$raw_poly_data
  if (!is.null(raw_data)) {
    process_filtered_poly_data(raw_data, poly_auto_update(), poly_host_bin(),
                              poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
  }
})

# observer for show elements checkbox
observeEvent(input$polyShowElementsInput, {
  poly_show_elements(input$polyShowElementsInput)
  cache_set("poly.show_elements", input$polyShowElementsInput)
  # re-filter data with new type selection
  raw_data <- state$raw_poly_data
  if (!is.null(raw_data)) {
    process_filtered_poly_data(raw_data, poly_auto_update(), poly_host_bin(),
                              poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
  }
})

# observer for assembly/segments/zoom changes (only if auto-update is enabled)
observeEvent(list(state$assembly, state$segments, state$zoom), {
  if (poly_auto_update() == "host") {
    # still load data in host mode, but don't auto-filter
    if (!is.null(state$assembly)) {
      raw_data <- load_unified_data_for_assembly(state$assembly)
      state$raw_poly_data <- raw_data
      
      # if host is already set, filter by it
      if (!is.null(poly_host_bin()) && poly_host_bin() != "") {
        process_filtered_poly_data(raw_data, "host", poly_host_bin(),
                                  poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
      }
    }
    return()
  }
  
  # get or load unified data
  raw_data <- state$raw_poly_data
  if (is.null(raw_data) && !is.null(state$assembly)) {
    raw_data <- load_unified_data_for_assembly(state$assembly)
    state$raw_poly_data <- raw_data
  }
  
  # process and filter based on auto-update mode (only if we have data)
  if (!is.null(raw_data)) {
    process_filtered_poly_data(raw_data, poly_auto_update(), poly_host_bin(),
                              poly_min_span(), get_state_contigs(), state$zoom, state$assembly)
  }
}, ignoreNULL = FALSE, ignoreInit = FALSE)

# ---- Plot Renderers ----

# host abundance plot
output$polyHostAbundancePlot <- plotly::renderPlotly({
  if (poly_auto_update() != "host") {
    return(NULL)
  }
  
  host_bin <- poly_host_bin()
  if (is.null(host_bin) || host_bin == "") {
    return(NULL)
  }
  
  abundance_data <- load_host_abundance_data(state$assembly, host_bin)
  
  if (is.null(abundance_data)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No abundance data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  plot_data <- data.frame(
    library = factor(abundance_data$libraries, levels = abundance_data$libraries),
    abundance = abundance_data$values
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = abundance, group = 1)) +
    ggplot2::geom_line(color = "#2c7bb6", linewidth = 1.5) +
    ggplot2::geom_point(color = "#2c7bb6", size = 3) +
    ggplot2::labs(x = "Library", y = "Abundance (%)", title = sprintf("Host %s Abundance", host_bin)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, NA)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  plotly::ggplotly(p)
})

# helper function to create temporal or scatter plot for support/coverage/percent
create_temporal_plot <- function(data, y_label, title, color_mode, selected_uids = NULL, plot_type = "temporal", x_lib = NULL, y_lib = NULL, jitter_enabled = FALSE, hover_enabled = TRUE) {
  if (is.null(data) || is.null(data$table) || nrow(data$table) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # get value matrix (support, coverage, or percent)
  if (y_label == "Support") {
    value_matrix <- data$support
  } else if (y_label == "Coverage") {
    value_matrix <- data$coverage
  } else if (y_label == "Percent") {
    value_matrix <- calculate_percent(data)
  } else {
    return(NULL)
  }
  
  if (is.null(value_matrix)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # build plot data
  plot_data <- data.frame()
  lib_ids <- data$library_ids
  table_df <- data$table
  for (i in seq_len(nrow(table_df))) {
    uid <- table_df$uid[i]
    values <- value_matrix[i, ]
    desc <- if ("desc" %in% names(table_df)) table_df$desc[i] else ""
    if (is.na(desc) || desc == "") {
      desc <- uid
    }
    # split description on ";" and join with line breaks for hover text
    desc_formatted <- gsub(";", "<br>", desc)
    
    # ensure values are numeric
    values_numeric <- as.numeric(values)
    
    variant_df <- data.frame(
      library = factor(lib_ids, levels = lib_ids),
      value = values_numeric,
      uid = uid,
      type = table_df$type[i],
      is_sweeping = table_df$is_sweeping[i],
      desc = desc_formatted,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, variant_df)
  }
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # handle scatter mode
  if (plot_type == "scatter" && !is.null(x_lib) && !is.null(y_lib)) {
    # get library indices
    x_lib_idx <- match(x_lib, lib_ids)
    y_lib_idx <- match(y_lib, lib_ids)
    
    if (is.na(x_lib_idx) || is.na(y_lib_idx)) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Invalid library selection", size = 5) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
      return(plotly::ggplotly(p))
    }
    
    # create scatter plot data (one point per variant)
    scatter_data <- data.frame()
    for (i in seq_len(nrow(table_df))) {
      uid <- table_df$uid[i]
      x_value <- suppressWarnings(as.numeric(value_matrix[i, x_lib_idx]))
      y_value <- suppressWarnings(as.numeric(value_matrix[i, y_lib_idx]))
      desc <- if ("desc" %in% names(table_df)) table_df$desc[i] else ""
      if (is.na(desc) || desc == "") {
        desc <- uid
      }
      desc_formatted <- gsub(";", "<br>", desc)
      
      # skip invalid values
      if (!is.finite(x_value) || is.na(x_value) || !is.finite(y_value) || is.na(y_value)) {
        next
      }
      
      variant_df <- data.frame(
        uid = uid,
        x_value = x_value,
        y_value = y_value,
        type = table_df$type[i],
        is_sweeping = table_df$is_sweeping[i],
        desc = desc_formatted,
        stringsAsFactors = FALSE
      )
      scatter_data <- rbind(scatter_data, variant_df)
    }
    
    if (nrow(scatter_data) == 0) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
      return(plotly::ggplotly(p))
    }
    
    # always apply y-direction jitter (5% of range) with fixed seed for reproducibility
    y_range <- max(scatter_data$y_value) - min(scatter_data$y_value)
    set.seed(42)
    y_jitter <- y_range * 0.05
    scatter_data$y_value <- scatter_data$y_value + runif(nrow(scatter_data), -y_jitter, y_jitter)
    
    # apply x-direction jitter if enabled (with fixed seed for reproducibility)
    if (jitter_enabled) {
      x_range <- max(scatter_data$x_value) - min(scatter_data$x_value)
      set.seed(42)
      x_jitter <- x_range * 0.05
      scatter_data$x_value <- scatter_data$x_value + runif(nrow(scatter_data), -x_jitter, x_jitter)
    }
    
    # clamp percent values to 0-100 range
    if (y_label == "Percent") {
      scatter_data$x_value <- pmax(0, pmin(100, scatter_data$x_value))
      scatter_data$y_value <- pmax(0, pmin(100, scatter_data$y_value))
    }
    
    # determine colors based on color mode
    unique_uids <- unique(scatter_data$uid)
    if (color_mode == "sweep") {
      colors <- ifelse(scatter_data$is_sweeping, "red", "gray")
      names(colors) <- scatter_data$uid
    } else if (color_mode == "type") {
      type_colors <- c("local" = "grey", "rearrange" = "orange", "element" = "blue")
      uid_to_type <- unique(scatter_data[, c("uid", "type")])
      uid_to_type <- setNames(uid_to_type$type, uid_to_type$uid)
      colors <- type_colors[uid_to_type[unique_uids]]
      colors[is.na(colors)] <- "black"
      names(colors) <- unique_uids
    } else {
      colors <- rainbow(length(unique_uids))
      names(colors) <- unique_uids
    }
    
    # create hover text only if enabled
    if (hover_enabled) {
      scatter_data$hover_text <- paste0(
        scatter_data$uid, "<br>",
        x_lib, " ", y_label, ": ", round(scatter_data$x_value, 3), "<br>",
        y_lib, " ", y_label, ": ", round(scatter_data$y_value, 3), "<br>",
        scatter_data$desc
      )
    } else {
      scatter_data$hover_text <- ""
    }
    
    # check for selected variants (handle both single UID and vector)
    scatter_data$is_selected <- FALSE
    if (!is.null(selected_uids)) {
      if (is.character(selected_uids) && length(selected_uids) == 1) {
        scatter_data$is_selected <- scatter_data$uid == selected_uids
      } else if (length(selected_uids) > 0) {
        scatter_data$is_selected <- scatter_data$uid %in% selected_uids
      }
    }
    
    # create plot with proper layering: locals (bottom), rearrangements (middle), elements (top)
    non_selected_data <- scatter_data[!scatter_data$is_selected, ]
    selected_data <- scatter_data[scatter_data$is_selected, ]
    
    # split non-selected by type for proper layering
    non_selected_locals <- non_selected_data[non_selected_data$type == "local", ]
    non_selected_rearrange <- non_selected_data[non_selected_data$type == "rearrange", ]
    non_selected_elements <- non_selected_data[non_selected_data$type == "element", ]
    
    p <- ggplot2::ggplot(scatter_data, ggplot2::aes(x = x_value, y = y_value, color = uid, key = uid, text = hover_text)) +
      ggplot2::scale_color_manual(values = colors)
    
    # draw in order: locals (bottom), rearrangements (middle), elements (top)
    if (nrow(non_selected_locals) > 0) {
      p <- p + ggplot2::geom_point(data = non_selected_locals, size = 2, alpha = 0.8)
    }
    if (nrow(non_selected_rearrange) > 0) {
      p <- p + ggplot2::geom_point(data = non_selected_rearrange, size = 2, alpha = 0.8)
    }
    if (nrow(non_selected_elements) > 0) {
      p <- p + ggplot2::geom_point(data = non_selected_elements, size = 2, alpha = 0.8)
    }
    
    if (nrow(selected_data) > 0) {
      p <- p + ggplot2::geom_point(data = selected_data, size = 4, alpha = 1, 
                                  stroke = 1.5, shape = 21, fill = "white")
    }
    
    p <- p +
      ggplot2::labs(
        x = paste(x_lib, y_label),
        y = paste(y_lib, y_label)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none"
      )
    
    # format axes
    if (y_label == "Percent") {
      p <- p + ggplot2::scale_x_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
               ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100))
    }
    
    cat("poly: creating ggplotly for", y_label, "scatter plot\n")
    plotly_obj <- plotly::ggplotly(p, source = paste0("poly_", tolower(y_label), "_plot"), tooltip = "text")
    cat("poly: finished creating ggplotly for", y_label, "scatter plot\n")
    
    return(plotly_obj)
  }
  
  # temporal mode (original code)
  # determine colors based on color mode
  unique_uids <- unique(plot_data$uid)
  if (color_mode == "sweep") {
    # sweeps in red, others in gray
    colors <- ifelse(plot_data$is_sweeping, "red", "gray")
    names(colors) <- plot_data$uid
  } else if (color_mode == "type") {
    # color by type: locals=grey (bottom), rearrangements=orange (middle), elements=blue (top)
    type_colors <- c("local" = "grey", "rearrange" = "orange", "event" = "blue")
    # get type for each unique uid (all rows for same uid have same type)
    uid_to_type <- unique(plot_data[, c("uid", "type")])
    uid_to_type <- setNames(uid_to_type$type, uid_to_type$uid)
    # assign colors based on type
    colors <- type_colors[uid_to_type[unique_uids]]
    colors[is.na(colors)] <- "black"  # default to black if type not found
    names(colors) <- unique_uids
    
    # order plot_data by type for proper layering: locals first (bottom), then rearrangements, then elements (top)
    type_order <- c("local" = 1, "rearrange" = 2, "element" = 3)
    plot_data$type_order <- type_order[plot_data$type]
    plot_data$type_order[is.na(plot_data$type_order)] <- 0
    plot_data <- plot_data[order(plot_data$type_order, plot_data$uid), ]
  } else {
    # rainbow colors
    colors <- rainbow(length(unique_uids))
    names(colors) <- unique_uids
  }
  
  # highlight selected variants (handle both single UID and vector)
  has_selection <- FALSE
  if (!is.null(selected_uids)) {
    if (is.character(selected_uids) && length(selected_uids) == 1) {
      has_selection <- TRUE
      selected_data <- plot_data[plot_data$uid == selected_uids, ]
      other_data <- plot_data[plot_data$uid != selected_uids, ]
    } else if (length(selected_uids) > 0) {
      has_selection <- TRUE
      selected_data <- plot_data[plot_data$uid %in% selected_uids, ]
      other_data <- plot_data[!plot_data$uid %in% selected_uids, ]
    }
  }
  
  if (has_selection) {
    # split other_data by type for proper layering
    other_locals <- other_data[other_data$type == "local", ]
    other_rearrange <- other_data[other_data$type == "rearrange", ]
    other_elements <- other_data[other_data$type == "element", ]
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = uid, group = uid, key = uid, text = desc)) +
      ggplot2::scale_color_manual(values = colors)
    
    # draw in order: locals (bottom), rearrangements (middle), elements (top)
    if (nrow(other_locals) > 0) {
      p <- p + ggplot2::geom_line(data = other_locals, linewidth = 1, alpha = 0.6) +
        ggplot2::geom_point(data = other_locals, size = 2)
    }
    if (nrow(other_rearrange) > 0) {
      p <- p + ggplot2::geom_line(data = other_rearrange, linewidth = 1, alpha = 0.6) +
        ggplot2::geom_point(data = other_rearrange, size = 2)
    }
    if (nrow(other_elements) > 0) {
      p <- p + ggplot2::geom_line(data = other_elements, linewidth = 1, alpha = 0.6) +
        ggplot2::geom_point(data = other_elements, size = 2)
    }
    
    if (nrow(selected_data) > 0) {
      p <- p + ggplot2::geom_line(data = selected_data, linewidth = 3, alpha = 0.9) +
        ggplot2::geom_point(data = selected_data, size = 4)
    }
  } else {
    # split plot_data by type for proper layering
    plot_locals <- plot_data[plot_data$type == "local", ]
    plot_rearrange <- plot_data[plot_data$type == "rearrange", ]
    plot_elements <- plot_data[plot_data$type == "element", ]
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = uid, group = uid, key = uid, text = desc)) +
      ggplot2::scale_color_manual(values = colors)
    
    # draw in order: locals (bottom), rearrangements (middle), elements (top)
    if (nrow(plot_locals) > 0) {
      p <- p + ggplot2::geom_line(data = plot_locals, linewidth = 1, alpha = 0.7) +
        ggplot2::geom_point(data = plot_locals, size = 2)
    }
    if (nrow(plot_rearrange) > 0) {
      p <- p + ggplot2::geom_line(data = plot_rearrange, linewidth = 1, alpha = 0.7) +
        ggplot2::geom_point(data = plot_rearrange, size = 2)
    }
    if (nrow(plot_elements) > 0) {
      p <- p + ggplot2::geom_line(data = plot_elements, linewidth = 1, alpha = 0.7) +
        ggplot2::geom_point(data = plot_elements, size = 2)
    }
  }
  
  p <- p +
    ggplot2::labs(x = "Library", y = y_label, color = "Variant") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  if (y_label == "Percent") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100))
  }
  
  # convert to plotly
  cat("poly: creating ggplotly for", y_label, "temporal plot\n")
  plotly_obj <- plotly::ggplotly(p, source = paste0("poly_", tolower(y_label), "_plot"), tooltip = c("text", "x", "y"))
  cat("poly: finished creating ggplotly for", y_label, "temporal plot\n")
  
  # customize hover text to show uid, value, and description
  # create mapping from (uid, library) to hover text
  # create hover text only if enabled
  if (hover_enabled) {
    # hover text format: "uid<br>value: X<br>description"
    hover_texts <- character(nrow(plot_data))
    for (i in seq_len(nrow(plot_data))) {
      uid <- plot_data$uid[i]
      value <- plot_data$value[i]
      desc <- plot_data$desc[i]
      # ensure value is numeric before rounding
      value_num <- as.numeric(value)
      if (is.na(value_num)) {
        value_str <- "N/A"
      } else {
        value_str <- as.character(round(value_num, 3))
      }
      hover_texts[i] <- paste0(uid, "<br>", y_label, ": ", value_str, "<br>", desc)
    }
    plot_data$hover_text <- hover_texts
  } else {
    plot_data$hover_text <- ""
  }
  
  # update hover text for each trace
  for (i in seq_along(plotly_obj$x$data)) {
    if (!is.null(plotly_obj$x$data[[i]]$key)) {
      # get uid from key (which is set in the aes)
      trace_uids <- unique(plotly_obj$x$data[[i]]$key)
      if (length(trace_uids) > 0) {
        trace_uid <- trace_uids[1]
        # find matching rows in plot_data for this uid
        matching_rows <- plot_data$uid == trace_uid
        if (any(matching_rows)) {
          # set hover text - plotly will use the text for each point
          # we need to match by the order of points in the trace
          trace_data <- plot_data[matching_rows, ]
          # plotly traces are in library order, so we need to match that
          if (nrow(trace_data) > 0) {
            # sort by library to match plotly's order
            trace_data <- trace_data[order(trace_data$library), ]
            plotly_obj$x$data[[i]]$text <- trace_data$hover_text
          }
        }
      }
    }
  }
  
  return(plotly_obj)
}

# support plot
output$polySupportPlot <- plotly::renderPlotly({
  # depend on host input to trigger re-render when host changes
  poly_host_bin()
  poly_auto_update()
  filtered_data <- state$filtered_poly_data
  create_temporal_plot(filtered_data, "Support", "", poly_color_mode(), 
                      poly_selected_uids(), poly_plot_type(), poly_x_lib(), poly_y_lib(), poly_jitter(), poly_hover())
})

# coverage plot
output$polyCoveragePlot <- plotly::renderPlotly({
  # depend on host input to trigger re-render when host changes
  poly_host_bin()
  poly_auto_update()
  filtered_data <- state$filtered_poly_data
  create_temporal_plot(filtered_data, "Coverage", "", poly_color_mode(),
                      poly_selected_uids(), poly_plot_type(), poly_x_lib(), poly_y_lib(), poly_jitter(), poly_hover())
})

# percent plot
output$polyPercentPlot <- plotly::renderPlotly({
  # depend on host input to trigger re-render when host changes
  poly_host_bin()
  poly_auto_update()
  filtered_data <- state$filtered_poly_data
  create_temporal_plot(filtered_data, "Percent", "", poly_color_mode(),
                      poly_selected_uids(), poly_plot_type(), poly_x_lib(), poly_y_lib(), poly_jitter(), poly_hover())
})

# ---- Table Renderers ----

# helper function to render a variant table
render_variant_table <- function(variant_data, table_type) {
  if (is.null(variant_data) || is.null(variant_data$table) || nrow(variant_data$table) == 0) {
    return(datatable(
      data.frame(Message = sprintf("No %s variants found", table_type)),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  display_df <- variant_data$table
  
  # add is_sweep column (derived from is_sweeping)
  display_df$is_sweep <- display_df$is_sweeping
  
  # create coords column from coord1 and coord2
  if ("coord1" %in% names(display_df) && "coord2" %in% names(display_df)) {
    display_df$coords <- sapply(seq_len(nrow(display_df)), function(i) {
      coord1 <- display_df$coord1[i]
      coord2 <- display_df$coord2[i]
      if (is.na(coord1) || is.na(coord2)) {
        return("")
      }
      if (coord1 == coord2) {
        return(as.character(coord1))
      } else {
        return(sprintf("%d-%d", coord1, coord2))
      }
    })
  }
  
  # select columns for display
  col_names <- c(
    "UID" = "uid",
    "AID" = "aid",
    "Type" = "type",
    "Host Bin" = "host_bin",
    "Contig" = "contig",
    "Coords" = "coords",
    "Variant Length" = "variant_length",
    "Is Sweep" = "is_sweep",
    "Description" = "desc",
    "Data" = "data"
  )
  
  # filter to available columns
  available_cols <- names(display_df)
  col_names <- col_names[col_names %in% available_cols]
  
  # add select button column
  display_df$Actions <- sapply(seq_len(nrow(display_df)), function(i) {
    uid <- display_df$uid[i]
    as.character(actionButton(
      paste0("select_uid_btn_", table_type, "_", i),
      "Select",
      class = "btn btn-secondary btn-sm",
      onclick = sprintf("Shiny.setInputValue('polySelectUidTrigger', '%s', {priority: 'event'})", uid)
    ))
  })
  
  # create datatable
  dt <- datatable(
    display_df[, c(col_names, "Actions"), drop = FALSE],
    rownames = FALSE,
    colnames = c(names(col_names), "Actions"),
    escape = FALSE,
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip",
      ordering = TRUE,
      columnDefs = list(
        list(targets = length(col_names), orderable = FALSE, width = "80px")
      )
    ),
    selection = list(mode = "single", target = "row"),
    filter = "top"
  )
  
  # add is_sweep column styling (red for sweeps)
  if ("is_sweep" %in% names(display_df)) {
    dt <- dt %>% formatStyle("is_sweep", 
                            backgroundColor = styleEqual(TRUE, "lightcoral"))
  }
  
  dt
}

# unified variants table
output$polyVariantsTable <- renderDT({
  filtered_data <- state$filtered_poly_data
  if (is.null(filtered_data)) {
    return(datatable(
      data.frame(Message = "No data loaded"),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  render_variant_table(filtered_data, "variant")
})

# ---- Plot Click Handlers ----

# helper function to handle plot clicks and select in table
handle_plot_click <- function(uid) {
  # set selected UID (single selection)
  poly_selected_uids(uid)
  
  filtered_data <- state$filtered_poly_data
  if (is.null(filtered_data) || is.null(filtered_data$table)) {
    return()
  }
  
  # find row in unified table
  matching_rows <- which(filtered_data$table$uid == uid)
  if (length(matching_rows) > 0) {
    matching_row <- matching_rows[1]
    # select the row (single selection)
    proxy <- DT::dataTableProxy("polyVariantsTable")
    DT::selectRows(proxy, matching_row)
  }
}

# click observer for support plot
observeEvent(plotly::event_data("plotly_click", source = "poly_support_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "poly_support_plot")
  if (!is.null(event_data) && !is.null(event_data$key)) {
    handle_plot_click(event_data$key)
  }
})

# click observer for coverage plot
observeEvent(plotly::event_data("plotly_click", source = "poly_coverage_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "poly_coverage_plot")
  if (!is.null(event_data) && !is.null(event_data$key)) {
    handle_plot_click(event_data$key)
  }
})

# click observer for percent plot
observeEvent(plotly::event_data("plotly_click", source = "poly_percent_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "poly_percent_plot")
  if (!is.null(event_data) && !is.null(event_data$key)) {
    handle_plot_click(event_data$key)
  }
})

# observer for select UID button
observeEvent(input$polySelectUidTrigger, {
  uid <- as.character(input$polySelectUidTrigger)
  if (!is.null(uid) && uid != "") {
    poly_selected_uids(uid)
    handle_plot_click(uid)
  }
})

# observer for unified table row selection
observeEvent(input$polyVariantsTable_rows_selected, {
  filtered_data <- state$filtered_poly_data
  selected_rows <- input$polyVariantsTable_rows_selected
  
  if (is.null(selected_rows) || length(selected_rows) == 0) {
    poly_selected_uids(NULL)
    return()
  }
  
  if (!is.null(filtered_data) && !is.null(filtered_data$table)) {
    valid_rows <- selected_rows[selected_rows <= nrow(filtered_data$table)]
    if (length(valid_rows) > 0) {
      selected_uid <- filtered_data$table$uid[valid_rows[1]]
      poly_selected_uids(selected_uid)
    }
  }
}, ignoreNULL = FALSE)

# helper function to get selected variant from any table
get_selected_poly_variant <- function() {
  selected_uid <- poly_selected_uids()
  if (is.null(selected_uid)) {
    return(NULL)
  }
  
  filtered_data <- state$filtered_poly_data
  if (is.null(filtered_data) || is.null(filtered_data$table)) {
    return(NULL)
  }
  
  matching_rows <- which(filtered_data$table$uid == selected_uid)
  if (length(matching_rows) == 0) {
    return(NULL)
  }
  
  return(filtered_data$table[matching_rows[1], ])
}

# goto site button handler
observeEvent(input$polyGotoSiteBtn, {
  variant <- get_selected_poly_variant()
  
  if (is.null(variant)) {
    showNotification("Please select a variant to navigate to", type = "warning")
    return()
  }
  
  if (!all(c("contig", "coord1", "coord2") %in% names(variant))) {
    showNotification("Variant missing required columns (contig, coord1, coord2)", type = "error")
    return()
  }
  
  contig_name <- as.character(variant$contig)
  coord1 <- as.numeric(variant$coord1)
  coord2 <- as.numeric(variant$coord2)
  
  # push current region to undo before changing
  regions_module_output$push_undo_state()
  
  # set segments to selected contig's segments FIRST (before converting coordinates)
  segments <- get_segments(state$assembly)
  selected_segments <- segments[segments$contig == contig_name, ]
  state$segments <- selected_segments
  
  # explicitly update the context before converting coordinates
  cxt_set_view(selected_segments)
  
  # convert to global coordinates (now safe since segments are set)
  gstart <- cxt_contig2global(contig_name, min(coord1, coord2))
  gend <- cxt_contig2global(contig_name, max(coord1, coord2))
  
  # add 5kb margin on both sides
  margin <- 5000
  zoom_start <- gstart - margin
  zoom_end <- gend + margin
  
  state$zoom <- c(zoom_start, zoom_end)
  
  showNotification(sprintf("Navigated to variant %s", variant$uid), type = "message")
})

# goto element button handler
observeEvent(input$polyGotoElementBtn, {
  variant <- get_selected_poly_variant()
  
  if (is.null(variant)) {
    showNotification("Please select a variant to navigate to", type = "warning")
    return()
  }
  
  # only works for elements
  if (variant$type != "element") {
    showNotification("Goto Element only works for element variants", type = "warning")
    return()
  }
  
  # extract csegment from data field (format: "csegment ; side")
  if (!"data" %in% names(variant) || is.na(variant$data) || variant$data == "") {
    showNotification("Variant missing data field", type = "error")
    return()
  }
  
  data_str <- as.character(variant$data)
  data_parts <- strsplit(data_str, " ; ", fixed = TRUE)[[1]]
  
  if (length(data_parts) < 1) {
    showNotification("Cannot extract csegment from variant data", type = "error")
    return()
  }
  
  csegment <- data_parts[1]
  
  # load cseg map and filter by assembly
  cseg_map <- get_data("BINNING_CSEG_MAP", null.on.missing = TRUE)
  if (is.null(cseg_map)) {
    showNotification("CSEG map not available", type = "error")
    return()
  }
  
  # filter by assembly (aid)
  aid <- variant$aid
  cseg_map <- cseg_map[cseg_map$assembly == aid, ]
  
  if (nrow(cseg_map) == 0) {
    showNotification(sprintf("No cseg map data for assembly %s", aid), type = "error")
    return()
  }
  
  # find all segments for this csegment
  matching_segments <- cseg_map[cseg_map$csegment == csegment, ]
  
  if (nrow(matching_segments) == 0) {
    showNotification(sprintf("No segments found for csegment %s", csegment), type = "error")
    return()
  }
  
  segment_ids <- matching_segments$segment
  
  # get segment table to get coordinates
  segments <- get_segments(state$assembly)
  if (is.null(segments)) {
    showNotification("Segment table not available", type = "error")
    return()
  }
  
  # filter to matching segments
  matching_seg_data <- segments[segments$segment %in% segment_ids, ]
  
  if (nrow(matching_seg_data) == 0) {
    showNotification("No segment data found for matching segments", type = "error")
    return()
  }
  
  # push current region to undo before changing
  regions_module_output$push_undo_state()
  
  # set segments FIRST (before converting coordinates)
  # get all unique contigs from matching segments
  unique_contigs <- unique(matching_seg_data$contig)
  selected_segments <- segments[segments$contig %in% unique_contigs, ]
  state$segments <- selected_segments
  
  # explicitly update the context before converting coordinates
  cxt_set_view(selected_segments)
  
  # now convert all segments to global coordinates and calculate spanning range
  all_gcoords <- numeric(0)
  
  for (i in seq_len(nrow(matching_seg_data))) {
    seg_contig <- matching_seg_data$contig[i]
    seg_start <- matching_seg_data$start[i]
    seg_end <- matching_seg_data$end[i]
    
    # now that segments are set, we can convert coordinates
    gstart <- cxt_contig2global(seg_contig, seg_start)
    gend <- cxt_contig2global(seg_contig, seg_end)
    
    all_gcoords <- c(all_gcoords, gstart, gend)
  }
  
  # calculate spanning range with 10% margin
  min_coord <- min(all_gcoords)
  max_coord <- max(all_gcoords)
  span <- max_coord - min_coord
  margin <- span * 0.1
  zoom_start <- min_coord - margin
  zoom_end <- max_coord + margin
  
  # set zoom
  state$zoom <- c(zoom_start, zoom_end)
  
  showNotification(sprintf("Navigated to %d segments for csegment %s", nrow(matching_seg_data), csegment), type = "message")
})

# clear selection button handler
observeEvent(input$polyClearSelectionBtn, {
  # clear table selection
  proxy <- DT::dataTableProxy("polyVariantsTable")
  DT::selectRows(proxy, NULL)
  
  # clear reactive selection
  poly_selected_uids(NULL)
  
  showNotification("Cleared selection", type = "message")
})

