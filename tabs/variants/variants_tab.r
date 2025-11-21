# load variant utilities (only loaded when variants tab is registered)
source("tabs/variants/variants_utils.r", local = TRUE)
source("core/frequency_plots.r")

# validate and extract tab parameters
tab <- get_tab_by_id("variants")
if (is.null(tab)) {
  stop("variants tab not found during loading")
}

# check if dynamic or static mode
is_dynamic <- tab$is.dynamic %||% TRUE  # default to dynamic for backward compatibility

if (is_dynamic) {
  # dynamic mode: requires alignment functions
  required_params <- c("min_reads", "min_coverage", "min_libraries", "get_aln_f", "library_ids")
  missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
  if (length(missing_params) > 0) {
    stop(sprintf("variants tab (dynamic mode) missing required parameters: %s", paste(missing_params, collapse = ", ")))
  }
} else {
  # static mode: requires file getter functions
  required_params <- c("library_ids", "get_variants_table_f", "get_variants_support_f", "get_variants_coverage_f")
  missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
  if (length(missing_params) > 0) {
    stop(sprintf("variants tab (static mode) missing required parameters: %s", paste(missing_params, collapse = ", ")))
  }
}

# extract common parameters
library_ids <- tab$library_ids
if (!is.character(library_ids) || length(library_ids) == 0) {
  stop("library_ids must be a non-empty character vector")
}

# optional gene parameters
get_gene_table_f <- tab$get_gene_table_f
codon_table_path <- tab$codon_table_path
get_fasta_f <- tab$get_fasta_f
use_genes <- tab$use_genes

# extract mode-specific parameters
if (is_dynamic) {
  min_reads <- tab$min_reads
  min_coverage <- tab$min_coverage
  min_libraries <- tab$min_libraries
  get_aln_f <- tab$get_aln_f
  
  if (!is.function(get_aln_f)) {
    stop(sprintf("get_aln_f must be a function, got: %s", class(get_aln_f)))
  }
} else {
  # static mode parameters
  get_variants_table_f <- tab$get_variants_table_f
  get_variants_support_f <- tab$get_variants_support_f
  get_variants_coverage_f <- tab$get_variants_coverage_f
  
  if (!is.function(get_variants_table_f)) {
    stop(sprintf("get_variants_table_f must be a function, got: %s", class(get_variants_table_f)))
  }
  if (!is.function(get_variants_support_f)) {
    stop(sprintf("get_variants_support_f must be a function, got: %s", class(get_variants_support_f)))
  }
  if (!is.function(get_variants_coverage_f)) {
    stop(sprintf("get_variants_coverage_f must be a function, got: %s", class(get_variants_coverage_f)))
  }
}

# set the tab panel UI
set_tab_panel_f(function() {
  tabPanel(
    "Variants",
    # top row: update controls (25%) + plot (75%)
    fluidRow(
      column(3,
        wellPanel(
          h5("Variant Controls"),
          if (is_dynamic) {
            list(
              actionButton("updateVariantsBtn", "Update Variants", class = "btn-primary", width = "100%"),
              br(), br()
            )
          },
          verbatimTextOutput("variantCountText", placeholder = TRUE),
          br(),
          checkboxInput("autoUpdateProfilesChk", "Auto-update profiles", 
                       value = cache_get_if_exists("auto_update_profiles", FALSE), width = "100%"),
          br(), br(),
          h5("Filtering"),
          numericInput("variantSpanFilter", "Min Span:", 
                      value = 0.5, min = 0, max = 1, step = 0.1, width = "100%")
        )
      ),
      column(9,
        create_frequency_plot_ui("variant", library_ids)
      )
    ),
    # bottom row: variants table (full width)
    fluidRow(
      column(12,
        h4("Variants"),
        fluidRow(
          column(12,
            actionButton("gotoVariantsBtn", "Goto", class = "btn-secondary"),
            actionButton("clearVariantsBtn", "Clear Selection", class = "btn-secondary"),
            br(), br()
          )
        ),
        DTOutput("variantsTable")
      )
    )
  )
})

# ---- Variant Data Management ----
# (utility functions moved to variants_utils.r)

# query function that uses appropriate loading method based on mode
query_variants <- function(assembly, contigs, zoom) {
  if (is_dynamic) {
    # dynamic mode: use alignment query
    tab_config <- list(
      min_reads = min_reads,
      min_coverage = min_coverage,
      min_libraries = min_libraries,
      get_aln_f = get_aln_f,
      library_ids = library_ids,
      use_genes = use_genes,
      get_gene_table_f = get_gene_table_f,
      get_fasta_f = get_fasta_f,
      codon_table_path = codon_table_path
    )
    
    return(query_variants_for_context(assembly, contigs, zoom, tab_config))
  } else {
    # static mode: load from files
    tab_config <- list(
      library_ids = library_ids,
      get_variants_table_f = get_variants_table_f,
      get_variants_support_f = get_variants_support_f,
      get_variants_coverage_f = get_variants_coverage_f
    )
    
    return(load_variants_from_files(assembly, contigs, zoom, tab_config))
  }
}

# ---- Event Handlers ----

# reactive value to track selected variants for highlighting
selected_variants <- reactiveVal(NULL)

# observer for table row selection
observeEvent(input$variantsTable_rows_selected, {
  variant_data <- state$filtered_variant_data
  selected_rows <- input$variantsTable_rows_selected
  
  # always check if we have any selected rows
  if (is.null(selected_rows) || length(selected_rows) == 0) {
    # explicitly clear selection when no rows are selected
    selected_variants(NULL)
    cache_set("variants.selected", NULL)
    
    # refresh plots to remove highlighting only if auto-update is enabled
    if (input$autoUpdateProfilesChk %||% FALSE) {
      if (exists("refresh_trigger")) {
        current_val <- refresh_trigger()
        refresh_trigger(current_val + 1)
      }
    } else {
      # mark plots as needing refresh when highlighting changes but auto-update is disabled
      # but only if there was a previous selection (don't invalidate on startup with no selection)
      if (exists("invalidate_plot") && is.function(invalidate_plot)) {
        previous_selection <- selected_variants()
        if (!is.null(previous_selection)) {
          invalidate_plot()
        }
      }
    }
    return()
  }
  
  if (!is.null(variant_data) && !is.null(variant_data$variants)) {
    # filter valid row indices
    valid_rows <- selected_rows[selected_rows <= nrow(variant_data$variants)]
    
    if (length(valid_rows) > 0) {
      # get the selected variants info
      selected_vars <- variant_data$variants[valid_rows, ]
      selected_variants(selected_vars)
      # store in cache for profile access
      cache_set("variants.selected", selected_vars)
      
      # refresh plots to show highlighting only if auto-update is enabled
      if (input$autoUpdateProfilesChk %||% FALSE) {
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
    } else {
      # clear selection when no valid rows are selected
      selected_variants(NULL)
      cache_set("variants.selected", NULL)
      
      # refresh plots to remove highlighting
      if (exists("refresh_trigger")) {
        current_val <- refresh_trigger()
        refresh_trigger(current_val + 1)
      }
    }
  }
}, ignoreNULL = FALSE)

# additional observer to ensure highlighting is cleared when table selection changes
observeEvent(input$variantsTable_state, {
  # get current state
  state_info <- input$variantsTable_state
  if (!is.null(state_info) && !is.null(state_info$selected)) {
    selected_rows <- state_info$selected
    
    # if no rows are selected, ensure highlighting is cleared
    if (is.null(selected_rows) || length(selected_rows) == 0) {
      current_selection <- selected_variants()
      if (!is.null(current_selection)) {
        selected_variants(NULL)
        cache_set("variants.selected", NULL)
        
        # refresh plots to remove highlighting
        if (exists("refresh_trigger")) {
          current_val <- refresh_trigger()
          refresh_trigger(current_val + 1)
        }
      }
    }
  }
}, ignoreNULL = FALSE)

# observer for auto-update profiles checkbox
observeEvent(input$autoUpdateProfilesChk, {
  cache_set("auto_update_profiles", input$autoUpdateProfilesChk)
  
  # if auto-update is enabled and plots were out of sync, refresh now
  if (input$autoUpdateProfilesChk && exists("plot_updated") && is.function(plot_updated)) {
    if (!plot_updated() && exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
  }
})

# goto button handler  
observeEvent(input$gotoVariantsBtn, {
  variant_data <- state$filtered_variant_data
  selected_rows <- input$variantsTable_rows_selected
  
  if (is.null(selected_rows) || length(selected_rows) == 0) {
    showNotification("Please select variants to navigate to", type = "warning")
    return()
  }
  
  if (is.null(variant_data) || is.null(variant_data$variants)) {
    showNotification("No variant data available", type = "error")
    return()
  }
  
  # filter valid row indices
  valid_rows <- selected_rows[selected_rows <= nrow(variant_data$variants)]
  if (length(valid_rows) == 0) {
    showNotification("Invalid variant selection", type = "error")
    return()
  }
  
  # get selected variants
  selected_vars <- variant_data$variants[valid_rows, ]
  
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
  selected_vars$gcoord <- cxt$mapper$l2g(selected_vars$contig, selected_vars$coord)
  
  # calculate spanning range with appropriate margin
  min_coord <- min(selected_vars$gcoord)
  max_coord <- max(selected_vars$gcoord)
  
  if (length(valid_rows) == 1) {
    # single variant: minimum 10kb window
    window_size <- 10000  # 10kb minimum window
    center <- selected_vars$gcoord[1]
    half_window <- window_size / 2
    zoom_start <- center - half_window
    zoom_end <- center + half_window
  } else {
    # multiple variants: 10% margin on each side
    span <- max_coord - min_coord
    margin <- span * 0.1
    zoom_start <- min_coord - margin
    zoom_end <- max_coord + margin
  }
  
  # push current region to undo before changing
  regions_module_output$push_undo_state()
  
  # set zoom to calculated range
  state$zoom <- c(zoom_start, zoom_end)
  
  showNotification(sprintf("Navigated to %d selected variants", length(valid_rows)), type = "message")
})

# clear selection button handler  
observeEvent(input$clearVariantsBtn, {
  # clear table selection by using DT proxy
  proxy <- DT::dataTableProxy("variantsTable")
  DT::selectRows(proxy, NULL)
  
  # also clear the reactive selection
  selected_variants(NULL)
  cache_set("variants.selected", NULL)
  
  showNotification("Cleared variant selection", type = "message")
})

# common function to update variants data
update_variants_data <- function() {
  # query variants
  variant_data <- query_variants(state$assembly, state$contigs, state$zoom)
  
  # store raw data
  state$raw_variant_data <- variant_data
  
  # clear any previous variant selection when updating
  selected_variants(NULL)
  cache_set("variants.selected", NULL)
  
  # apply span filter and store filtered data
  if (!is.null(variant_data)) {
    filtered_data <- filter_variants_by_span(variant_data, input$variantSpanFilter %||% 0.5)
    state$filtered_variant_data <- filtered_data
    
    # add colors and store variants dataframe in global state and cache for profiles
    if (!is.null(filtered_data$variants)) {
      colored_variants <- add_variant_colors(filtered_data$variants)
      state$variants <- colored_variants
      cache_set("variants.current", colored_variants)
    } else {
      state$variants <- NULL
      cache_set("variants.current", NULL)
    }
    
    # refresh profile plots only if auto-update is enabled
    if (input$autoUpdateProfilesChk %||% FALSE) {
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
  } else {
    state$filtered_variant_data <- NULL
    state$variants <- NULL
    cache_set("variants.current", NULL)
  }
}

# update button handler (only for dynamic mode)
if (is_dynamic) {
  observeEvent(input$updateVariantsBtn, {
    update_variants_data()
  })
} else {
  # static mode: auto-load when assembly, contigs, or zoom change
  observeEvent(list(state$assembly, state$contigs, state$zoom), {
    # load if we have a valid assembly
    if (!is.null(state$assembly)) {
      update_variants_data()
    }
  }, ignoreNULL = FALSE, ignoreInit = FALSE)
}

# span filter change handler
observeEvent(input$variantSpanFilter, {
  # reapply filter when span changes
  if (!is.null(state$raw_variant_data)) {
    # clear selection when filter changes
    selected_variants(NULL)
    cache_set("variants.selected", NULL)
    filtered_data <- filter_variants_by_span(state$raw_variant_data, input$variantSpanFilter)
    state$filtered_variant_data <- filtered_data
    
    # add colors and update variants dataframe in global state and cache for profiles
    if (!is.null(filtered_data$variants)) {
      colored_variants <- add_variant_colors(filtered_data$variants)
      state$variants <- colored_variants
      cache_set("variants.current", colored_variants)
    } else {
      state$variants <- NULL
      cache_set("variants.current", NULL)
    }
    
    # refresh profile plots only if auto-update is enabled
    if (input$autoUpdateProfilesChk %||% FALSE) {
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
})

# ---- Output Renderers ----

# variant count text output
output$variantCountText <- renderText({
  variant_data <- state$filtered_variant_data
  raw_data <- state$raw_variant_data
  
  if (is.null(variant_data) || is.null(variant_data$variants)) {
    return("No variants loaded")
  }
  
  filtered_count <- nrow(variant_data$variants)
  
  # always show total and filtered when we have raw data
  if (!is.null(raw_data) && !is.null(raw_data$variants)) {
    total_count <- nrow(raw_data$variants)
    return(sprintf("Total: %d, Filtered: %d", total_count, filtered_count))
  }
  
  # fallback when no raw data
  if (filtered_count == 0) {
    return("0 variants found")
  }
  
  return(sprintf("%d variants found", filtered_count))
})

output$variantsTable <- renderDT({
  variant_data <- state$filtered_variant_data
  
  if (is.null(variant_data) || is.null(variant_data$variants) || nrow(variant_data$variants) == 0) {
    message <- if (is_dynamic) {
      "Click 'Update' to load variants for the current view"
    } else {
      "Select contigs to load variants"
    }
    return(datatable(
      data.frame(Message = message),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # format the variants table for display
  display_df <- variant_data$variants

  # truncate long sequences for display
  if ("sequence" %in% names(display_df)) {
    display_df$sequence <- sapply(display_df$sequence, function(seq) {
      if (is.na(seq) || nchar(seq) <= 10) {
        return(seq)
      } else {
        first_5 <- substr(seq, 1, 5)
        last_5 <- substr(seq, nchar(seq) - 4, nchar(seq))
        return(paste0(first_5, "...", last_5))
      }
    })
  }
  
  # store original descriptions for hover and truncate for display
  if ("desc" %in% names(display_df)) {
    # store original descriptions
    display_df$desc_full <- display_df$desc
    
    # truncate for display
    display_df$desc <- sapply(display_df$desc, function(desc) {
      if (is.na(desc) || nchar(desc) <= 15) {
        return(desc)
      } else {
        # truncate middle for long descriptions (keep first and last 5 characters)
        first_5 <- substr(desc, 1, 5)
        last_5 <- substr(desc, nchar(desc) - 4, nchar(desc))
        return(paste0(first_5, "...", last_5))
      }
    })
  }
  
  # round frequency to 3 decimal places
  if ("frequency" %in% names(display_df)) {
    display_df$frequency <- round(display_df$frequency, 3)
  }
  
  # if is_genic is missing, derive it from gene_desc when available
  if (!("is_genic" %in% names(display_df)) && ("gene_desc" %in% names(display_df))) {
    display_df$is_genic <- display_df$gene_desc != "none"
  }
  
  # format genic column for better display
  if ("is_genic" %in% names(display_df)) {
    display_df$is_genic <- ifelse(display_df$is_genic, "Genic", "Intergenic")
  }
  
  # create column name mapping with logical grouping
  col_names <- c(
    # Basic variant info
    "Variant ID" = "variant_id",
    "Contig" = "contig", 
    "Position" = "coord"
  )
  
  # Add gene info columns if available
  if ("is_genic" %in% names(display_df)) {
    col_names <- c(col_names, "Location" = "is_genic")
  }
  if ("gene_desc" %in% names(display_df)) {
    col_names <- c(col_names, "Gene" = "gene_desc")
  }
  
  # Add variant details
  col_names <- c(col_names, 
    "Type" = "type",
    "Description" = "desc"
  )
  
  # Add mutation description if available
  if ("mutation_desc" %in% names(display_df)) {
    col_names <- c(col_names, "AA Change" = "mutation_desc")
  }
  
  # Add statistics
  col_names <- c(col_names,
    "Libraries" = "library_count",
    "Support" = "total_support",
    "Coverage" = "total_coverage",
    "Frequency" = "frequency"
  )
  
  # filter to existing columns
  available_cols <- names(display_df)
  col_names <- col_names[col_names %in% available_cols]
  
  # create tooltips for columns
  tooltips <- c(
    "Variant ID" = "Unique identifier for the variant",
    "Contig" = "Reference sequence name",
    "Position" = "1-based position in the reference",
    "Location" = "Whether the variant is within a gene (Genic) or between genes (Intergenic)",
    "Gene" = "Description of the gene containing the variant",
    "Type" = "Type of variant (substitution, insertion, deletion)",
    "Description" = "Human-readable description of the change (hover for full text)",
    "AA Change" = "Amino acid change in the protein sequence",
    "Libraries" = "Number of libraries containing this variant",
    "Support" = "Total number of reads supporting this variant",
    "Coverage" = "Total read coverage at this position",
    "Frequency" = "Fraction of reads supporting this variant"
  )
  
  dt <- datatable(
    display_df[, col_names, drop = FALSE],
    rownames = FALSE,
    colnames = names(col_names),
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip"
    ),
    selection = list(mode = "multiple", target = "row"),
    filter = "top"
  )
  
  # style gene-related columns if present (use actual column names, not display names)
  gene_cols <- c()
  if ("is_genic" %in% names(display_df)) gene_cols <- c(gene_cols, "is_genic")
  if ("gene_desc" %in% names(display_df)) gene_cols <- c(gene_cols, "gene_desc")
  if ("mutation_desc" %in% names(display_df)) gene_cols <- c(gene_cols, "mutation_desc")
  
  if (length(gene_cols) > 0) {
    dt <- dt %>% formatStyle(
      gene_cols,
      backgroundColor = "rgba(240, 248, 255, 0.5)"
    )
  }
  
  if ("is_genic" %in% names(display_df)) {
    dt <- dt %>% formatStyle(
      "is_genic",
      color = styleEqual(
        c("Genic", "Intergenic"),
        c("#2c7bb6", "#888888")
      )
    )
  }
  
  dt
})

# ---- Plot Renderer ----

# helper function to get selected variants for highlighting
get_selected_variants <- function() {
  selected <- selected_variants()
  if (!is.null(selected) && nrow(selected) > 0) {
    # ensure the selected items have an 'id' field that matches items_df$id
    selected$id <- selected$variant_id
  }
  return(selected)
}

output$variantFrequencyPlot <- plotly::renderPlotly({
  
  # get current data
  raw_data <- state$filtered_variant_data
  
  # get current settings with defaults
  plot_type <- input$variantPlotType %||% "temporal"
  plot_value <- input$variantPlotValue %||% "frequency"
  x_lib <- input$variantXLib %||% library_ids[1]
  y_lib <- input$variantYLib %||% (if(length(library_ids) > 1) library_ids[2] else library_ids[1])
  jitter_enabled <- input$variantJitter %||% FALSE
  
  # get selected items
  selected_items <- get_selected_variants()
  
  # check if we have data
  has_data <- !is.null(raw_data) && !is.null(raw_data$variants)
  
  # prepare items dataframe if we have data
  items_df <- NULL
  if (has_data) {
    items_df <- add_variant_colors(raw_data$variants)
    items_df$id <- items_df$variant_id
    items_df$label <- paste(items_df$type, items_df$contig, items_df$coord, sep = " ")
  }
  
  # render the plot using the cleaned internal function
  no_data_message <- if (is_dynamic) {
    "Click 'Update' to load variants"
  } else {
    "Select contigs to load variants"
  }
  
  render_frequency_plot_internal(has_data, items_df, raw_data$support, raw_data$coverage,
                                plot_type, plot_value, x_lib, y_lib, 
                                jitter_enabled, selected_items, library_ids, 
                                no_data_message)
})

# frequency plot observers for caching
observeEvent(input$variantPlotType, {
  cache_set("variant_frequency_plot_type", input$variantPlotType)
})

observeEvent(input$variantPlotValue, {
  cache_set("variant_frequency_plot_value", input$variantPlotValue)
})

observeEvent(input$variantXLib, {
  cache_set("variant_frequency_plot_x_lib", input$variantXLib)
})

observeEvent(input$variantYLib, {
  cache_set("variant_frequency_plot_y_lib", input$variantYLib)
})

observeEvent(input$variantJitter, {
  cache_set("variant_frequency_plot_jitter", input$variantJitter)
})

# click observers for frequency plot interaction
observeEvent(plotly::event_data("plotly_click", source = "scatter_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "scatter_plot")
  if (!is.null(event_data) && !is.null(event_data$key)) {
    variant_data <- state$filtered_variant_data
    if (!is.null(variant_data) && !is.null(variant_data$variants)) {
      matching_rows <- which(variant_data$variants$variant_id == event_data$key)
      if (length(matching_rows) > 0) {
        proxy <- DT::dataTableProxy("variantsTable")
        DT::selectRows(proxy, matching_rows[1])
      }
    }
  }
})

observeEvent(plotly::event_data("plotly_click", source = "temporal_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "temporal_plot")
  if (!is.null(event_data) && !is.null(event_data$key)) {
    variant_data <- state$filtered_variant_data
    if (!is.null(variant_data) && !is.null(variant_data$variants)) {
      matching_rows <- which(variant_data$variants$variant_id == event_data$key)
      if (length(matching_rows) > 0) {
        proxy <- DT::dataTableProxy("variantsTable")
        DT::selectRows(proxy, matching_rows[1])
      }
    }
  }
})

# export function for PDF generation
variants_export_pdf <- function(region_info) {
  # load data fresh for export region (bypass state to avoid stale data)
  if (is_dynamic) {
    # dynamic mode: query fresh data from alntools
    tab_config <- list(
      min_reads = min_reads,
      min_coverage = min_coverage,
      min_libraries = min_libraries,
      get_aln_f = get_aln_f,
      library_ids = library_ids,
      use_genes = use_genes,
      get_gene_table_f = get_gene_table_f,
      get_fasta_f = get_fasta_f,
      codon_table_path = codon_table_path
    )
    raw_data <- query_variants_for_context(region_info$assembly, region_info$contigs, region_info$context_zoom, tab_config)
  } else {
    # static mode: load fresh data from files
    tab_config <- list(
      library_ids = library_ids,
      get_variants_table_f = get_variants_table_f,
      get_variants_support_f = get_variants_support_f,
      get_variants_coverage_f = get_variants_coverage_f
    )
    raw_data <- load_variants_from_files(region_info$assembly, region_info$contigs, NULL, tab_config)
    # filter to zoom coordinates
    raw_data <- filter_variants_by_region(raw_data, region_info$contigs, region_info$context_zoom, region_info$assembly)
  }
  
  # apply current filters
  span_filter <- input$variantSpanFilter %||% cache_get_if_exists("variant.span_filter", 0.5)
  
  # filter the data
  filtered_data <- filter_variants_by_span(raw_data, span_filter)
  
  # cache filtered variants for profile access during export
  if (!is.null(filtered_data) && !is.null(filtered_data$variants)) {
    colored_variants <- add_variant_colors(filtered_data$variants)
    cache_set("variants.current", colored_variants)
  } else {
    cache_set("variants.current", NULL)
  }
  
  # get current plot settings
  plot_type <- input$variantPlotType %||% "temporal"
  plot_value <- input$variantPlotValue %||% "frequency"
  x_lib <- input$variantXLib %||% library_ids[1]
  y_lib <- input$variantYLib %||% (if(length(library_ids) > 1) library_ids[2] else library_ids[1])
  jitter_enabled <- input$variantJitter %||% FALSE
  
  # check if we have data
  has_data <- !is.null(filtered_data) && !is.null(filtered_data$variants)
  
  # prepare items dataframe if we have data
  items_df <- NULL
  if (has_data) {
    items_df <- add_variant_colors(filtered_data$variants)
    items_df$id <- items_df$variant_id
    items_df$label <- paste(items_df$type, items_df$contig, items_df$coord, sep = " ")
  }
  
  # create the plot using the unified export function
  return(create_frequency_plot_for_export(has_data, items_df, filtered_data$support, filtered_data$coverage,
                                         plot_type, plot_value, x_lib, y_lib, 
                                         jitter_enabled, library_ids, 
                                         title = "Variants"))
}

# export function for table generation
variants_export_table <- function(region_info) {
  # load data fresh for export region (bypass state to avoid stale data)
  if (is_dynamic) {
    # dynamic mode: query fresh data from alntools
    tab_config <- list(
      min_reads = min_reads,
      min_coverage = min_coverage,
      min_libraries = min_libraries,
      get_aln_f = get_aln_f,
      library_ids = library_ids,
      use_genes = use_genes,
      get_gene_table_f = get_gene_table_f,
      get_fasta_f = get_fasta_f,
      codon_table_path = codon_table_path
    )
    raw_data <- query_variants_for_context(region_info$assembly, region_info$contigs, region_info$context_zoom, tab_config)
  } else {
    # static mode: load fresh data from files
    tab_config <- list(
      library_ids = library_ids,
      get_variants_table_f = get_variants_table_f,
      get_variants_support_f = get_variants_support_f,
      get_variants_coverage_f = get_variants_coverage_f
    )
    raw_data <- load_variants_from_files(region_info$assembly, region_info$contigs, NULL, tab_config)
    # filter to zoom coordinates
    raw_data <- filter_variants_by_region(raw_data, region_info$contigs, region_info$context_zoom, region_info$assembly)
  }
  
  # apply current filters
  span_filter <- input$variantSpanFilter %||% cache_get_if_exists("variant.span_filter", 0.5)
  
  # filter the data
  filtered_data <- filter_variants_by_span(raw_data, span_filter)
  
  # return the filtered variants dataframe (or NULL if no data)
  if (!is.null(filtered_data) && !is.null(filtered_data$variants)) {
    return(filtered_data$variants)
  }
  
  return(NULL)
}

# register the export functions for this tab
register_tab_export_function("variants", "pdf", variants_export_pdf)
register_tab_export_function("variants", "table", variants_export_table)
