# load variant utilities (only loaded when variants tab is registered)
source("tabs/variants/variants_utils.r", local = TRUE)
source("profiles/align/align_utils.r")
source("core/frequency_plots.r")

# validate and extract tab parameters
tab <- get_tab_by_id("variants")
if (is.null(tab)) {
  stop("variants tab not found during loading")
}

# check if dynamic or static mode
is_dynamic <- tab$is.dynamic %||% TRUE  # default to dynamic for backward compatibility

if (is_dynamic) {
  required_params <- c("min_reads", "min_coverage", "min_libraries", "get_aln_f", "library_ids")
  missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
  if (length(missing_params) > 0) {
    stop(sprintf("variants tab (dynamic mode) missing required parameters: %s", paste(missing_params, collapse = ", ")))
  }
} else {
  required_params <- c("get_variants_table_f", "get_variants_support_f", "get_variants_coverage_f")
  missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
  if (length(missing_params) > 0) {
    stop(sprintf("variants tab (static mode) missing required parameters: %s", paste(missing_params, collapse = ", ")))
  }
  if (is.null(tab$library_ids) && is.null(tab$library_id_map)) {
    stop("variants tab (static mode) requires either library_ids or library_id_map")
  }
}

# extract library configuration
use_library_id_map <- !is.null(tab$library_id_map)
if (use_library_id_map) {
  library_id_map <- tab$library_id_map
} else {
  library_ids_param <- tab$library_ids
  if (!is.character(library_ids_param) || length(library_ids_param) == 0) {
    stop("library_ids must be a non-empty character vector")
  }
  library_id_map <- data.frame(aid = "*", set_name = "default", stringsAsFactors = FALSE)
  library_id_map$set_ids <- list(library_ids_param)
}

has_multiple_sets <- nrow(library_id_map) > 1
all_library_ids <- unique(unlist(library_id_map$set_ids))
library_ids <- library_id_map$set_ids[[1]]

# optional gene parameters
get_gene_table_f <- tab$get_gene_table_f
codon_table_path <- tab$codon_table_path
get_fasta_f <- tab$get_fasta_f
use_genes <- tab$use_genes

# optional sample map for matrix individual indicators and dividers
sample_map <- tab$sample_map

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
          if (has_multiple_sets) {
            selectInput("variantSampleSet", "Sample Set:",
                       choices = setNames(library_id_map$set_name, library_id_map$set_name),
                       selected = cache_get_if_exists("variant_sample_set", library_id_map$set_name[1]),
                       width = "100%")
          },
          if (is_dynamic) {
            list(
              actionButton("updateVariantsBtn", "Update Variants", class = "btn-primary", width = "100%"),
              br(), br()
            )
          },
          verbatimTextOutput("variantCountText", placeholder = TRUE),
          br(),
          h5("Selected Variant"),
          textInput("selectedVariantId", label = NULL,
                    value = "", placeholder = "variant ID", width = "100%"),
          br(),
          checkboxInput("autoUpdateProfilesChk", "Auto-update profiles", 
                       value = cache_get_if_exists("auto_update_profiles", FALSE), width = "100%"),
          br(), br(),
          h5("Filtering"),
          numericInput("variantSpanFilter", "Min Span:", 
                      value = cache_get_if_exists("variant.span_filter", 0), min = 0, max = 1, step = 0.1, width = "100%"),
          numericInput("variantMinSupportFilter", "Min Support:",
                      value = cache_get_if_exists("variant.min_support_filter", 2), min = 0, step = 1, width = "100%"),
          br(),
          h6("Types"),
          checkboxGroupInput("variantTypeFilter", label = NULL,
            choices  = c("Substitution" = "sub", "Insertion" = "ins", "Deletion" = "del", "Clip" = "clip"),
            selected = cache_get_if_exists("variant.type_filter", c("sub", "ins", "del", "clip")),
            inline   = FALSE),
          br(),
          h5("Sort By"),
          selectInput("variantSortBy", label = NULL,
            choices  = c("ID" = "id", "Frequency" = "frequency"),
            selected = cache_get_if_exists("variant.sort_by", "id"),
            width    = "100%"),
          conditionalPanel(
            condition = "input.variantPlotType == 'matrix'",
            h5("Sample Order"),
            selectInput("variantColSortBy", label = NULL,
              choices  = c("ID" = "id", "Frequency" = "frequency"),
              selected = cache_get_if_exists("variant.col_sort_by", "id"),
              width    = "100%"),
            h5("Max Variants"),
            numericInput("variantMatrixMaxItems", label = NULL,
              value = cache_get_if_exists("variant.matrix_max_items", 100),
              min = 1, step = 10, width = "100%")
          )
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
            actionButton("selectVariantFromTableBtn", "Select from Table", class = "btn-secondary"),
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
    tab_config <- list(
      get_variants_table_f = get_variants_table_f,
      get_variants_support_f = get_variants_support_f,
      get_variants_coverage_f = get_variants_coverage_f
    )
    if (use_library_id_map) {
      tab_config$all_library_ids <- all_library_ids
    } else {
      tab_config$library_ids <- library_ids
    }
    return(load_variants_from_files(assembly, contigs, zoom, tab_config))
  }
}

# ---- Active Library Set ----

get_active_library_ids <- reactive({
  if (!has_multiple_sets) return(library_ids)
  selected_set <- input$variantSampleSet
  if (is.null(selected_set)) return(library_ids)
  assembly <- state$assembly %||% "*"
  valid_sets <- library_id_map[library_id_map$aid == "*" | library_id_map$aid == assembly, ]
  idx <- which(valid_sets$set_name == selected_set)
  if (length(idx) == 0) return(valid_sets$set_ids[[1]])
  valid_sets$set_ids[[idx[1]]]
})

# subset columns to active set and apply span filter
apply_and_store_filters <- function() {
  raw_data <- state$raw_variant_data
  if (is.null(raw_data) || is.null(raw_data$support)) {
    return()
  }
  
  active_ids <- get_active_library_ids()
  span_filter <- input$variantSpanFilter %||% cache_get_if_exists("variant.span_filter", 0)
  min_support_filter <- input$variantMinSupportFilter %||% cache_get_if_exists("variant.min_support_filter", 2)
  type_filter <- input$variantTypeFilter %||% cache_get_if_exists("variant.type_filter", c("sub", "ins", "del", "clip"))

  available_cols <- intersect(active_ids, colnames(raw_data$support))
  if (length(available_cols) == 0) {
    state$filtered_variant_data <- NULL
    state$variants <- NULL
    cache_set("variants.current", NULL)
    return()
  }
  
  subset_data <- list(
    variants = raw_data$variants,
    support = raw_data$support[, available_cols, drop = FALSE],
    coverage = raw_data$coverage[, available_cols, drop = FALSE],
    library_ids = available_cols
  )
  
  filtered_data <- filter_variants_by_span(subset_data, span_filter)
  filtered_data <- filter_variants_by_min_support(filtered_data, min_support_filter)
  filtered_data <- filter_variants_by_types(filtered_data, type_filter)
  state$filtered_variant_data <- filtered_data
  
  if (!is.null(filtered_data$variants) && nrow(filtered_data$variants) > 0) {
    colored_variants <- add_variant_colors(filtered_data$variants)
    state$variants <- colored_variants
    cache_set("variants.current", colored_variants)
  } else {
    state$variants <- NULL
    cache_set("variants.current", NULL)
  }

  # restore selection from the text input after filtering
  current_id <- trimws(input$selectedVariantId %||% "")
  if (nchar(current_id) > 0 && !is.null(filtered_data$variants)) {
    matching_rows <- which(filtered_data$variants$variant_id == current_id)
    if (length(matching_rows) > 0) {
      selected_vars <- filtered_data$variants[matching_rows[1], ]
      selected_vars$id <- selected_vars$variant_id
      selected_variants(selected_vars)
      cache_set("variants.selected", selected_vars)
    } else {
      selected_variants(NULL)
      cache_set("variants.selected", NULL)
    }
  } else {
    selected_variants(NULL)
    cache_set("variants.selected", NULL)
  }
  
  if (input$autoUpdateProfilesChk %||% FALSE) {
    if (exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
  } else {
    if (exists("invalidate_plot") && is.function(invalidate_plot)) {
      invalidate_plot()
    }
  }
}

# observer for sample set changes
if (has_multiple_sets) {
  observeEvent(input$variantSampleSet, {
    cache_set("variant_sample_set", input$variantSampleSet)
    active_ids <- get_active_library_ids()
    
    updateSelectInput(session, "variantXLib",
                      choices = setNames(active_ids, active_ids),
                      selected = active_ids[1])
    updateSelectInput(session, "variantYLib",
                      choices = setNames(active_ids, active_ids),
                      selected = if (length(active_ids) > 1) active_ids[2] else active_ids[1])
    
    if (!is.null(state$raw_variant_data)) {
      apply_and_store_filters()
    }
  })
}

# ---- Event Handlers ----

# reactive value to track selected variants for highlighting
selected_variants <- reactiveVal(NULL)

# apply selection by variant ID; update_text=TRUE when called from button/plot (not text input)
apply_variant_selection <- function(variant_id, update_text = TRUE) {
  variant_data <- state$filtered_variant_data
  if (is.null(variant_id) || nchar(trimws(variant_id)) == 0) {
    selected_variants(NULL)
    cache_set("variants.selected", NULL)
    if (update_text) updateTextInput(session, "selectedVariantId", value = "")
  } else {
    if (!is.null(variant_data) && !is.null(variant_data$variants)) {
      matching_rows <- which(variant_data$variants$variant_id == variant_id)
      if (length(matching_rows) > 0) {
        selected_vars <- variant_data$variants[matching_rows[1], ]
        selected_vars$id <- selected_vars$variant_id
        selected_variants(selected_vars)
        cache_set("variants.selected", selected_vars)
        if (update_text) updateTextInput(session, "selectedVariantId", value = variant_id)
        # sync table highlight
        proxy <- DT::dataTableProxy("variantsTable")
        DT::selectRows(proxy, matching_rows[1])
      } else {
        selected_variants(NULL)
        cache_set("variants.selected", NULL)
      }
    }
  }
  if (input$autoUpdateProfilesChk %||% FALSE) {
    if (exists("refresh_trigger")) refresh_trigger(refresh_trigger() + 1)
  } else {
    if (exists("invalidate_plot") && is.function(invalidate_plot)) invalidate_plot()
  }
}

# select from table button
observeEvent(input$selectVariantFromTableBtn, {
  selected_rows <- input$variantsTable_rows_selected
  variant_data  <- state$filtered_variant_data
  if (!is.null(selected_rows) && length(selected_rows) > 0 &&
      !is.null(variant_data) && !is.null(variant_data$variants)) {
    valid_row <- selected_rows[selected_rows <= nrow(variant_data$variants)][1]
    if (!is.na(valid_row)) {
      apply_variant_selection(variant_data$variants$variant_id[valid_row], update_text = TRUE)
    }
  }
})

# text input: user typed a variant ID
observeEvent(input$selectedVariantId, {
  apply_variant_selection(input$selectedVariantId, update_text = FALSE)
}, ignoreInit = TRUE)

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
  current_id   <- trimws(input$selectedVariantId %||% "")

  if (nchar(current_id) == 0) {
    showNotification("Please select a variant to navigate to", type = "warning")
    return()
  }

  if (is.null(variant_data) || is.null(variant_data$variants)) {
    showNotification("No variant data available", type = "error")
    return()
  }

  matching_rows <- which(variant_data$variants$variant_id == current_id)
  if (length(matching_rows) == 0) {
    showNotification("Selected variant not found in current data", type = "warning")
    return()
  }

  # get selected variants
  selected_vars <- variant_data$variants[matching_rows, ]
  
  # convert to global coordinates using context services
  selected_vars$gcoord <- cxt_contig2global(selected_vars$contig, selected_vars$coord)
  if (any(is.na(selected_vars$gcoord))) {
    showNotification("Selected variant(s) are not in the current view", type = "warning")
    return()
  }
  
  # calculate spanning range with appropriate margin
  min_coord <- min(selected_vars$gcoord)
  max_coord <- max(selected_vars$gcoord)
  
  if (length(matching_rows) == 1) {
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
  proxy <- DT::dataTableProxy("variantsTable")
  DT::selectRows(proxy, NULL)
  selected_variants(NULL)
  cache_set("variants.selected", NULL)
  updateTextInput(session, "selectedVariantId", value = "")
  if (exists("invalidate_plot") && is.function(invalidate_plot)) invalidate_plot()
  showNotification("Cleared variant selection", type = "message")
})

# common function to update variants data
update_variants_data <- function() {
  variant_data <- query_variants(state$assembly, get_state_contigs(), state$zoom)
  state$raw_variant_data <- variant_data
  
  selected_variants(NULL)
  cache_set("variants.selected", NULL)
  
  if (!is.null(variant_data)) {
    apply_and_store_filters()
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
  # static mode: auto-load when assembly, segments, or zoom change
  observeEvent(list(state$assembly, state$segments, state$zoom), {
    # load if we have a valid assembly
    if (!is.null(state$assembly)) {
      update_variants_data()
    }
  }, ignoreNULL = FALSE, ignoreInit = FALSE, priority = -1)
}

# span filter change handler
observeEvent(input$variantSpanFilter, {
  cache_set("variant.span_filter", input$variantSpanFilter)
  if (!is.null(state$raw_variant_data)) {
    apply_and_store_filters()
  }
})

# min support filter change handler
observeEvent(input$variantMinSupportFilter, {
  cache_set("variant.min_support_filter", input$variantMinSupportFilter)
  if (!is.null(state$raw_variant_data)) {
    apply_and_store_filters()
  }
})

# variant type filter change handler
observeEvent(input$variantTypeFilter, {
  cache_set("variant.type_filter", input$variantTypeFilter)
  if (!is.null(state$raw_variant_data)) {
    apply_and_store_filters()
  }
}, ignoreNULL = FALSE)

# sort-by change handlers
observeEvent(input$variantSortBy, {
  cache_set("variant.sort_by", input$variantSortBy)
})

observeEvent(input$variantColSortBy, {
  cache_set("variant.col_sort_by", input$variantColSortBy)
})

observeEvent(input$variantMatrixColorBy, {
  cache_set("variant_frequency_plot_matrix_color_by", input$variantMatrixColorBy)
})

observeEvent(input$variantMatrixMaxItems, {
  cache_set("variant.matrix_max_items", input$variantMatrixMaxItems)
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
    display_df$sequence <- format_long_sequence_for_display(display_df$sequence)
  }
  
  # store original descriptions for hover and truncate for display
  if ("desc" %in% names(display_df)) {
    # store original descriptions
    display_df$desc_full <- display_df$desc
    
    # truncate for display (indels: prefix + count; other long: middle ellipsis)
    display_df$desc <- sapply(display_df$desc, function(desc) {
      if (is.na(desc)) {
        return(desc)
      }
      if (grepl("^\\+|^-", desc)) {
        return(format_indel_desc_for_hover(desc))
      }
      if (nchar(desc) <= 15) {
        return(desc)
      }
      first_5 <- substr(desc, 1, 5)
      last_5 <- substr(desc, nchar(desc) - 4, nchar(desc))
      paste0(first_5, "...", last_5)
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
  raw_data <- state$filtered_variant_data
  active_ids <- get_active_library_ids()
  
  plot_type <- input$variantPlotType %||% "temporal"
  plot_value <- input$variantPlotValue %||% "frequency"
  x_lib <- input$variantXLib %||% active_ids[1]
  y_lib <- input$variantYLib %||% (if(length(active_ids) > 1) active_ids[2] else active_ids[1])
  jitter_enabled <- input$variantJitter %||% FALSE
  
  selected_items <- get_selected_variants()
  
  has_data <- !is.null(raw_data) && !is.null(raw_data$variants)
  
  items_df <- NULL
  if (has_data) {
    items_df <- add_variant_colors(raw_data$variants)
    items_df$id <- items_df$variant_id
    items_df$label <- paste(items_df$type, items_df$contig, items_df$coord, sep = " ")
  }
  
  no_data_message <- if (is_dynamic) {
    "Click 'Update' to load variants"
  } else {
    "Select contigs to load variants"
  }
  
  sort_by          <- input$variantSortBy %||% cache_get_if_exists("variant.sort_by", "id")
  col_sort_by      <- input$variantColSortBy %||% cache_get_if_exists("variant.col_sort_by", "id")
  matrix_color_by  <- input$variantMatrixColorBy %||%
                        cache_get_if_exists("variant_frequency_plot_matrix_color_by", "value")
  matrix_max_items <- input$variantMatrixMaxItems %||%
                        cache_get_if_exists("variant.matrix_max_items", 100)

  render_frequency_plot_internal(has_data, items_df, raw_data$support, raw_data$coverage,
                                plot_type, plot_value, x_lib, y_lib,
                                jitter_enabled, selected_items, active_ids,
                                no_data_message, max_items = matrix_max_items, sort_by = sort_by,
                                col_sort_by = col_sort_by, sample_map = sample_map,
                                matrix_color_by = matrix_color_by)
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
  if (!is.null(event_data) && !is.null(event_data$key))
    apply_variant_selection(event_data$key, update_text = TRUE)
})

observeEvent(plotly::event_data("plotly_click", source = "temporal_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "temporal_plot")
  if (!is.null(event_data) && !is.null(event_data$key))
    apply_variant_selection(event_data$key, update_text = TRUE)
})

observeEvent(plotly::event_data("plotly_click", source = "matrix_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "matrix_plot")
  if (!is.null(event_data) && !is.null(event_data$y))
    apply_variant_selection(event_data$y, update_text = TRUE)
})

# export function for PDF generation
variants_export_pdf <- function(region_info) {
  active_ids <- get_active_library_ids()
  
  if (is_dynamic) {
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
    tab_config <- list(
      get_variants_table_f = get_variants_table_f,
      get_variants_support_f = get_variants_support_f,
      get_variants_coverage_f = get_variants_coverage_f
    )
    if (use_library_id_map) {
      tab_config$all_library_ids <- all_library_ids
    } else {
      tab_config$library_ids <- library_ids
    }
    raw_data <- load_variants_from_files(region_info$assembly, region_info$contigs, NULL, tab_config)
    raw_data <- filter_variants_by_region(raw_data, region_info$contigs, region_info$context_zoom, region_info$assembly)
  }
  
  # subset to active set and apply span, support, and type filters
  span_filter <- input$variantSpanFilter %||% cache_get_if_exists("variant.span_filter", 0)
  min_support_filter <- input$variantMinSupportFilter %||% cache_get_if_exists("variant.min_support_filter", 2)
  type_filter <- input$variantTypeFilter %||% cache_get_if_exists("variant.type_filter", c("sub", "ins", "del", "clip"))

  if (!is.null(raw_data) && !is.null(raw_data$support)) {
    available_cols <- intersect(active_ids, colnames(raw_data$support))
    if (length(available_cols) > 0) {
      raw_data$support <- raw_data$support[, available_cols, drop = FALSE]
      raw_data$coverage <- raw_data$coverage[, available_cols, drop = FALSE]
    }
  }
  
  filtered_data <- filter_variants_by_span(raw_data, span_filter)
  filtered_data <- filter_variants_by_min_support(filtered_data, min_support_filter)
  filtered_data <- filter_variants_by_types(filtered_data, type_filter)
  
  if (!is.null(filtered_data) && !is.null(filtered_data$variants)) {
    colored_variants <- add_variant_colors(filtered_data$variants)
    cache_set("variants.current", colored_variants)
  } else {
    cache_set("variants.current", NULL)
  }
  
  plot_type <- input$variantPlotType %||% "temporal"
  plot_value <- input$variantPlotValue %||% "frequency"
  x_lib <- input$variantXLib %||% active_ids[1]
  y_lib <- input$variantYLib %||% (if(length(active_ids) > 1) active_ids[2] else active_ids[1])
  jitter_enabled <- input$variantJitter %||% FALSE
  
  has_data <- !is.null(filtered_data) && !is.null(filtered_data$variants)
  
  items_df <- NULL
  if (has_data) {
    items_df <- add_variant_colors(filtered_data$variants)
    items_df$id <- items_df$variant_id
    items_df$label <- paste(items_df$type, items_df$contig, items_df$coord, sep = " ")
  }
  
  col_sort_by     <- cache_get_if_exists("variant.col_sort_by", "id")
  sort_by         <- cache_get_if_exists("variant.sort_by", "id")
  matrix_color_by <- cache_get_if_exists("variant_frequency_plot_matrix_color_by", "value")

  matrix_max_items <- cache_get_if_exists("variant.matrix_max_items", 100)

  selected_vid <- cache_get_if_exists("variants.selected", NULL)
  selected_vid <- if (!is.null(selected_vid)) selected_vid$variant_id[1] else NULL

  return(create_frequency_plot_for_export(has_data, items_df, filtered_data$support, filtered_data$coverage,
                                         plot_type, plot_value, x_lib, y_lib,
                                         jitter_enabled, active_ids,
                                         title = "Variants", max_items = matrix_max_items,
                                         col_sort_by = col_sort_by, sort_by = sort_by,
                                         matrix_color_by = matrix_color_by,
                                         sample_map = sample_map,
                                         selected_variant_id = selected_vid))
}

# export function for table generation
variants_export_table <- function(region_info) {
  active_ids <- get_active_library_ids()
  
  if (is_dynamic) {
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
    tab_config <- list(
      get_variants_table_f = get_variants_table_f,
      get_variants_support_f = get_variants_support_f,
      get_variants_coverage_f = get_variants_coverage_f
    )
    if (use_library_id_map) {
      tab_config$all_library_ids <- all_library_ids
    } else {
      tab_config$library_ids <- library_ids
    }
    raw_data <- load_variants_from_files(region_info$assembly, region_info$contigs, NULL, tab_config)
    raw_data <- filter_variants_by_region(raw_data, region_info$contigs, region_info$context_zoom, region_info$assembly)
  }
  
  # subset to active set and apply span, support, and type filters
  span_filter <- input$variantSpanFilter %||% cache_get_if_exists("variant.span_filter", 0)
  min_support_filter <- input$variantMinSupportFilter %||% cache_get_if_exists("variant.min_support_filter", 2)
  type_filter <- input$variantTypeFilter %||% cache_get_if_exists("variant.type_filter", c("sub", "ins", "del", "clip"))

  if (!is.null(raw_data) && !is.null(raw_data$support)) {
    available_cols <- intersect(active_ids, colnames(raw_data$support))
    if (length(available_cols) > 0) {
      raw_data$support <- raw_data$support[, available_cols, drop = FALSE]
      raw_data$coverage <- raw_data$coverage[, available_cols, drop = FALSE]
    }
  }
  
  filtered_data <- filter_variants_by_span(raw_data, span_filter)
  filtered_data <- filter_variants_by_min_support(filtered_data, min_support_filter)
  filtered_data <- filter_variants_by_types(filtered_data, type_filter)
  
  if (!is.null(filtered_data) && !is.null(filtered_data$variants)) {
    return(filtered_data$variants)
  }
  
  return(NULL)
}

# register the export functions for this tab
register_tab_export_function("variants", "pdf", variants_export_pdf)
register_tab_export_function("variants", "table", variants_export_table)
