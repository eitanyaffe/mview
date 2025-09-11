# load variant utilities (only loaded when variants tab is registered)
source("tabs/variants/variants_utils.r", local = TRUE)
source("tabs/variants/variants_plots.r", local = TRUE)

# validate and extract tab parameters
tab <- get_tab_by_id("variants")
if (is.null(tab)) {
  stop("variants tab not found during loading")
}

required_params <- c("min_reads", "min_coverage", "min_libraries", "get_aln_f", "library_ids")
missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
if (length(missing_params) > 0) {
  stop(sprintf("variants tab missing required parameters: %s", paste(missing_params, collapse = ", ")))
}

# optional gene parameters
get_gene_table_f <- tab$get_gene_table_f
codon_table_path <- tab$codon_table_path
get_fasta_f <- tab$get_fasta_f
use_genes <- tab$use_genes

min_reads <- tab$min_reads
min_coverage <- tab$min_coverage
min_libraries <- tab$min_libraries
get_aln_f <- tab$get_aln_f
library_ids <- tab$library_ids

if (!is.function(get_aln_f)) {
  stop(sprintf("get_aln_f must be a function, got: %s", class(get_aln_f)))
}
if (!is.character(library_ids) || length(library_ids) == 0) {
  stop("library_ids must be a non-empty character vector")
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
          actionButton("updateVariantsBtn", "Update Variants", class = "btn-primary", width = "100%"),
          br(), br(),
          verbatimTextOutput("variantCountText", placeholder = TRUE),
          br(),
          checkboxInput("autoUpdateProfilesChk", "Auto-update profiles", 
                       value = cache_get_if_exists("auto_update_profiles", FALSE), width = "100%"),
          br(),
          actionButton("gotoVariantsBtn", "Goto", class = "btn-secondary", width = "100%"),
          br(), br(),
          h5("Filtering"),
          numericInput("variantSpanFilter", "Min Span:", 
                      value = 0.5, min = 0, max = 1, step = 0.1, width = "100%"),
          br(),
          selectInput("variantPlotMode", "Plot Mode:", 
                     choices = c("Frequency" = "frequency", 
                                "Variant Support" = "variant_support", 
                                "Variant Coverage" = "variant_coverage"),
                     selected = "frequency", width = "100%")
        )
      ),
      column(9,
        h4("Variant Plot"),
        plotly::plotlyOutput("variantFrequencyPlot", height = "400px")
      )
    ),
    # bottom row: variants table (full width)
    fluidRow(
      column(12,
        h4("Variants"),
        DTOutput("variantsTable")
      )
    )
  )
})

# ---- Variant Data Management ----
# (utility functions moved to variants_utils.r)

# simplified query function that uses the utility function
query_variants <- function(assembly, contigs, zoom) {
  # create tab config object for the utility function
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


# update button handler
observeEvent(input$updateVariantsBtn, {
  
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
    # still mark as needing refresh even if no data
    if (!(input$autoUpdateProfilesChk %||% FALSE)) {
      if (exists("invalidate_plot") && is.function(invalidate_plot)) {
        invalidate_plot()
      }
    }
  }
})

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
    return(datatable(
      data.frame(Message = "Click 'Update' to load variants for the current view"),
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

output$variantFrequencyPlot <- plotly::renderPlotly({
  variant_data <- state$filtered_variant_data
  
  if (is.null(variant_data) || is.null(variant_data$variants) || nrow(variant_data$variants) == 0) {
    # empty plot with message
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Click 'Update' to load variants", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # prepare data for plotting
  support_matrix <- variant_data$support
  coverage_matrix <- variant_data$coverage
  variants_df <- variant_data$variants
  
  # always add color information to variants dataframe
  variants_df <- add_variant_colors(variants_df)
  
  # no sampling - use all variants for better plot interaction
    
  # reorder matrices to match configured library order
  lib_indices <- match(library_ids, colnames(support_matrix))
  support_matrix <- support_matrix[, lib_indices, drop = FALSE]
  coverage_matrix <- coverage_matrix[, lib_indices, drop = FALSE]
    
  # get plot mode
  plot_mode <- input$variantPlotMode %||% "frequency"
  
  # calculate data matrix based on plot mode
  if (plot_mode == "variant_support") {
    data_matrix <- support_matrix
    y_label <- "Variant Support"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else if (plot_mode == "variant_coverage") {
    data_matrix <- coverage_matrix
    y_label <- "Variant Coverage"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else { # frequency mode
    data_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
    y_label <- "Frequency"
    y_limits <- c(0, 1)
    y_format <- scales::percent_format()
  }
  
  # convert to long format for plotting
  plot_data <- data.frame()
  for (i in seq_len(nrow(data_matrix))) {
    for (j in seq_len(ncol(data_matrix))) {
      plot_data <- rbind(plot_data, data.frame(
        variant_id = variants_df$variant_id[i],
        variant_type = variants_df$type[i],
        contig = variants_df$contig[i],
        coord = variants_df$coord[i],
        library = library_ids[j],
        value = data_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$value) & !is.na(plot_data$value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # use colors from variants dataframe for plot
  plot_data$variant_color <- variants_df$color[match(plot_data$variant_id, variants_df$variant_id)]
  
  # check if we have selected variants for highlighting
  # always check current table selection state directly
  current_table_selection <- input$variantsTable_rows_selected
  plot_data$is_selected <- FALSE
  
  if (!is.null(current_table_selection) && length(current_table_selection) > 0 && 
      !is.null(variant_data$variants) && nrow(variant_data$variants) > 0) {
    # filter valid row indices
    valid_selection <- current_table_selection[current_table_selection <= nrow(variant_data$variants)]
    if (length(valid_selection) > 0) {
      selected_variant_ids <- variant_data$variants$variant_id[valid_selection]
      plot_data$is_selected <- plot_data$variant_id %in% selected_variant_ids
    }
  }
  
  # set library order as factor using our configured library_ids order
  plot_data$library <- factor(plot_data$library, levels = library_ids)
  
  # create position label and hover text
  plot_data$position_label <- paste0(plot_data$contig, ":", plot_data$coord)
  
  # format value for hover based on plot mode
  if (plot_mode == "frequency") {
    value_text <- paste0(round(plot_data$value * 100, 1), "%")
  } else {
    value_text <- as.character(plot_data$value)
  }
  
  # get descriptions for hover text (truncate if too long)
  plot_data$description_display <- sapply(variants_df$desc[match(plot_data$variant_id, variants_df$variant_id)], function(desc) {
    if (is.na(desc) || nchar(desc) <= 15) {
      return(desc)
    } else {
      # truncate middle for hover (keep first and last 5 characters)
      first_5 <- substr(desc, 1, 5)
      last_5 <- substr(desc, nchar(desc) - 4, nchar(desc))
      return(paste0(first_5, "...", last_5))
    }
  })
  
  # add mutation description to hover if available
  plot_data$mutation_desc <- variants_df$mutation_desc[match(plot_data$variant_id, variants_df$variant_id)]
  
  plot_data$hover_text <- paste0(
    "Variant: ", plot_data$variant_id, "<br>",
    "Position: ", plot_data$position_label, "<br>",
    "Description: ", plot_data$description_display, "<br>",
    ifelse(!is.na(plot_data$mutation_desc) & plot_data$mutation_desc != "", 
           paste0("AA Change: ", plot_data$mutation_desc, "<br>"), ""),
    "Library: ", plot_data$library, "<br>",
    y_label, ": ", value_text, "<br>"
  )
  
  # create plot with conditional styling for selected variant
  non_selected_data <- plot_data[!plot_data$is_selected, ]
  selected_data <- plot_data[plot_data$is_selected, ]
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = variant_color, group = variant_id, text = hover_text))
  
  # draw non-selected variants first with lower alpha
  if (nrow(non_selected_data) > 0) {
    p <- p + 
      ggplot2::geom_line(data = non_selected_data, alpha = 0.3, size = 0.5) +
      ggplot2::geom_point(data = non_selected_data, size = 2, alpha = 0.8)
  }
  
  # draw selected variant on top with higher visibility
  if (nrow(selected_data) > 0) {
    p <- p + 
      ggplot2::geom_line(data = selected_data, alpha = 0.9, size = 2) +
      ggplot2::geom_point(data = selected_data, size = 4, alpha = 1, stroke = 1.5, shape = 21, fill = "white")
  }
  
  p <- p + ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Library",
      y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # add y-axis formatting based on plot mode
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = y_format)
  }
  
  plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(hovermode = "closest", showlegend = FALSE)
})
