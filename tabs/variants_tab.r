# validate and extract tab parameters
tab <- get_tab_by_id("variants")
if (is.null(tab)) {
  stop("variants tab not found during loading")
}

required_params <- c("min_reads", "min_coverage", "min_libraries", "get_aln_f", "library_ids", "max_variants")
missing_params <- required_params[!sapply(required_params, function(p) p %in% names(tab))]
if (length(missing_params) > 0) {
  stop(sprintf("variants tab missing required parameters: %s", paste(missing_params, collapse = ", ")))
}

min_reads <- tab$min_reads
min_coverage <- tab$min_coverage
min_libraries <- tab$min_libraries
get_aln_f <- tab$get_aln_f
library_ids <- tab$library_ids
max_variants_default <- tab$max_variants

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
          h5("Update Variants"),
          actionButton("updateVariantsBtn", "Update", class = "btn-primary", width = "100%"),
          br(), br(),
          selectInput("variantPlotMode", "Plot Mode:", 
                     choices = c("Frequency" = "frequency", 
                                "Variant Support" = "variant_support", 
                                "Variant Coverage" = "variant_coverage"),
                     selected = "frequency", width = "100%"),
          numericInput("variantSpanFilter", "Min Span:", 
                      value = 0.5, min = 0, max = 1, step = 0.1, width = "100%"),
          numericInput("maxVariantsPlot", "Max Variants:", 
                      value = max_variants_default, min = 100, max = 10000, step = 100, width = "100%")
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

# get alignment stores for variant querying using configured get_aln_f function
get_variant_alignment_stores <- function(assembly, contigs, zoom) {
  if (is.null(assembly) || length(contigs) == 0) {
    return(NULL)
  }
  
  # get alignment stores from all configured library_ids
  stores <- list()
  
  for (library_id in library_ids) {
    tryCatch({
      aln <- get_aln_f(assembly, library_id)
      if (!is.null(aln) && inherits(aln, "externalptr")) {
        stores[[library_id]] <- aln
      } else {
        warning(sprintf("get_aln_f returned invalid alignment for library_id '%s'", library_id))
      }
    }, error = function(e) {
      warning(sprintf("error getting alignment for library_id '%s': %s", library_id, e$message))
    })
  }
  
  if (length(stores) == 0) {
    return(NULL)
  }
  
  return(stores)
}

# query variants using alntools
query_variants <- function(assembly, contigs, zoom) {
  # early return if no context
  if (is.null(assembly) || length(contigs) == 0) {
    return(NULL)
  }
  
  # get alignment stores using the new function
  stores <- get_variant_alignment_stores(assembly, contigs, zoom)
  if (is.null(stores) || length(stores) == 0) {
    return(NULL)
  }
  
  # get contigs data and build context for intervals
  contigs_table <- get_contigs(assembly)
  if (is.null(contigs_table)) {
    return(NULL)
  }
  
  # build context for intervals
  cxt <- build_context(contigs, contigs_table, zoom, assembly)
  if (is.null(cxt) || is.null(cxt$intervals) || nrow(cxt$intervals) == 0) {
    return(NULL)
  }
  
  # get alignment filter parameters from UI
  get_filter_param <- function(param_id, default_value, type_converter = identity) {
    param_key <- paste0("align_filter_", param_id)
    if (!param_key %in% names(param_registry)) {
      warning(sprintf("alignment filter parameter '%s' not found, using default: %s", param_id, default_value))
      return(default_value)
    }
    
    param_accessor <- get_param("align_filter", param_id)
    value <- param_accessor()
    return(type_converter(value))
  }
  
  clip_mode <- get_filter_param("clip_mode", "all")
  clip_margin <- get_filter_param("clip_margin", 10L, as.integer)
  min_mutations_percent <- get_filter_param("min_mutations_percent", 0.0, as.numeric)
  max_mutations_percent <- get_filter_param("max_mutations_percent", 10.0, as.numeric)
  min_alignment_length <- get_filter_param("min_alignment_length", 0L, as.integer)
  max_alignment_length <- get_filter_param("max_alignment_length", 0L, as.integer)
  min_indel_length <- get_filter_param("min_indel_length", 3L, as.integer)

  # print filter parameters before query
  cat(sprintf("variant query filters: clip_mode=%s, clip_margin=%d, min_mutations_percent=%.1f%%, max_mutations_percent=%.1f%%, min_alignment_length=%d, max_alignment_length=%d, min_indel_length=%d\n",
              clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, 
              min_alignment_length, max_alignment_length, min_indel_length))
  
  tryCatch({
    # call aln_query_variants with alignment filter parameters from UI
    variants_results <- aln_query_variants(
      store_list = stores,
      intervals_df = cxt$intervals,
      min_variants_variant_support = min_reads,
      min_variants_library_support = min_libraries,
      min_variants_coverage_support = min_coverage,
      clip_mode_str = clip_mode,
      clip_margin = clip_margin,
      min_mutations_percent = min_mutations_percent,
      max_mutations_percent = max_mutations_percent,
      min_alignment_length = min_alignment_length,
      max_alignment_length = max_alignment_length,
      min_indel_length = min_indel_length
    )
    return(variants_results)
  }, error = function(e) {
    cat(sprintf("error querying variants: %s\n", e$message))
    return(NULL)
  })
}

# filter variants by span (frequency range)
filter_variants_by_span <- function(variant_data, min_span) {
  if (is.null(variant_data) || is.null(variant_data$support) || is.null(variant_data$coverage)) {
    return(variant_data)
  }
  
  # calculate frequency for each variant across libraries
  support_matrix <- variant_data$support
  coverage_matrix <- variant_data$coverage
  
  # avoid division by zero
  freq_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
  
  # calculate span (max - min frequency) for each variant
  variant_spans <- apply(freq_matrix, 1, function(row) {
    valid_freqs <- row[!is.na(row) & is.finite(row)]
    if (length(valid_freqs) == 0) return(0)
    max(valid_freqs) - min(valid_freqs)
  })
  
  # filter variants meeting span threshold
  keep_variants <- variant_spans >= min_span
  
  # filter all components
  filtered_data <- list(
    variants = variant_data$variants[keep_variants, ],
    support = variant_data$support[keep_variants, , drop = FALSE],
    coverage = variant_data$coverage[keep_variants, , drop = FALSE],
    library_ids = variant_data$library_ids
  )
  
  return(filtered_data)
}

# ---- Event Handlers ----

# update button handler
observeEvent(input$updateVariantsBtn, {
  
  # query variants
  variant_data <- query_variants(state$assembly, state$contigs, state$zoom)
  
  # store raw data
  state$raw_variant_data <- variant_data
  
  # apply span filter and store filtered data
  if (!is.null(variant_data)) {
    filtered_data <- filter_variants_by_span(variant_data, input$variantSpanFilter %||% 0.5)
    state$filtered_variant_data <- filtered_data
  } else {
    state$filtered_variant_data <- NULL
  }
})

# span filter change handler
observeEvent(input$variantSpanFilter, {
  # reapply filter when span changes
  if (!is.null(state$raw_variant_data)) {
    filtered_data <- filter_variants_by_span(state$raw_variant_data, input$variantSpanFilter)
    state$filtered_variant_data <- filtered_data
  }
})

# ---- Table Renderer ----

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
  
  # round frequency to 3 decimal places
  if ("frequency" %in% names(display_df)) {
    display_df$frequency <- round(display_df$frequency, 3)
  }
  
  # create column name mapping
  col_names <- c(
    "Variant ID" = "variant_id",
    "Contig" = "contig", 
    "Position" = "coord",
    "Type" = "type",
    "Sequence" = "sequence",
    "Libraries" = "library_count",
    "Support" = "total_support",
    "Coverage" = "total_coverage",
    "Frequency" = "frequency"
  )
  
  # filter to existing columns
  available_cols <- names(display_df)
  col_names <- col_names[col_names %in% available_cols]
  
  datatable(
    display_df[, col_names, drop = FALSE],
    rownames = FALSE,
    colnames = names(col_names),
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip"
    ),
    selection = list(mode = "single", target = "row"),
    filter = "top"
  )
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
  
  # apply max variants limit before plotting
  max_variants_limit <- input$maxVariantsPlot %||% max_variants_default
  if (nrow(variants_df) > max_variants_limit) {
    keep_indices <- sample(nrow(variants_df), max_variants_limit)
    variants_df <- variants_df[keep_indices, ]
    support_matrix <- support_matrix[keep_indices, , drop = FALSE]
    coverage_matrix <- coverage_matrix[keep_indices, , drop = FALSE]
  }
    
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
  
  # load alignment color functions
  if (!exists("get_variant_type_colors") || !exists("get_mutation_type_colors")) {
    source("profiles/align/align_utils.r")
  }
  
  # get mutation color mode from alignment profiles
  mutation_color_mode <- "detailed"
  tryCatch({
    profiles <- profiles_get_all()
    for (profile in profiles) {
      if (profile$type == "align" && !is.null(profile$mutation_color_mode)) {
        mutation_color_mode <- profile$mutation_color_mode
        break
      }
    }
  }, error = function(e) {
    # use default
  })
  
  # apply coloring based on mutation mode
  if (mutation_color_mode == "detailed") {
    plot_data$mutation_desc <- variants_df$desc[match(plot_data$variant_id, variants_df$variant_id)]
    unique_descs <- unique(plot_data$mutation_desc)
    type_colors <- get_variant_type_colors(unique_descs)
    names(type_colors) <- unique_descs
    
    missing_descs <- unique_descs[!unique_descs %in% names(type_colors)]
    if (length(missing_descs) > 0) {
      black_colors <- rep("#000000", length(missing_descs))
      names(black_colors) <- missing_descs
      type_colors <- c(type_colors, black_colors)
    }
    plot_data$color_type <- plot_data$mutation_desc
  } else {
    type_mapping <- c("sub" = "substitution", "ins" = "addition", "del" = "deletion")
    plot_data$color_type <- type_mapping[plot_data$variant_type]
    plot_data$color_type[is.na(plot_data$color_type)] <- "clip"
    
    color_types <- unique(plot_data$color_type)
    type_colors <- get_mutation_type_colors(color_types)
    names(type_colors) <- color_types
    type_colors["clip"] <- "#000000"
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
  
  # truncate sequences for hover text
  plot_data$sequence_display <- sapply(variants_df$sequence[match(plot_data$variant_id, variants_df$variant_id)], function(seq) {
    if (is.na(seq) || nchar(seq) <= 10) {
      return(seq)
    } else {
      first_5 <- substr(seq, 1, 5)
      last_5 <- substr(seq, nchar(seq) - 4, nchar(seq))
      return(paste0(first_5, "...", last_5))
    }
  })
  
  plot_data$hover_text <- paste0(
    "Variant: ", plot_data$variant_id, "<br>",
    "Position: ", plot_data$position_label, "<br>",
    "Sequence: ", plot_data$sequence_display, "<br>",
    "Library: ", plot_data$library, "<br>",
    y_label, ": ", value_text, "<br>",
    "Type: ", plot_data$color_type
  )
  
  # create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = color_type, group = variant_id, text = hover_text)) +
    ggplot2::geom_line(alpha = 0.3, size = 0.5) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = type_colors, na.value = "#000000") +
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
