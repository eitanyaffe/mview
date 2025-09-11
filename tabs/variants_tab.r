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
          h5("Update Variants"),
          actionButton("updateVariantsBtn", "Update", class = "btn-primary", width = "100%"),
          br(), br(),
          selectInput("variantPlotMode", "Plot Mode:", 
                     choices = c("Frequency" = "frequency", 
                                "Variant Support" = "variant_support", 
                                "Variant Coverage" = "variant_coverage"),
                     selected = "frequency", width = "100%"),
          numericInput("variantSpanFilter", "Min Span:", 
                      value = 0.5, min = 0, max = 1, step = 0.1, width = "100%")
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

# ---- Gene Object Creation ----

# create genes object for variant annotation
create_genes_object <- function(assembly, use_genes, get_gene_table_f, get_fasta_f, codon_table_path) {
  genes_object <- NULL

  if (use_genes) {
    if (!is.function(get_gene_table_f) || !is.function(get_fasta_f)) {
      stop("get_gene_table_f and get_fasta_f must be functions")
    }

    # create genes object with caching
    genes_object <- cache(paste0(assembly, "_genes_object"), {
      # get gene table data directly
      gene_table_data <- get_gene_table_f(assembly)
      if (!is.null(gene_table_data) && nrow(gene_table_data) > 0) {
        # get reference sequences directly
        reference_sequences_data <- get_fasta_f(assembly)
        
        genes_obj <- aln_genes(
          gene_table = gene_table_data,
          reference_sequences = reference_sequences_data,
          codon_table_path = codon_table_path
        )
        
        cat(sprintf("using gene annotation: %d genes loaded, %d reference sequences\n", 
                   nrow(gene_table_data), length(reference_sequences_data)))
        
        genes_obj
      } else {
        cat("gene table function returned no data, proceeding without gene annotation\n")
        NULL
      }
    })
  }

  return(genes_object)
}

# enhance variants with gene annotation data
enhance_variants_with_gene_annotation <- function(variants_results, use_genes) {
  # if genes were used, enhance variants data with gene annotation from alntools output
  if (use_genes) {
    variants_df <- variants_results$variants
    
    # only add gene annotation columns if there are variants
    if (nrow(variants_df) > 0) {
      # check if gene annotation data is available
      if (!is.null(variants_results$genic) && nrow(variants_results$genic) > 0) {
        genic_data <- variants_results$genic
        # match variant_id with row_id in genic data
        ix <- match(variants_df$variant_id, genic_data$row_id)
        
        # add gene columns using vectorized operations
        variants_df$gene_desc <- ifelse(!is.na(ix), genic_data$gene_desc[ix], "none")
        variants_df$mutation_desc <- ifelse(!is.na(ix), genic_data$mutation_desc[ix], "")
        
        cat(sprintf("enhanced %d variants with gene annotation data\n", sum(!is.na(ix))))
      } else {
        # no gene data available, add default columns
        variants_df$gene_desc <- "none"
        variants_df$mutation_desc <- ""
        cat("no gene annotation data available\n")
      }
      
      # update the variants results
      variants_results$variants <- variants_df
    } else {
      cat("no variants found, skipping gene annotation\n")
    }
  }
  
  return(variants_results)
}

# ---- Variant Data Management ----

# add color information to variants dataframe
add_variant_colors <- function(variants_df) {
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    return(variants_df)
  }
  
  # load alignment color functions if not available
  if (!exists("get_variant_type_colors") || !exists("get_mutation_type_colors")) {
    source("profiles/align/align_utils.r")
  }
  
  # get mutation color mode from alignment profiles
  mutation_color_mode <- profile_get_param("mutation_color_mode", "align", "detailed")
  
  # apply coloring based on mutation mode and add colors to variants dataframe
  if (mutation_color_mode == "detailed") {
    unique_descs <- unique(variants_df$desc)
    type_colors <- get_variant_type_colors(unique_descs)
    names(type_colors) <- unique_descs
    
    missing_descs <- unique_descs[!unique_descs %in% names(type_colors)]
    if (length(missing_descs) > 0) {
      black_colors <- rep("#000000", length(missing_descs))
      names(black_colors) <- missing_descs
      type_colors <- c(type_colors, black_colors)
    }
    
    # add color column to variants dataframe
    variants_df$color <- type_colors[variants_df$desc]
  } else {
    type_mapping <- c("sub" = "substitution", "ins" = "addition", "del" = "deletion")
    color_type <- type_mapping[variants_df$type]
    color_type[is.na(color_type)] <- "clip"
    
    color_types <- unique(color_type)
    type_colors <- get_mutation_type_colors(color_types)
    names(type_colors) <- color_types
    type_colors["clip"] <- "#000000"
    
    # add color column to variants dataframe
    variants_df$color <- type_colors[color_type]
  }
  
  return(variants_df)
}

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
  
  # prepare gene annotation data if genes are enabled
  genes_object <- create_genes_object(assembly, use_genes, get_gene_table_f, get_fasta_f, codon_table_path)

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
    min_indel_length = min_indel_length,
    genes = genes_object
  )
  
  # enhance variants with gene annotation if genes were used
  variants_results <- enhance_variants_with_gene_annotation(variants_results, use_genes)
  
  return(variants_results)
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
    
    # refresh plots to remove highlighting
    if (exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
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
      
      # refresh plots to show highlighting
      if (exists("refresh_trigger")) {
        current_val <- refresh_trigger()
        refresh_trigger(current_val + 1)
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
    
    # refresh plots to show updated variants
    if (exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
  } else {
    state$filtered_variant_data <- NULL
    state$variants <- NULL
    cache_set("variants.current", NULL)
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
    
    # refresh plots to show updated variants
    if (exists("refresh_trigger")) {
      current_val <- refresh_trigger()
      refresh_trigger(current_val + 1)
    }
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
