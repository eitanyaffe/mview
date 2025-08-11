# ---- Gene Table Functions ----

# set the tab panel UI
set_tab_panel_f(function() {
  tabPanel(
    "Genes",
    div(
      style = "margin-bottom: 10px;",
      actionButton("showGeneDetailsBtn", "Show Details"),
      actionButton("zoomToGeneBtn", "Zoom to Gene")
    ),
    DTOutput("genesTable")
  )
})

# function to get genes for the current context (assembly, contigs, zoom)
get_genes_for_context <- function(assembly, contigs, zoom) {
  # early return if no assembly or contigs
  if (is.null(assembly) || length(contigs) == 0) {
    return(NULL)
  }

  # get contigs data and build context
  contigs_table <- get_contigs(assembly)
  if (is.null(contigs_table)) {
    return(NULL)
  }

  tab <- get_tab_by_id("genes")
  if (is.null(tab$get_genes_f)) {
    stop("get_genes_f is not defined for tab ", tab$tab_id)
  }
  # get gene data using the tab's gene function if provided
  cxt <- build_context(contigs, contigs_table, zoom, assembly)
  if (is.null(cxt)) return(NULL)
  genes <- tab$get_genes_f(cxt)

  if (is.null(genes) || nrow(genes) == 0) {
    return(NULL)
  }

  # ensure required columns exist
  required_cols <- c("contig", "start", "end")
  missing_cols <- required_cols[!required_cols %in% names(genes)]
  if (length(missing_cols) > 0) {
    warning(paste("missing required columns in genes table:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }

  # process genes for filtering
  genes$contig <- as.character(genes$contig)
  genes$start <- as.numeric(genes$start)
  genes$end <- as.numeric(genes$end)

  # build context for filtering
  cxt <- build_context(contigs, contigs_table, zoom, assembly)
  if (is.null(cxt)) {
    return(NULL)
  }

  # filter to visible range using the context
  filtered_genes <- filter_segments(genes, cxt, cxt$mapper$xlim)

  # verify we have the expected columns after filtering
  if (!is.null(filtered_genes)) {
    expected_cols <- c("gstart", "gend")
    missing_after_filter <- expected_cols[!expected_cols %in% names(filtered_genes)]
    if (length(missing_after_filter) > 0) {
      warning(paste("missing columns after filtering:", paste(missing_after_filter, collapse = ", ")))
      return(NULL)
    }
  }

  return(filtered_genes)
}

# ---- DataTable Renderer ----

output$genesTable <- renderDT({
  genes_df <- get_genes_for_context(state$assembly, state$contigs, state$zoom)

  if (is.null(genes_df) || nrow(genes_df) == 0) {
    return(datatable(
      data.frame(Message = "No genes found in the current view"),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }

  # Explicitly verify gene data structure before proceeding
  required_gene_cols <- c("gene", "contig", "start", "end", "gstart", "gend")
  missing_gene_cols <- required_gene_cols[!required_gene_cols %in% names(genes_df)]
  if (length(missing_gene_cols) > 0) {
    cat("ERROR: Missing required columns in genes data:", paste(missing_gene_cols, collapse = ", "), "\n")
    return(datatable(
      data.frame(Error = paste("Data error: Missing required columns:", paste(missing_gene_cols, collapse = ", "))),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }

  # Select and order columns for display
  display_cols <- c("gene", "contig", "start", "end", "strand")

  # Add additional columns if they exist
  additional_cols <- c("prot_desc", "tax", "identity", "coverage")
  for (col in additional_cols) {
    if (col %in% names(genes_df)) {
      display_cols <- c(display_cols, col)
    }
  }

  # Ensure all display columns exist in the data
  display_cols <- display_cols[display_cols %in% names(genes_df)]

  # Format the numeric columns
  formatted_df <- genes_df[, display_cols, drop = FALSE]

  # Round numeric values for better display
  numeric_cols <- c("identity", "coverage", "evalue", "bitscore")
  for (col in numeric_cols) {
    if (col %in% names(formatted_df) && is.numeric(formatted_df[[col]])) {
      if (col == "evalue") {
        # Format e-values in scientific notation
        formatted_df[[col]] <- formatC(formatted_df[[col]], format = "e", digits = 2)
      } else {
        # Round other numeric values to 2 decimal places
        formatted_df[[col]] <- round(formatted_df[[col]], 2)
      }
    }
  }

  # Create a better column name map for display
  available_cols <- names(formatted_df)
  col_names <- c()

  # Corrected mapping: 'New Display Name' = 'Original Column Name'
  if ("gene" %in% available_cols) col_names["Gene"] <- "gene"
  if ("contig" %in% available_cols) col_names["Contig"] <- "contig"
  if ("start" %in% available_cols) col_names["Start"] <- "start"
  if ("end" %in% available_cols) col_names["End"] <- "end"
  if ("strand" %in% available_cols) col_names["Strand"] <- "strand"
  if ("prot_desc" %in% available_cols) col_names["Description"] <- "prot_desc"
  if ("tax" %in% available_cols) col_names["Taxonomy"] <- "tax"
  if ("identity" %in% available_cols) col_names["Identity"] <- "identity"
  if ("coverage" %in% available_cols) col_names["Coverage"] <- "coverage"

  # Create DT options
  dt_options <- list(
    pageLength = 15,
    lengthMenu = c(5, 10, 15, 25, 50),
    scrollX = TRUE,
    dom = "lftip"
  )

  # Create the datatable WITHOUT specifying the escape parameter at all
  datatable(
    formatted_df,
    rownames = FALSE,
    colnames = col_names,
    options = dt_options,
    selection = list(mode = "single", target = "row"),
    filter = "top"
  )
})

# ---- Button Event Observers ----

observeEvent(input$showGeneDetailsBtn, {
  selected_row <- input$genesTable_rows_selected
  if (length(selected_row) > 0) {
    genes_df <- get_genes_for_context(state$assembly, state$contigs, state$zoom)
    if (!is.null(genes_df) && nrow(genes_df) > 0) {
      selected_gene <- genes_df[selected_row, ]

      # Create a formatted string with gene details
      gene_details <- paste0(
        "Gene: ", selected_gene$gene, "\n",
        "Contig: ", selected_gene$contig, "\n",
        "Position: ", selected_gene$start, " - ", selected_gene$end, " (", selected_gene$strand, ")\n",
        if ("prot_desc" %in% names(selected_gene)) paste0("Description: ", selected_gene$prot_desc, "\n") else "",
        if ("tax" %in% names(selected_gene)) paste0("Taxonomy: ", selected_gene$tax, "\n") else ""
      )

      # Show the details in a modal dialog
      showModal(modalDialog(
        title = paste("Gene Details:", selected_gene$gene),
        tags$pre(gene_details),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    }
  }
})

observeEvent(input$zoomToGeneBtn, {
  selected_row <- input$genesTable_rows_selected
  if (length(selected_row) > 0) {
    genes_df <- get_genes_for_context(state$assembly, state$contigs, state$zoom)
    if (!is.null(genes_df) && nrow(genes_df) > 0) {
      selected_gene <- genes_df[selected_row, ]

      # Save current state before changing zoom
      states_module_output$push_state()

      # Set zoom to the gene with some padding
      padding <- max(500, (selected_gene$gend - selected_gene$gstart) * 0.2)
      state$zoom <- c(
        selected_gene$gstart - padding,
        selected_gene$gend + padding
      )      
    }
  }
})
