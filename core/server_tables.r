# ---- DataTable Renderers ----

get.highlight.options <- function(contigs, index, enable_highlighting = TRUE) {
    if (!enable_highlighting) {
        # return empty options when highlighting is disabled
        return(list())
    }
    
    str <- sprintf(
        "function(row, data) {
       const selected = %s;
       if (selected.includes(data[%s])) {
         $(row).css('background', '#ccffcc');
       } else {
         $(row).css('background', '');
       } }",
        jsonlite::toJSON(contigs, auto_unbox = TRUE), index - 1
    )

    list(rowCallback = JS(str))
}

output$contigTable <- renderDT({
    
    dat <- get_contigs(state$assembly)
    
    # apply length filtering based on options (convert kb to bp)
    if (!is.null(dat) && "length" %in% names(dat)) {
        min_length_kb <- if (is.null(input$min_contig_length)) 0 else input$min_contig_length
        max_length_kb <- input$max_contig_length
        
        # convert kb to bp for comparison with contig lengths
        min_length_bp <- min_length_kb * 1000
        max_length_bp <- if (!is.null(max_length_kb) && !is.na(max_length_kb)) max_length_kb * 1000 else NULL
        
        # filter by minimum length
        dat <- dat[dat$length >= min_length_bp, ]
        
        # filter by maximum length if specified
        if (!is.null(max_length_bp)) {
            dat <- dat[dat$length <= max_length_bp, ]
        }
    }
    
    index <- which(names(dat) == "contig")
    # check if highlighting is enabled in options
    enable_highlighting <- if (is.null(input$enable_contig_highlighting)) TRUE else input$enable_contig_highlighting
    
    # determine selection mode based on checkbox state
    selection_mode <- if (!is.null(input$allowMultipleContigsChk) && !input$allowMultipleContigsChk) "single" else "multiple"
    
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = selection_mode, target = "row"),
        options = get.highlight.options(get_state_contigs(), index, enable_highlighting)
    )
})

output$genomeTable <- renderDT({
    # depend on field selection trigger to ensure updates
    genome_fields_trigger()
    
    # determine selection mode based on checkbox state
    selection_mode <- if (!is.null(input$allowMultipleGenomesChk) && !input$allowMultipleGenomesChk) "single" else "multiple"
    
    # get genome data
    genomes_df <- get_genomes(state$assembly)
    
    # filter fields based on selection
    selected_fields <- get_selected_genome_fields(state$assembly)
    genomes_df <- filter_genome_table_fields(genomes_df, selected_fields)
    
    # format column names for display
    col_names <- get_genome_column_names(names(genomes_df))
    
    datatable(
        genomes_df, 
        selection = list(mode = selection_mode, target = "row"),
        colnames = col_names,
        rownames = FALSE,
        escape = FALSE,
        options = list(
          columnDefs = list(
            list(className = "dt-left", targets = "_all")
          )
        )
    )
})

output$segmentTable <- renderDT({
    dat <- get_segments(state$assembly)
    index <- which(names(dat) == "segment")
    enable_highlighting <- if (is.null(input$enable_contig_highlighting)) TRUE else input$enable_contig_highlighting
    current_segments <- get_state_segments()$segment
    selection_mode <- if (!is.null(input$allowMultipleSegmentsChk) && !input$allowMultipleSegmentsChk) "single" else "multiple"
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = selection_mode, target = "row"),
        options = get.highlight.options(current_segments, index, enable_highlighting)
    )
})

output$csegmentTable <- renderDT({
  cseg_table <- get_csegments(state$assembly)
  mapping <- get_cluster_mapping(state$assembly)
  
  if (is.null(cseg_table)) {
    return(datatable(
      data.frame(Message = "No csegment data available"),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # add segments column
  if (!is.null(mapping)) {
    mapping_split <- split(mapping$segment, mapping$csegment)
    cseg_table$segments <- sapply(cseg_table$csegment, function(cs) {
      segs <- mapping_split[[as.character(cs)]]
      if (is.null(segs)) "" else paste(segs, collapse = ", ")
    })
  } else {
    cseg_table$segments <- ""
  }
  
  index <- which(names(cseg_table) == "csegment")
  enable_highlighting <- if (is.null(input$enable_contig_highlighting)) TRUE else input$enable_contig_highlighting
  selection_mode <- if (!is.null(input$allowMultipleCsegmentsChk) && !input$allowMultipleCsegmentsChk) "single" else "multiple"
  
  datatable(
    cseg_table,
    rownames = FALSE,
    selection = list(mode = selection_mode, target = "row"),
    options = get.highlight.options(character(), index, enable_highlighting)
  )
})

output$mapTable <- renderDT({
    dat <- get_segment_map(state$assembly)
    index <- which(names(dat) == "segment")
    # check if highlighting is enabled in options
    enable_highlighting <- if (is.null(input$enable_contig_highlighting)) TRUE else input$enable_contig_highlighting
    current_segments <- get_state_segments()$segment
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row"),
        options = get.highlight.options(current_segments, index, enable_highlighting)
    )
})

