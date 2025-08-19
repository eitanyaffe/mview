# ---- DataTable Renderers ----

get.highlight.options <- function(contigs, index) {
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
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row"),
        options = get.highlight.options(state$contigs, index)
    )
})

output$genomeTable <- renderDT({
    datatable(get_genomes(state$assembly), selection = list(mode = "multiple", target = "row"))
})

output$mapTable <- renderDT({
    dat <- get_contig_map(state$assembly)
    index <- which(names(dat) == "contig")
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row"),
        options = get.highlight.options(state$contigs, index)
    )
})

output$selectedTable <- renderDT({
    contigs_data <- get_contigs(state$assembly)
    contig_map_data <- get_contig_map(state$assembly)
    selected_df <- contigs_data[contigs_data$contig %in% state$contigs, ]
    selected_df$gid_list <- sapply(selected_df$contig, function(contig) {
        paste(contig_map_data$gid[contig_map_data$contig == contig], collapse = ", ")
    })
    selected_df <- selected_df[, c("contig", "length", "gid_list")]
    datatable(selected_df,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row")
    )
})
