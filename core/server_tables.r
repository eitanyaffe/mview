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
