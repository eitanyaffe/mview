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
    dat <- contigs
    index <- which(names(dat) == "cid")
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row"),
        options = get.highlight.options(state$contigs, index)
    )
})

output$genomeTable <- renderDT({
    datatable(genomes, selection = list(mode = "multiple", target = "row"))
})

output$mapTable <- renderDT({
    dat <- contig_map
    index <- which(names(dat) == "cid")
    datatable(
        dat,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row"),
        options = get.highlight.options(state$contigs, index)
    )
})

output$selectedTable <- renderDT({
    selected_df <- contigs[contigs$cid %in% state$contigs, ]
    selected_df$gid_list <- sapply(selected_df$cid, function(cid) {
        paste(contig_map$gid[contig_map$cid == cid], collapse = ", ")
    })
    selected_df <- selected_df[, c("cid", "length", "gid_list")]
    datatable(selected_df,
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row")
    )
})
