# ---- Button Event Observers ----

selected_cids_from_table <- function(
    state_cids,
    full_contig_table,
    selected_rows) {
  visible_df <- full_contig_table[full_contig_table$cid %in% state_cids, ]
  visible_df <- visible_df[match(state_cids, visible_df$cid), ]
  visible_df$cid[selected_rows]
}

observeEvent(input$addContigsBtn, {
  rows <- input$contigTable_rows_selected
  if (length(rows) > 0) {
    contigs_data <- get_contigs(state$assembly)
    new_cids <- contigs_data$cid[rows]
    new_cids <- new_cids[!new_cids %in% state$contigs]
    if (length(new_cids) > 0) {
      states_module_output$push_state()
      state$contigs <- c(state$contigs, new_cids)
      addLog(paste("Add:", paste(new_cids, collapse = ",")))
    }
  }
})

observeEvent(input$addGenomesBtn, {
  rows <- input$genomeTable_rows_selected
  if (length(rows) > 0) {
    genomes_data <- get_genomes(state$assembly)
    contig_map_data <- get_contig_map(state$assembly)
    selected_gids <- genomes_data$gid[rows]
    new_cids <- contig_map_data$cid[contig_map_data$gid %in% selected_gids]
    new_cids <- new_cids[!new_cids %in% state$contigs]
    if (length(new_cids) > 0) {
      states_module_output$push_state()
      state$contigs <- c(state$contigs, new_cids)
      addLog(paste("From genome:", paste(selected_gids, collapse = ",")))
    }
  }
})

observeEvent(input$removeContigsBtn, {
  rows <- input$selectedTable_rows_selected
  contigs_data <- get_contigs(state$assembly)
  removed_cids <- selected_cids_from_table(state$contigs, contigs_data, rows)
  if (length(removed_cids) > 0) {
    states_module_output$push_state()
    state$contigs <- setdiff(state$contigs, removed_cids)
    addLog(paste("Del:", paste(removed_cids, collapse = ",")))
  }
})
