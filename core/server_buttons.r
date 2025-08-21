# ---- Button Event Observers ----

selected_contigs_from_table <- function(
    state_contigs,
    full_contig_table,
    selected_rows) {
  visible_df <- full_contig_table[full_contig_table$contig %in% state_contigs, ]
  visible_df <- visible_df[match(state_contigs, visible_df$contig), ]
  visible_df$contig[selected_rows]
}

observeEvent(input$addContigsBtn, {
  rows <- input$contigTable_rows_selected
  if (length(rows) > 0) {
    contigs_data <- get_contigs(state$assembly)
    new_contigs <- contigs_data$contig[rows]
    new_contigs <- new_contigs[!new_contigs %in% state$contigs]
    if (length(new_contigs) > 0) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      state$contigs <- c(state$contigs, new_contigs)
      # reset zoom to see full range of all contigs
      state$zoom <- NULL
    }
  }
})

observeEvent(input$setContigsBtn, {
  rows <- input$contigTable_rows_selected
  if (length(rows) > 0) {
    contigs_data <- get_contigs(state$assembly)
    selected_contigs <- contigs_data$contig[rows]
    if (!identical(selected_contigs, state$contigs)) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      state$contigs <- selected_contigs
      # reset zoom to see full range of new contigs
      state$zoom <- NULL
    }
  }
})

observeEvent(input$addGenomesBtn, {
  rows <- input$genomeTable_rows_selected
  if (length(rows) > 0) {
    genomes_data <- get_genomes(state$assembly)
    contig_map_data <- get_contig_map(state$assembly)
    selected_gids <- genomes_data$gid[rows]
    new_contigs <- contig_map_data$contig[contig_map_data$gid %in% selected_gids]
    new_contigs <- new_contigs[!new_contigs %in% state$contigs]
    if (length(new_contigs) > 0) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      state$contigs <- c(state$contigs, new_contigs)
      # reset zoom to see full range of all contigs
      state$zoom <- NULL
    }
  }
})

observeEvent(input$removeContigsBtn, {
  rows <- input$selectedTable_rows_selected
  contigs_data <- get_contigs(state$assembly)
  removed_contigs <- selected_contigs_from_table(state$contigs, contigs_data, rows)
  if (length(removed_contigs) > 0) {
    # push current state to undo before changing
    regions_module_output$push_undo_state()
    state$contigs <- setdiff(state$contigs, removed_contigs)
  }
})
