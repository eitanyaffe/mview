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
    new_cids <- contigs$cid[rows]
    new_cids <- new_cids[!new_cids %in% state$contigs]
    if (length(new_cids) > 0) {
      push_state_trigger(isolate(push_state_trigger()) + 1)
      state$contigs <- c(state$contigs, new_cids)
      addLog(paste("Add:", paste(new_cids, collapse = ",")))
    }
  }
})

observeEvent(input$addGenomesBtn, {
  rows <- input$genomeTable_rows_selected
  if (length(rows) > 0) {
    selected_gids <- genomes$gid[rows]
    new_cids <- contig_map$cid[contig_map$gid %in% selected_gids]
    new_cids <- new_cids[!new_cids %in% state$contigs]
    if (length(new_cids) > 0) {
      push_state_trigger(isolate(push_state_trigger()) + 1)
      state$contigs <- c(state$contigs, new_cids)
      addLog(paste("From genome:", paste(selected_gids, collapse = ",")))
    }
  }
})

observeEvent(input$removeContigsBtn, {
  rows <- input$selectedTable_rows_selected
  removed_cids <- selected_cids_from_table(state$contigs, contigs, rows)
  if (length(removed_cids) > 0) {
    push_state_trigger(isolate(push_state_trigger()) + 1)
    state$contigs <- setdiff(state$contigs, removed_cids)
    addLog(paste("Del:", paste(removed_cids, collapse = ",")))
  }
})

observeEvent(input$helpBtn, {
  showModal(modalDialog(
    title = "Keyboard Shortcuts",
    "Double-click window to zoom in",
    "Double-click without window to zoom out",
    easyClose = TRUE
  ))
})
