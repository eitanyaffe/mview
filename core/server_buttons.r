# ---- Button Event Observers ----

# helper to add segments for selected contigs
add_segments_for_contigs <- function(contigs_to_add, assembly) {
  if (length(contigs_to_add) == 0) return(data.frame())
  
  segments <- get_segments(assembly)
  new_segments <- segments[segments$contig %in% contigs_to_add, ]
  
  # filter out segments already in state
  current_segments <- get_state_segments()
  if (nrow(current_segments) > 0) {
    new_segments <- new_segments[!new_segments$segment %in% current_segments$segment, ]
  }
  
  return(new_segments)
}

# helper to remove segments for contigs
remove_segments_for_contigs <- function(contigs_to_remove) {
  current_segments <- get_state_segments()
  if (nrow(current_segments) == 0) return(current_segments)
  
  remaining_segments <- current_segments[!current_segments$contig %in% contigs_to_remove, ]
  return(remaining_segments)
}

observeEvent(input$addContigsBtn, {
  rows <- input$contigTable_rows_selected
  if (length(rows) > 0) {
    contigs_data <- get_contigs(state$assembly)
    selected_contigs <- contigs_data$contig[rows]
    current_contigs <- get_state_contigs()
    new_contigs <- selected_contigs[!selected_contigs %in% current_contigs]
    if (length(new_contigs) > 0) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      
      # add segments for selected contigs
      new_segments <- add_segments_for_contigs(new_contigs, state$assembly)
      if (nrow(new_segments) > 0) {
        state$segments <- rbind(get_state_segments(), new_segments)
      }
      
      # reset zoom to see full range of all segments
      state$zoom <- NULL
    }
  }
})

observeEvent(input$gotoContigsBtn, {
  rows <- input$contigTable_rows_selected
  if (length(rows) > 0) {
    contigs_data <- get_contigs(state$assembly)
    selected_contigs <- contigs_data$contig[rows]
    current_contigs <- get_state_contigs()
    if (!identical(selected_contigs, current_contigs)) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      
      # get segments for selected contigs
      segments <- get_segments(state$assembly)
      new_segments <- segments[segments$contig %in% selected_contigs, ]
      state$segments <- new_segments
      
      # reset zoom to see full range of new segments
      state$zoom <- NULL
    }
  }
})

observeEvent(input$addGenomesBtn, {
  rows <- input$genomeTable_rows_selected
  if (length(rows) > 0) {
    genomes_data <- get_genomes(state$assembly)
    segment_map_data <- get_segment_map(state$assembly)
    selected_gids <- genomes_data$gid[rows]
    selected_segment_ids <- segment_map_data$segment[segment_map_data$gid %in% selected_gids]
    current_segment_ids <- get_state_segments()$segment
    new_segment_ids <- selected_segment_ids[!selected_segment_ids %in% current_segment_ids]
    if (length(new_segment_ids) > 0) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      
      # add selected segments
      all_segments <- get_segments(state$assembly)
      new_segments <- all_segments[all_segments$segment %in% new_segment_ids, ]
      if (nrow(new_segments) > 0) {
        state$segments <- rbind(get_state_segments(), new_segments)
      }
      
      # reset zoom to see full range of all segments
      state$zoom <- NULL
    }
  }
})

observeEvent(input$gotoGenomesBtn, {
  rows <- input$genomeTable_rows_selected
  if (length(rows) > 0) {
    genomes_data <- get_genomes(state$assembly)
    segment_map_data <- get_segment_map(state$assembly)
    selected_gids <- genomes_data$gid[rows]
    selected_segment_ids <- segment_map_data$segment[segment_map_data$gid %in% selected_gids]
    current_segment_ids <- get_state_segments()$segment
    if (!identical(sort(selected_segment_ids), sort(current_segment_ids))) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      
      # get selected segments
      all_segments <- get_segments(state$assembly)
      new_segments <- all_segments[all_segments$segment %in% selected_segment_ids, ]
      state$segments <- new_segments
      
      # reset zoom to see full range of new segments
      state$zoom <- NULL
    }
  }
})

observeEvent(input$removeContigsFromListBtn, {
  rows <- input$contigTable_rows_selected
  if (length(rows) > 0) {
    contigs_data <- get_contigs(state$assembly)
    selected_contigs <- contigs_data$contig[rows]
    current_contigs <- get_state_contigs()
    contigs_to_remove <- selected_contigs[selected_contigs %in% current_contigs]
    if (length(contigs_to_remove) > 0) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      
      # remove segments for contigs
      state$segments <- remove_segments_for_contigs(contigs_to_remove)
      
      # reset zoom to see full range of remaining segments
      if (nrow(get_state_segments()) > 0) {
        state$zoom <- NULL
      }
    }
  }
})

observeEvent(input$removeGenomesFromListBtn, {
  rows <- input$genomeTable_rows_selected
  if (length(rows) > 0) {
    genomes_data <- get_genomes(state$assembly)
    segment_map_data <- get_segment_map(state$assembly)
    selected_gids <- genomes_data$gid[rows]
    segments_to_remove <- segment_map_data$segment[segment_map_data$gid %in% selected_gids]
    current_segment_ids <- get_state_segments()$segment
    segments_to_remove <- segments_to_remove[segments_to_remove %in% current_segment_ids]
    if (length(segments_to_remove) > 0) {
      # push current region to undo before changing
      regions_module_output$push_undo_state()
      
      # remove selected segments
      current_segments <- get_state_segments()
      state$segments <- current_segments[!current_segments$segment %in% segments_to_remove, ]
      
      # reset zoom to see full range of remaining segments
      if (nrow(get_state_segments()) > 0) {
        state$zoom <- NULL
      }
    }
  }
})

# create proxies for tables to manipulate selections
contigTableProxy <- DT::dataTableProxy("contigTable")
genomeTableProxy <- DT::dataTableProxy("genomeTable")

# persist checkbox states in cache
observeEvent(input$allowMultipleContigsChk, {
  cache_set("allow_multiple_contigs", input$allowMultipleContigsChk)
})

observeEvent(input$allowMultipleGenomesChk, {
  cache_set("allow_multiple_genomes", input$allowMultipleGenomesChk)
})


observeEvent(input$clearContigSelectionBtn, {
  # clear the selected rows in the contig table
  DT::selectRows(contigTableProxy, NULL)
})

observeEvent(input$clearGenomeSelectionBtn, {
  # clear the selected rows in the genome table
  DT::selectRows(genomeTableProxy, NULL)
})


