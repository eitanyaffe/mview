# core/genome_fields.r
# Genome table field selection functionality

# reactive trigger to force table update when field selection changes
genome_fields_trigger <- reactiveVal(0)

# get available fields from genome table
get_genome_fields <- function(assembly) {
  genomes <- get_genomes(assembly)
  if (is.null(genomes) || nrow(genomes) == 0) {
    return(character())
  }
  # exclude total_length since we already have length
  # exclude var_sweep_density
  fields <- names(genomes)
  fields <- fields[fields != "total_length"]
  fields <- fields[fields != "var_sweep_density"]
  return(fields)
}

# get selected fields from cache (defaults to all fields if not set)
get_selected_genome_fields <- function(assembly) {
  cache_key <- "genome_fields_selected"
  cached <- cache_get_if_exists(cache_key, NULL)
  
  # if no cache, return all fields
  if (is.null(cached)) {
    return(NULL)  # NULL means show all fields
  }
  
  # get available fields
  available_fields <- get_genome_fields(assembly)
  if (length(available_fields) == 0) {
    return(NULL)
  }
  
  # filter cached selection to only include fields that exist
  selected <- intersect(cached, available_fields)
  
  # if no valid fields selected, return NULL to show all
  if (length(selected) == 0) {
    return(NULL)
  }
  
  return(selected)
}

# set selected fields in cache
set_selected_genome_fields <- function(fields) {
  cache_key <- "genome_fields_selected"
  cache_set(cache_key, fields)
}

# format column name for display
format_genome_column_name <- function(col_name) {
  # handle special cases first
  if (col_name == "total_element_length") {
    return("element<br>length")
  }
  if (col_name == "total_length") {
    return("length")
  }
  
  # default: replace underscores with <br>
  gsub("_", "<br>", col_name)
}

# create column name vector for display (in same order as columns)
get_genome_column_names <- function(column_names) {
  display_names <- sapply(column_names, format_genome_column_name, USE.NAMES = FALSE)
  return(display_names)
}

# filter genome table to show only selected fields
filter_genome_table_fields <- function(genomes_df, selected_fields) {
  if (is.null(genomes_df) || nrow(genomes_df) == 0) {
    return(genomes_df)
  }
  
  # exclude total_length since we already have length
  # exclude var_sweep_density
  available_fields <- names(genomes_df)
  available_fields <- available_fields[available_fields != "total_length"]
  available_fields <- available_fields[available_fields != "var_sweep_density"]
  
  # if no selection or NULL, show all fields (except excluded ones)
  if (is.null(selected_fields) || length(selected_fields) == 0) {
    return(genomes_df[, available_fields, drop = FALSE])
  }
  
  # ensure gid is always included
  required_fields <- c("gid")
  fields_to_show <- unique(c(required_fields, selected_fields))
  
  # filter to only show selected fields that exist (and exclude excluded ones)
  fields_to_show <- fields_to_show[fields_to_show %in% available_fields]
  
  if (length(fields_to_show) == 0) {
    return(genomes_df[, available_fields, drop = FALSE])
  }
  
  return(genomes_df[, fields_to_show, drop = FALSE])
}

# observer for field selection dialog
observeEvent(input$genomeFieldsBtn, {
  assembly <- state$assembly
  if (is.null(assembly)) {
    showNotification("No assembly selected", type = "error")
    return()
  }
  
  # get available fields
  available_fields <- get_genome_fields(assembly)
  if (length(available_fields) == 0) {
    showNotification("No genome data available", type = "error")
    return()
  }
  
  # get currently selected fields (or all if none selected)
  current_selected <- get_selected_genome_fields(assembly)
  if (is.null(current_selected)) {
    current_selected <- available_fields
  }
  
  # create checkboxes for each field
  field_checkboxes <- lapply(available_fields, function(field) {
    checkboxInput(
      inputId = paste0("genome_field_", gsub("[^A-Za-z0-9]", "_", field)),
      label = field,
      value = field %in% current_selected
    )
  })
  
  # show modal dialog
  showModal(modalDialog(
    title = "Select Genome Table Fields",
    div(
      style = "max-height: 400px; overflow-y: auto;",
      do.call(div, field_checkboxes)
    ),
    footer = tagList(
      actionButton("genomeFieldsSelectAll", "Select All", class = "btn-sm"),
      actionButton("genomeFieldsDeselectAll", "Deselect All", class = "btn-sm"),
      modalButton("Cancel"),
      actionButton("genomeFieldsApply", "Apply", class = "btn-primary")
    ),
    easyClose = TRUE,
    size = "m"
  ))
})

# select all fields
observeEvent(input$genomeFieldsSelectAll, {
  assembly <- state$assembly
  if (is.null(assembly)) return()
  
  available_fields <- get_genome_fields(assembly)
  for (field in available_fields) {
    field_id <- paste0("genome_field_", gsub("[^A-Za-z0-9]", "_", field))
    updateCheckboxInput(session, field_id, value = TRUE)
  }
})

# deselect all fields
observeEvent(input$genomeFieldsDeselectAll, {
  assembly <- state$assembly
  if (is.null(assembly)) return()
  
  available_fields <- get_genome_fields(assembly)
  for (field in available_fields) {
    field_id <- paste0("genome_field_", gsub("[^A-Za-z0-9]", "_", field))
    updateCheckboxInput(session, field_id, value = FALSE)
  }
})

# apply field selection
observeEvent(input$genomeFieldsApply, {
  assembly <- state$assembly
  if (is.null(assembly)) {
    removeModal()
    return()
  }
  
  available_fields <- get_genome_fields(assembly)
  selected_fields <- character()
  
  # collect selected fields
  for (field in available_fields) {
    field_id <- paste0("genome_field_", gsub("[^A-Za-z0-9]", "_", field))
    if (!is.null(input[[field_id]]) && input[[field_id]]) {
      selected_fields <- c(selected_fields, field)
    }
  }
  
  # ensure gid is always included
  if (!"gid" %in% selected_fields) {
    selected_fields <- c("gid", selected_fields)
  }
  
  # save to cache
  set_selected_genome_fields(selected_fields)
  
  # trigger table update
  genome_fields_trigger(genome_fields_trigger() + 1)
  
  removeModal()
  showNotification("Field selection updated", type = "message")
})

