# core/export_dialogs.r
# enhanced PDF export functionality with context options and variant plots

library(ggplot2)
library(shinyjs)

# helper function to sanitize filenames
sanitize_filename <- function(name) {
  # remove or replace problematic characters
  name <- gsub("[^a-zA-Z0-9._-]", "_", name)
  # remove leading/trailing underscores
  name <- gsub("^_+|_+$", "", name)
  # ensure not empty
  if (nchar(name) == 0) name <- "unnamed"
  return(name)
}

# get context window options
get_context_options <- function() {
  c("Precise", "-10kb margin", "1kb margin", "10kb margin", "100kb margin", "1mb margin", "All")
}

# convert context option to folder name
context_to_folder_name <- function(context_option) {
  folder_names <- c(
    "Precise" = "precise",
    "-10kb margin" = "minus_10kb_margin",
    "1kb margin" = "1kb_margin", 
    "10kb margin" = "10kb_margin",
    "100kb margin" = "100kb_margin",
    "1mb margin" = "1mb_margin",
    "All" = "all"
  )
  return(folder_names[context_option])
}

# generate filenames for single region export (simple names)
generate_single_region_filenames <- function(include_variants) {
  if (include_variants) {
    return(list(profiles = "profiles.pdf", variants = "variants.pdf"))
  } else {
    return(list(profiles = "profiles.pdf"))
  }
}

# generate filenames for multi-region export (with region prefix)
generate_multi_region_filenames <- function(region_name, include_variants) {
  profiles_filename <- paste0(region_name, "_profiles.pdf")
  
  if (include_variants) {
    variants_filename <- paste0(region_name, "_variants.pdf")
    return(list(profiles = profiles_filename, variants = variants_filename))
  } else {
    return(list(profiles = profiles_filename))
  }
}

# calculate context zoom based on selected option
calculate_context_zoom <- function(current_zoom, context_option, contigs_table) {
  if (is.null(current_zoom)) {
    # if no current zoom, return full range for any option
    return(NULL)
  }
  
  switch(context_option,
    "Precise" = current_zoom,
    "-10kb margin" = {
      # negative margin makes region smaller
      margin <- 10000
      new_start <- current_zoom[1] + margin
      new_end <- current_zoom[2] - margin
      
      # skip regions that would become zero or negative
      if (new_end <= new_start) {
        return(NULL)
      }
      
      c(new_start, new_end)
    },
    "1kb margin" = {
      margin <- 1000
      c(current_zoom[1] - margin, current_zoom[2] + margin)
    },
    "10kb margin" = {
      margin <- 10000
      c(current_zoom[1] - margin, current_zoom[2] + margin)
    },
    "100kb margin" = {
      margin <- 100000
      c(current_zoom[1] - margin, current_zoom[2] + margin)
    },
    "1mb margin" = {
      margin <- 1000000
      c(current_zoom[1] - margin, current_zoom[2] + margin)
    },
    "All" = NULL  # NULL zoom shows all contigs
  )
}


# load required utilities for export
load_export_utilities <- function() {
  # load frequency plots (for PDF export)
  if (!exists("create_frequency_plot_for_export")) {
    source("core/frequency_plots.r", local = parent.frame())
  }
}

# create export directories for all contexts
create_export_dirs <- function(output_dir, context_options, single_region_mode = FALSE) {
  context_dirs <- list()
  exportable_tabs <- get_exportable_tabs()
  
  for (context_option in context_options) {
    context_folder <- context_to_folder_name(context_option)
    context_dir <- file.path(output_dir, context_folder)
    if (!dir.exists(context_dir)) {
      dir.create(context_dir, recursive = TRUE)
    }
    
    if (single_region_mode) {
      # for single region export, use flat structure - all files go directly in context_dir
      dir_list <- list(
        context_dir = context_dir,
        regions_dir = context_dir,  # same as context_dir for flat structure
        profiles_dir = context_dir
      )
      
      # add tab directories (all point to context_dir for flat structure)
      for (tab_id in names(exportable_tabs)) {
        dir_list[[paste0(tab_id, "_dir")]] <- context_dir
      }
    } else {
      # for multi-region export, use organized folder structure
      # create standard subfolders
      regions_dir <- file.path(context_dir, "regions")
      profiles_dir <- file.path(context_dir, "profiles")
      
      subfolders <- c(regions_dir, profiles_dir)
      
      # create dynamic subfolders for each exportable tab
      tab_dirs <- list()
      for (tab_id in names(exportable_tabs)) {
        tab_dir <- file.path(context_dir, tab_id)
        tab_dirs[[tab_id]] <- tab_dir
        subfolders <- c(subfolders, tab_dir)
      }
      
      # create all directories
      for (subfolder in subfolders) {
        if (!dir.exists(subfolder)) {
          dir.create(subfolder, recursive = TRUE)
        }
      }
      
      # build directory list (include tab dirs dynamically)
      dir_list <- list(
        context_dir = context_dir,
        regions_dir = regions_dir,
        profiles_dir = profiles_dir
      )
      
      # add tab directories dynamically
      for (tab_id in names(tab_dirs)) {
        dir_list[[paste0(tab_id, "_dir")]] <- tab_dirs[[tab_id]]
      }
    }
    
    context_dirs[[context_option]] <- dir_list
  }
  
  return(context_dirs)
}

# export a single region in a specific context using new tab-based system
export_region <- function(region_data, context_option, dirs, export_params, use_simple_names = FALSE) {
  # set state for this region
  if (!is.null(region_data$contigs) && region_data$contigs != "") {
    contig_list <- trimws(strsplit(region_data$contigs, ",")[[1]])
    # get segments for these contigs
    segments <- get_segments(state$assembly)
    selected_segments <- segments[segments$contig %in% contig_list, ]
    state$segments <- selected_segments
  }

  region_zoom <- NULL
  if (!is.null(region_data$zoom_start) && !is.null(region_data$zoom_end) && 
      !is.na(region_data$zoom_start) && !is.na(region_data$zoom_end)) {
    region_zoom <- c(region_data$zoom_start, region_data$zoom_end)
  }
  state$zoom <- region_zoom

  if (!is.null(region_data$assembly) && region_data$assembly != "") {
    state$assembly <- region_data$assembly
  }

  # calculate context zoom for this region and context
  context_zoom <- calculate_context_zoom(region_zoom, context_option, get_contigs(state$assembly))
  
  # skip regions that would have zero or negative size (e.g., with negative margins)
  if (context_option == "-10kb margin" && is.null(context_zoom)) {
    return(FALSE)
  }
  
  # save current zoom and set export zoom temporarily
  saved_zoom <- state$zoom
  if (!is.null(context_zoom)) {
    cxt_set_zoom(context_zoom)
  }
  
  # generate region name
  region_name <- if (!is.null(region_data$id)) region_data$id else "region"
  region_name <- sanitize_filename(region_name)

  # create region_info for tabs
  region_info <- list(
    assembly = state$assembly,
    contigs = get_state_contigs(),
    zoom = state$zoom,
    context_zoom = context_zoom,
    region_name = region_name
  )

  # export tab plots and tables FIRST (to populate cache for profile rendering)
  exportable_tabs <- get_exportable_tabs()
  for (tab_id in names(exportable_tabs)) {
    # check if this tab was requested for export
    include_param <- paste0("include_", tab_id)
    if (!isTRUE(export_params[[include_param]])) {
      next  # skip this tab
    }
    
    # export PDF if available
    pdf_export_function <- get_tab_export_function(tab_id, "pdf")
    if (!is.null(pdf_export_function)) {
      # call the tab's PDF export function
      plot_result <- pdf_export_function(region_info)
      if (!is.null(plot_result) && !is.null(plot_result$plot)) {
        # save the plot
        if (use_simple_names) {
          # for single region export, use simple filename directly in tab folder
          tab_dir_key <- paste0(tab_id, "_dir")
          tab_file <- file.path(dirs[[tab_dir_key]], paste0(tab_id, ".pdf"))
        } else {
          # for multi-region export, use region prefix in regions directory
          tab_file <- file.path(dirs$regions_dir, paste0(region_name, "_", tab_id, ".pdf"))
        }
        ggplot2::ggsave(
          filename = tab_file,
          plot = plot_result$plot,
          width = plot_result$width_inches,
          height = plot_result$height_inches,
          units = "in",
          dpi = 300,
          device = "pdf"
        )
        
        # copy to tab-specific folder with simple name (only for multi-region)
        if (!use_simple_names) {
          tab_dir_key <- paste0(tab_id, "_dir")
          if (!is.null(dirs[[tab_dir_key]])) {
            tab_simple_file <- file.path(dirs[[tab_dir_key]], paste0(region_name, ".pdf"))
            file.copy(tab_file, tab_simple_file)
          }
        }
      } else {
        cat(sprintf("warning: PDF export function for tab '%s' returned no plot, skipping PDF\n", tab_id))
      }
    }
    
    # export table if available
    table_export_function <- get_tab_export_function(tab_id, "table")
    if (!is.null(table_export_function)) {
      # call the tab's table export function
      table_data <- table_export_function(region_info)
      
      if (!is.null(table_data)) {
        # save the table
        if (use_simple_names) {
          # for single region export, use simple filename directly in tab folder
          tab_dir_key <- paste0(tab_id, "_dir")
          table_file <- file.path(dirs[[tab_dir_key]], paste0(tab_id, "_table.txt"))
        } else {
          # for multi-region export, use region prefix in regions directory
          table_file <- file.path(dirs$regions_dir, paste0(region_name, "_", tab_id, "_table.txt"))
        }
        
        write.table(table_data, table_file, sep = "\t", row.names = FALSE, quote = FALSE)
        cat(sprintf("saved %s table (%d rows) to: %s\n", tab_id, nrow(table_data), table_file))
        
        # copy to tab-specific folder with simple name (only for multi-region)
        if (!use_simple_names) {
          tab_dir_key <- paste0(tab_id, "_dir")
          if (!is.null(dirs[[tab_dir_key]])) {
            tab_simple_table_file <- file.path(dirs[[tab_dir_key]], paste0(region_name, "_table.txt"))
            file.copy(table_file, tab_simple_table_file)
          }
        }
      } else {
        cat(sprintf("warning: table export function for tab '%s' returned no data, skipping table\n", tab_id))
      }
    }
    
    # warn if no export functions are available
    if (is.null(pdf_export_function) && is.null(table_export_function)) {
      cat(sprintf("warning: no export functions registered for tab '%s', skipping\n", tab_id))
    }
  }
  
  # NOW plot profiles (after cache has been populated by tab export functions)
  plot_result <- plot_profiles_ggplot()
  
  # restore original zoom
  if (!is.null(saved_zoom)) {
    cxt_set_zoom(saved_zoom)
  }
  
  if (is.null(plot_result) || is.null(plot_result$plot)) return(FALSE)

  # save profiles
  if (use_simple_names) {
    # for single region export, use simple filename directly
    profiles_file <- file.path(dirs$profiles_dir, "profiles.pdf")
  } else {
    # for multi-region export, use region prefix
    profiles_file <- file.path(dirs$regions_dir, paste0(region_name, "_profiles.pdf"))
  }

  ggplot2::ggsave(
    filename = profiles_file,
    plot = plot_result$plot,
    width = export_params$width,
    height = export_params$height,
    units = "in",
    dpi = 300,
    device = "pdf"
  )
  
  # copy to profiles folder with simple name (only for multi-region)
  if (!use_simple_names) {
    profiles_simple_file <- file.path(dirs$profiles_dir, paste0(region_name, ".pdf"))
    file.copy(profiles_file, profiles_simple_file)
  }
  
  return(TRUE)  # success
}

# show enhanced export current view dialog
observeEvent(input$plotViewBtn, {
  exportable_tabs <- get_exportable_tabs()
  
  # create tab checkboxes dynamically
  tab_checkboxes <- if (length(exportable_tabs) > 0) {
    tab_inputs <- lapply(names(exportable_tabs), function(tab_id) {
      tab <- exportable_tabs[[tab_id]]
      checkboxInput(paste0("pdf_include_", tab_id), 
                   paste("Include", tolower(tab$tab_label), "plots"),
                   value = cache_get_if_exists(paste0("export_include_", tab_id), FALSE), 
                   width = "100%")
    })
    do.call(tagList, tab_inputs)
  } else {
    p("No exportable tabs available", style = "color: #888;")
  }
  
  showModal(modalDialog(
    title = "Export Current View to PDF",
    size = "l",
    fluidRow(
      column(4,
        h5("Context"),
        checkboxGroupInput("pdf_contexts", "Windows (select multiple):",
                          choices = get_context_options(),
                          selected = cache_get_if_exists("export_contexts", "10kb margin"), width = "100%"),
        br(),
        h5("Dimensions"),
        numericInput("pdf_width", "Width (inches):", 
                    value = cache_get_if_exists("export_width", 10), min = 1, max = 20, step = 0.5),
        numericInput("pdf_height", "Height (inches):", 
                    value = cache_get_if_exists("export_height", 6), min = 1, max = 20, step = 0.5)
      ),
      column(4,
        h5("Output"),
        fluidRow(
          column(8,
            textInput("pdf_output_dir", "Folder name:", 
                     value = cache_get_if_exists("export_output_dir", "current_view"), placeholder = "folder name")
          ),
          column(4,
            br(),
            actionButton("pdf_browse_dir", "Browse...", class = "btn-secondary", width = "100%")
          )
        ),
        br(),
        p(id = "pdf_full_path", paste("Full path:", file.path(getwd(), "plots", cache_get_if_exists("export_output_dir", "current_view"))), 
          style = "color: #666; font-size: 0.9em; word-break: break-all;"),
        br(),
        p("Creates folders for each context and tab type.", style = "color: #666; font-size: 0.9em;")
      ),
      column(4,
        h5("Tabs to Export"),
        tab_checkboxes
      )
    ),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_export_view", "Export PDF", class = "btn-primary")
    )
  ))
})

# show enhanced export all regions dialog  
observeEvent(input$plotRegionsBtn, {
  exportable_tabs <- get_exportable_tabs()
  
  # create tab checkboxes dynamically
  tab_checkboxes <- if (length(exportable_tabs) > 0) {
    tab_inputs <- lapply(names(exportable_tabs), function(tab_id) {
      tab <- exportable_tabs[[tab_id]]
      checkboxInput(paste0("pdf_regions_include_", tab_id), 
                   paste("Include", tolower(tab$tab_label), "plots"),
                   value = cache_get_if_exists(paste0("export_include_", tab_id), FALSE), 
                   width = "100%")
    })
    do.call(tagList, tab_inputs)
  } else {
    p("No exportable tabs available", style = "color: #888;")
  }
  
  showModal(modalDialog(
    title = "Export All Regions to PDF",
    size = "l",
    fluidRow(
      column(4,
        h5("Context"),
        checkboxGroupInput("pdf_regions_contexts", "Windows (select multiple):",
                          choices = get_context_options(),
                          selected = cache_get_if_exists("export_contexts", "10kb margin"), width = "100%"),
        br(),
        h5("Dimensions"),
        numericInput("pdf_regions_width", "Width (inches):", 
                    value = cache_get_if_exists("export_width", 10), min = 1, max = 20, step = 0.5),
        numericInput("pdf_regions_height", "Height (inches):", 
                    value = cache_get_if_exists("export_height", 6), min = 1, max = 20, step = 0.5)
      ),
      column(4,
        h5("Output"),
        fluidRow(
          column(8,
            textInput("pdf_regions_output_dir", "Folder name:", 
                     value = cache_get_if_exists("export_regions_output_dir", "all_regions"), placeholder = "folder name")
          ),
          column(4,
            br(),
            actionButton("pdf_regions_browse_dir", "Browse...", class = "btn-secondary", width = "100%")
          )
        ),
        br(),
        p(id = "pdf_regions_full_path", paste("Full path:", file.path(getwd(), "plots", cache_get_if_exists("export_regions_output_dir", "all_regions"))), 
          style = "color: #666; font-size: 0.9em; word-break: break-all;"),
        br(),
        p("Creates context folders, each with subfolders for profiles and tab types.", style = "color: #666; font-size: 0.9em;")
      ),
      column(4,
        h5("Tabs to Export"),
        tab_checkboxes
      )
    ),
    p("This will export all regions from the current regions table."),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_export_regions", "Export PDFs", class = "btn-primary")
    )
  ))
})

# update full path display when user types in single export dialog
observeEvent(input$pdf_output_dir, {
  if (!is.null(input$pdf_output_dir)) {
    full_path <- file.path(getwd(), "plots", input$pdf_output_dir)
    runjs(sprintf("$('#pdf_full_path').text('Full path: %s');", full_path))
  }
})

# update full path display when user types in regions export dialog  
observeEvent(input$pdf_regions_output_dir, {
  if (!is.null(input$pdf_regions_output_dir)) {
    full_path <- file.path(getwd(), "plots", input$pdf_regions_output_dir)
    runjs(sprintf("$('#pdf_regions_full_path').text('Full path: %s');", full_path))
  }
})

# browse button for single export
observeEvent(input$pdf_browse_dir, {
  tryCatch({
    # use tcltk to open folder dialog
    if (requireNamespace("tcltk", quietly = TRUE)) {
      selected_dir <- tcltk::tk_choose.dir(default = file.path(getwd(), "plots"))
      if (!is.na(selected_dir) && selected_dir != "") {
        # extract relative path from plots directory or use absolute path
        plots_dir <- file.path(getwd(), "plots")
        if (startsWith(selected_dir, plots_dir)) {
          # relative to plots directory
          rel_path <- sub(paste0(plots_dir, .Platform$file.sep), "", selected_dir)
          rel_path <- sub(paste0("^", plots_dir, "$"), ".", rel_path)  # if exact match to plots dir
        } else {
          # use absolute path
          rel_path <- selected_dir
        }
        updateTextInput(session, "pdf_output_dir", value = rel_path)
      }
    } else {
      showNotification("Folder browser requires tcltk package", type = "warning")
    }
  }, error = function(e) {
    showNotification("Could not open folder browser", type = "warning")
  })
})

# browse button for regions export
observeEvent(input$pdf_regions_browse_dir, {
  tryCatch({
    # use tcltk to open folder dialog
    if (requireNamespace("tcltk", quietly = TRUE)) {
      selected_dir <- tcltk::tk_choose.dir(default = file.path(getwd(), "plots"))
      if (!is.na(selected_dir) && selected_dir != "") {
        # extract relative path from plots directory or use absolute path
        plots_dir <- file.path(getwd(), "plots")
        if (startsWith(selected_dir, plots_dir)) {
          # relative to plots directory
          rel_path <- sub(paste0(plots_dir, .Platform$file.sep), "", selected_dir)
          rel_path <- sub(paste0("^", plots_dir, "$"), ".", rel_path)  # if exact match to plots dir
        } else {
          # use absolute path
          rel_path <- selected_dir
        }
        updateTextInput(session, "pdf_regions_output_dir", value = rel_path)
      }
    } else {
      showNotification("Folder browser requires tcltk package", type = "warning")
    }
  }, error = function(e) {
    showNotification("Could not open folder browser", type = "warning")
  })
})

# export single view with enhanced options
observeEvent(input$confirm_export_view, {
  removeModal()
  
  tryCatch({
    # get export parameters
    context_options <- input$pdf_contexts %||% c("10kb margin")
    output_dir_name <- input$pdf_output_dir %||% "current_view"
    width <- input$pdf_width %||% 10
    height <- input$pdf_height %||% 6
    
    # get tab inclusion settings dynamically
    exportable_tabs <- get_exportable_tabs()
    tab_inclusions <- list()
    for (tab_id in names(exportable_tabs)) {
      include_param <- paste0("include_", tab_id)
      input_param <- paste0("pdf_include_", tab_id)
      tab_inclusions[[include_param]] <- input[[input_param]] %||% FALSE
      # cache the setting
      cache_set(paste0("export_include_", tab_id), tab_inclusions[[include_param]])
    }
    
    # cache other settings
    cache_set("export_contexts", context_options)
    cache_set("export_output_dir", output_dir_name)
    cache_set("export_width", width)
    cache_set("export_height", height)
    
    # validate that at least one context is selected
    if (length(context_options) == 0) {
      showNotification("Please select at least one context window", type = "warning")
      return()
    }
    
    # create output directory
    output_dir <- file.path("plots", sanitize_filename(output_dir_name))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # generate region name for filenames
    current_contigs <- get_state_contigs()
    region_name <- if (!is.null(state$zoom)) {
      paste0(paste(current_contigs, collapse = "_"), "_", 
             round(state$zoom[1]), "_", round(state$zoom[2]))
    } else {
      paste(current_contigs, collapse = "_")
    }
    region_name <- sanitize_filename(region_name)
    
    # store original state
    original_segments <- get_state_segments()
    original_zoom <- state$zoom
    original_assembly <- state$assembly
    
    # create context directories (use flat structure for single region)
    context_dirs <- create_export_dirs(output_dir, context_options, single_region_mode = TRUE)
    
    # create export parameters (include tab inclusions and basic settings)
    export_params <- c(
      list(
        width = width,
        height = height
      ),
      tab_inclusions
    )
    
    # create region data from current state
    current_region <- list(
      id = region_name,
      contigs = paste(get_state_contigs(), collapse = ","),
      zoom_start = if (!is.null(state$zoom)) state$zoom[1] else NULL,
      zoom_end = if (!is.null(state$zoom)) state$zoom[2] else NULL,
      assembly = state$assembly
    )
    
    # export for each selected context
    success_count <- 0
    for (context_option in context_options) {
      dirs <- context_dirs[[context_option]]
      
      if (export_region(current_region, context_option, dirs, export_params, use_simple_names = TRUE)) {
        success_count <- success_count + 1
      }
    }
    
    # restore original state
    state$contigs <- original_contigs
    state$zoom <- original_zoom
    state$assembly <- original_assembly
    
    showNotification(paste("Exported", success_count, "of", length(context_options), "contexts to:", output_dir), type = "message", duration = 5)
    
  }, error = function(e) {
    showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
  })
})

# export all regions with enhanced options
observeEvent(input$confirm_export_regions, {
  removeModal()
  
  tryCatch({
    # get export parameters
    context_options <- input$pdf_regions_contexts %||% c("10kb margin")
    output_dir_name <- input$pdf_regions_output_dir %||% "all_regions"
    width <- input$pdf_regions_width %||% 10
    height <- input$pdf_regions_height %||% 6
    
    # get tab inclusion settings dynamically
    exportable_tabs <- get_exportable_tabs()
    tab_inclusions <- list()
    for (tab_id in names(exportable_tabs)) {
      include_param <- paste0("include_", tab_id)
      input_param <- paste0("pdf_regions_include_", tab_id)
      tab_inclusions[[include_param]] <- input[[input_param]] %||% FALSE
      # cache the setting
      cache_set(paste0("export_include_", tab_id), tab_inclusions[[include_param]])
    }
    
    # cache other settings
    cache_set("export_contexts", context_options)
    cache_set("export_regions_output_dir", output_dir_name)
    cache_set("export_width", width)
    cache_set("export_height", height)
    
    # validate that at least one context is selected
    if (length(context_options) == 0) {
      showNotification("Please select at least one context window", type = "warning")
      return()
    }
    
    # get regions table from the regions module
    regions_data <- regions_module_output$get_region_table()
    if (is.null(regions_data) || nrow(regions_data) == 0) {
      showNotification("No regions available for export", type = "warning")
      return()
    }
    
    # create output directory and subfolders for multi-region export
    output_dir <- file.path("plots", sanitize_filename(output_dir_name))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # store original state
    original_segments <- get_state_segments()
    original_zoom <- state$zoom
    original_assembly <- state$assembly
    
    # show progress
    total_regions <- nrow(regions_data)
    total_exports <- total_regions * length(context_options)
    showNotification(paste("Exporting", total_regions, "regions with", length(context_options), "contexts..."), type = "message", duration = 2)
    
    success_count <- 0
    
    # create context directories
    context_dirs <- create_export_dirs(output_dir, context_options)
    
    # create export parameters (include tab inclusions and basic settings)
    export_params <- c(
      list(
        width = width,
        height = height
      ),
      tab_inclusions
    )
    
    # iterate through regions, then contexts (main loop is regions)
    for (i in seq_len(nrow(regions_data))) {
      region <- regions_data[i, ]
      
      cat(sprintf("---------------------------------------- Region %d/%d: %s ----------------------------------------\n", 
                  i, nrow(regions_data), region$id))
      
      # iterate through contexts for this region
      for (context_option in context_options) {
        dirs <- context_dirs[[context_option]]
        
        if (export_region(region, context_option, dirs, export_params, use_simple_names = FALSE)) {
          success_count <- success_count + 1
        }
      }
    }
    
    # restore original state
    state$segments <- original_segments
    state$zoom <- original_zoom
    state$assembly <- original_assembly
    
    showNotification(paste("Exported", success_count, "of", total_exports, "region-context combinations to:", output_dir), 
                     type = "message", duration = 8)
    
  }, error = function(e) {
    showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
  })
})
