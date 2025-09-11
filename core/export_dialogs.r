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
  c("Precise", "1kb margin", "10kb margin", "100kb margin", "1mb margin", "All")
}

# convert context option to folder name
context_to_folder_name <- function(context_option) {
  folder_names <- c(
    "Precise" = "precise",
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


# check if variants tab is registered and available
is_variants_tab_available <- function() {
  variants_tab <- get_tab_by_id("variants")
  return(!is.null(variants_tab))
}

# load variant utilities if needed (only when variants are being exported)
load_variant_utilities <- function() {
  # check if functions are already loaded
  if (!exists("query_variants_for_context")) {
    source("tabs/variants/variants_utils.r", local = parent.frame())
  }
  if (!exists("create_variant_frequency_plot")) {
    source("tabs/variants/variants_plots.r", local = parent.frame())
  }
}

# create export directories for all contexts
create_export_dirs <- function(output_dir, context_options) {
  context_dirs <- list()
  
  for (context_option in context_options) {
    context_folder <- context_to_folder_name(context_option)
    context_dir <- file.path(output_dir, context_folder)
    if (!dir.exists(context_dir)) {
      dir.create(context_dir, recursive = TRUE)
    }
    
    # create subfolders for this context
    regions_dir <- file.path(context_dir, "regions")
    profiles_dir <- file.path(context_dir, "profiles") 
    variants_dir <- file.path(context_dir, "variants")
    tables_dir <- file.path(context_dir, "tables")
    
    for (subfolder in c(regions_dir, profiles_dir, variants_dir, tables_dir)) {
      if (!dir.exists(subfolder)) {
        dir.create(subfolder, recursive = TRUE)
      }
    }
    
    context_dirs[[context_option]] <- list(
      context_dir = context_dir,
      regions_dir = regions_dir,
      profiles_dir = profiles_dir,
      variants_dir = variants_dir,
      tables_dir = tables_dir
    )
  }
  
  return(context_dirs)
}

# export a single region in a specific context
export_region <- function(region_data, context_option, dirs, export_params, use_simple_names = FALSE) {
  # set state for this region
  if (!is.null(region_data$contigs) && region_data$contigs != "") {
    state$contigs <- trimws(strsplit(region_data$contigs, ",")[[1]])
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
  
  # build context and plot profiles
  cxt <- build_context(
    state_contigs = state$contigs,
    contig_table = get_contigs(state$assembly),
    zoom = context_zoom,
    assembly = state$assembly
  )

  # early return if no context
  if (is.null(cxt)) return(FALSE)
  
  # generate region name and filenames
  region_name <- if (!is.null(region_data$id)) region_data$id else "region"
  region_name <- sanitize_filename(region_name)
  
  filenames <- if (use_simple_names) {
    generate_single_region_filenames(export_params$include_variants)
  } else {
    generate_multi_region_filenames(region_name, export_params$include_variants)
  }
  
  # query and cache variants FIRST (if requested)
  variant_data <- NULL
  if (export_params$include_variants && !is.null(export_params$variants_tab)) {
    load_variant_utilities()
    variant_data <- query_variants_for_context(
      assembly = state$assembly,
      contigs = state$contigs,
      zoom = context_zoom,
      tab_config = export_params$variants_tab
    )
    
    # cache variants for variant profile to access
    colored_variants <- if (!is.null(variant_data)) add_variant_colors(variant_data$variants) else NULL
    cache_set("variants.current", colored_variants)
  } else {
    cache_set("variants.current", NULL)
  }
  
  # export profiles
  plot_result <- plot_profiles_ggplot(cxt)
  if (is.null(plot_result) || is.null(plot_result$plot)) return(FALSE)
  
  profiles_file <- file.path(dirs$regions_dir, filenames$profiles)
  ggplot2::ggsave(
    filename = profiles_file,
    plot = plot_result$plot,
    width = export_params$width,
    height = export_params$height,
    units = "in",
    dpi = 300,
    device = "pdf"
  )
  
  # copy to profiles folder with simple name
  if (!use_simple_names) {
    profiles_simple_file <- file.path(dirs$profiles_dir, paste0(region_name, ".pdf"))
    file.copy(profiles_file, profiles_simple_file)
  }
  
  # export variants files if requested
  if (!export_params$include_variants || is.null(export_params$variants_tab)) return(TRUE)
  
  if (!is.null(variant_data)) {
    # save variants table
    colored_variants <- add_variant_colors(variant_data$variants)
    table_file <- save_variants_table(colored_variants, dirs$regions_dir, region_name, use_simple_names)
    
    if (!use_simple_names) {
      table_simple_file <- file.path(dirs$tables_dir, paste0(region_name, ".txt"))
      file.copy(table_file, table_simple_file)
    }
    
    # create and save variant plot
    plot_result <- create_variant_frequency_plot(
      variant_data = variant_data,
      library_ids = export_params$variants_tab$library_ids,
      plot_mode = export_params$plot_mode,
      variant_span_filter = export_params$variant_min_span
    )
    
    variants_file <- file.path(dirs$regions_dir, filenames$variants)
    ggplot2::ggsave(
      filename = variants_file,
      plot = plot_result$plot,
      width = plot_result$width_inches,
      height = plot_result$height_inches,
      units = "in",
      dpi = 300,
      device = "pdf"
    )
    
    if (!use_simple_names) {
      variants_simple_file <- file.path(dirs$variants_dir, paste0(region_name, ".pdf"))
      file.copy(variants_file, variants_simple_file)
    }
  } else {
    # create empty variants files
    empty_plot <- create_variant_frequency_plot(NULL, export_params$variants_tab$library_ids)
    
    variants_file <- file.path(dirs$regions_dir, filenames$variants)
    ggplot2::ggsave(
      filename = variants_file,
      plot = empty_plot$plot,
      width = empty_plot$width_inches,
      height = empty_plot$height_inches,
      units = "in",
      dpi = 300,
      device = "pdf"
    )
    
    table_file <- save_variants_table(NULL, dirs$regions_dir, region_name, use_simple_names)
    
    if (!use_simple_names) {
      variants_simple_file <- file.path(dirs$variants_dir, paste0(region_name, ".pdf"))
      file.copy(variants_file, variants_simple_file)
      
      table_simple_file <- file.path(dirs$tables_dir, paste0(region_name, ".txt"))
      file.copy(table_file, table_simple_file)
    }
  }
  
  return(TRUE)  # success
}

# show enhanced export current view dialog
observeEvent(input$plotViewBtn, {
  variants_available <- is_variants_tab_available()
  
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
        p("Creates folders for each context: precise/, 10kb_margin/, etc.", style = "color: #666; font-size: 0.9em;")
      ),
      column(4,
        h5("Options"),
        if (variants_available) {
          tagList(
            checkboxInput("pdf_include_variants", "Include variant plots", 
                         value = cache_get_if_exists("export_include_variants", FALSE), width = "100%"),
            br(),
            numericInput("pdf_variant_min_span", "Min variant span:", 
                        value = cache_get_if_exists("export_variant_min_span", 0.5), min = 0, max = 1, step = 0.1, width = "100%")
          )
        } else {
          p("Variant plots not available", style = "color: #888;")
        }
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
  variants_available <- is_variants_tab_available()
  
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
        p("Creates context folders, each with subfolders: regions/, profiles/, variants/, tables/", style = "color: #666; font-size: 0.9em;")
      ),
      column(4,
        h5("Options"),
        if (variants_available) {
          tagList(
            checkboxInput("pdf_regions_include_variants", "Include variant plots", 
                         value = cache_get_if_exists("export_include_variants", FALSE), width = "100%"),
            br(),
            numericInput("pdf_regions_variant_min_span", "Min variant span:", 
                        value = cache_get_if_exists("export_variant_min_span", 0.5), min = 0, max = 1, step = 0.1, width = "100%")
          )
        } else {
          p("Variant plots not available", style = "color: #888;")
        }
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
    include_variants <- input$pdf_include_variants %||% FALSE
    variant_min_span <- input$pdf_variant_min_span %||% 0.5
    width <- input$pdf_width %||% 10
    height <- input$pdf_height %||% 6
    
    # cache user choices for next time
    cache_set("export_contexts", context_options)
    cache_set("export_output_dir", output_dir_name)
    cache_set("export_include_variants", include_variants)
    cache_set("export_variant_min_span", variant_min_span)
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
    region_name <- if (!is.null(state$zoom)) {
      paste0(paste(state$contigs, collapse = "_"), "_", 
             round(state$zoom[1]), "_", round(state$zoom[2]))
    } else {
      paste(state$contigs, collapse = "_")
    }
    region_name <- sanitize_filename(region_name)
    
    # store original state
    original_contigs <- state$contigs
    original_zoom <- state$zoom
    original_assembly <- state$assembly
    
    # create context directories
    context_dirs <- create_export_dirs(output_dir, context_options)
    
    # get variants tab if needed
    variants_tab <- if (include_variants && is_variants_tab_available()) {
      get_tab_by_id("variants")
    } else {
      NULL
    }
    
    # get current variant settings
    current_plot_mode <- if (exists("input") && !is.null(input$variantPlotMode)) {
      input$variantPlotMode
    } else {
      "frequency"  # default
    }
    
    # create export parameters
    export_params <- list(
      include_variants = include_variants,
      variants_tab = variants_tab,
      variant_min_span = variant_min_span,
      width = width,
      height = height,
      plot_mode = current_plot_mode
    )
    
    # create region data from current state
    current_region <- list(
      id = region_name,
      contigs = paste(state$contigs, collapse = ","),
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
    include_variants <- input$pdf_regions_include_variants %||% FALSE
    variant_min_span <- input$pdf_regions_variant_min_span %||% 0.5
    width <- input$pdf_regions_width %||% 10
    height <- input$pdf_regions_height %||% 6
    
    # cache user choices for next time
    cache_set("export_contexts", context_options)
    cache_set("export_regions_output_dir", output_dir_name)
    cache_set("export_include_variants", include_variants)
    cache_set("export_variant_min_span", variant_min_span)
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
    original_contigs <- state$contigs
    original_zoom <- state$zoom
    original_assembly <- state$assembly
    
    # show progress
    total_regions <- nrow(regions_data)
    total_exports <- total_regions * length(context_options)
    showNotification(paste("Exporting", total_regions, "regions with", length(context_options), "contexts..."), type = "message", duration = 2)
    
    success_count <- 0
    
    # get variants tab if needed
    variants_tab <- if (include_variants && is_variants_tab_available()) {
      # load variant utilities only when needed
      load_variant_utilities()
      get_tab_by_id("variants")
    } else {
      NULL
    }
    
    # get current variant settings
    current_plot_mode <- if (exists("input") && !is.null(input$variantPlotMode)) {
      input$variantPlotMode
    } else {
      "frequency"  # default
    }
    
    # create context directories
    context_dirs <- create_export_dirs(output_dir, context_options)
    
    # create export parameters
    export_params <- list(
      include_variants = include_variants,
      variants_tab = variants_tab,
      variant_min_span = variant_min_span,
      width = width,
      height = height,
      plot_mode = current_plot_mode
    )
    
    # iterate through regions, then contexts (main loop is regions)
    for (i in seq_len(nrow(regions_data))) {
      region <- regions_data[i, ]
      
      # iterate through contexts for this region
      for (context_option in context_options) {
        dirs <- context_dirs[[context_option]]
        
        if (export_region(region, context_option, dirs, export_params, use_simple_names = FALSE)) {
          success_count <- success_count + 1
        }
      }
    }
    
    # restore original state
    state$contigs <- original_contigs
    state$zoom <- original_zoom
    state$assembly <- original_assembly
    
    showNotification(paste("Exported", success_count, "of", total_exports, "region-context combinations to:", output_dir), 
                     type = "message", duration = 8)
    
  }, error = function(e) {
    showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
  })
})
