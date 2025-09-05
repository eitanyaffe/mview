# core/server_pdf.r
# Server-side logic for PDF export functionality

# PDF export using native ggplot2 + patchwork


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

# show PDF export dialog
observeEvent(input$plotViewBtn, {
  showModal(modalDialog(
    title = "Export Current View to PDF",
    size = "m",
    fluidRow(
      column(6,
        numericInput("pdf_width", "Width (inches):", value = 10, min = 1, max = 20, step = 0.5),
        numericInput("pdf_height", "Height (inches):", value = 6, min = 1, max = 20, step = 0.5)
      ),
      column(6,
        textInput("pdf_filename", "Filename:", value = "current_view.pdf", placeholder = "filename.pdf")
      )
    ),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_export_view", "Export PDF", class = "btn-primary")
    )
  ))
})

observeEvent(input$plotRegionsBtn, {
  showModal(modalDialog(
    title = "Export All Regions to PDF",
    size = "m",
    fluidRow(
      column(6,
        numericInput("pdf_regions_width", "Width (inches):", value = 10, min = 1, max = 20, step = 0.5),
        numericInput("pdf_regions_height", "Height (inches):", value = 6, min = 1, max = 20, step = 0.5)
      ),
      column(6,
        textInput("pdf_regions_output_dir", "Output directory:", value = "all_regions", placeholder = "subdirectory under plots/")
      )
    ),
    p("This will export all regions from the current regions table."),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_export_regions", "Export PDFs", class = "btn-primary")
    )
  ))
})

# handle current view export
observeEvent(input$confirm_export_view, {
  removeModal()
  
  # PDF export now uses native ggplot2::ggsave() - no Python dependencies needed
  
  tryCatch({
    # get current plot
    cxt <- build_context(
      state_contigs = state$contigs,
      contig_table = get_contigs(state$assembly),
      zoom = state$zoom,
      assembly = state$assembly
    )
    
    if (is.null(cxt)) {
      showNotification("No context available for plotting", type = "error")
      return()
    }
    
    plot_result <- plot_profiles_ggplot(cxt)
    if (is.null(plot_result) || is.null(plot_result$plot)) {
      showNotification("No plot available for export", type = "error")
      return()
    }
    
    # create plots directory if needed
    plots_dir <- "plots"
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir, recursive = TRUE)
    }
    
    # create output file path
    filename <- input$pdf_filename
    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }
    output_file <- file.path(plots_dir, sanitize_filename(filename))
    
    # export to PDF using ggsave
    tryCatch({
      ggplot2::ggsave(
        filename = output_file,
        plot = plot_result$plot,
        width = input$pdf_width,
        height = input$pdf_height,
        units = "in",
        dpi = 300,
        device = "pdf"
      )
      showNotification(paste("PDF exported to:", output_file), type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Export failed:", e$message), type = "error")
    })
    
  }, error = function(e) {
    showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
  })
})

# handle all regions export
observeEvent(input$confirm_export_regions, {
  removeModal()
  
  # PDF export now uses native ggplot2::ggsave() - no Python dependencies needed
  
  tryCatch({
    # get regions table from the regions module
    regions_data <- regions_module_output$get_region_table()
    if (is.null(regions_data) || nrow(regions_data) == 0) {
      showNotification("No regions available for export", type = "warning")
      return()
    }
    
    # prepare dimensions for regions export
    width_px <- input$pdf_regions_width * 96
    height_px <- input$pdf_regions_height * 96
    
    # create output directory
    output_dir <- file.path("plots", sanitize_filename(input$pdf_regions_output_dir))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # store original state
    original_contigs <- state$contigs
    original_zoom <- state$zoom
    original_assembly <- state$assembly
    
    # show progress
    total_regions <- nrow(regions_data)
    showNotification(paste("Exporting", total_regions, "regions..."), type = "message", duration = 2)
    
    success_count <- 0
    
    # iterate through regions
    for (i in seq_len(nrow(regions_data))) {
      region <- regions_data[i, ]
      
      tryCatch({
        # set state for this region
        if (!is.null(region$contigs) && region$contigs != "") {
          state$contigs <- trimws(strsplit(region$contigs, ",")[[1]])
        }
        if (!is.null(region$zoom_start) && !is.null(region$zoom_end) && 
            !is.na(region$zoom_start) && !is.na(region$zoom_end)) {
          state$zoom <- c(region$zoom_start, region$zoom_end)
        }
        if (!is.null(region$assembly) && region$assembly != "") {
          state$assembly <- region$assembly
        }
        
        # build context and plot
        cxt <- build_context(
          state_contigs = state$contigs,
          contig_table = get_contigs(state$assembly),
          zoom = state$zoom,
          assembly = state$assembly
        )
        
        if (!is.null(cxt)) {
          plot_result <- plot_profiles_ggplot(cxt)
          if (!is.null(plot_result) && !is.null(plot_result$plot)) {
            # create filename from region id
            region_id <- if (!is.null(region$id)) region$id else paste0("region_", i)
            filename <- paste0(sanitize_filename(region_id), ".pdf")
            output_file <- file.path(output_dir, filename)
            
            # export to PDF using ggsave
            tryCatch({
              ggplot2::ggsave(
                filename = output_file,
                plot = plot_result$plot,
                width = input$pdf_regions_width,
                height = input$pdf_regions_height,
                units = "in",
                dpi = 300,
                device = "pdf"
              )
              success_count <- success_count + 1
            }, error = function(e) {
              warning(paste("Failed to save region", region_id, ":", e$message))
            })
          }
        }
        
      }, error = function(e) {
        warning(paste("Failed to export region", i, ":", e$message))
      })
    }
    
    # restore original state
    state$contigs <- original_contigs
    state$zoom <- original_zoom
    state$assembly <- original_assembly
    
    showNotification(paste("Exported", success_count, "of", total_regions, "regions to:", output_dir), 
                     type = "message", duration = 8)
    
  }, error = function(e) {
    showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
  })
})
