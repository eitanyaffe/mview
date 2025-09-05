# core/gene_regions.r
# Module for creating region files from genes in the current view

# ---- UI Definition ----

gene_regions_ui <- function(id = "gene_regions_module") {
  ns <- shiny::NS(id)
  shiny::tagList()
}

# ---- Server Logic ----

gene_regions_server <- function(id = "gene_regions_module", main_state_rv, session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
      # reactive value to store current genes for region creation
      current_genes <- shiny::reactiveVal(NULL)
      
      # show create regions dialog
      show_create_regions_dialog <- function(genes_df) {
        if (is.null(genes_df) || nrow(genes_df) == 0) {
          shiny::showNotification("no genes in current view", type = "warning")
          return()
        }
        
        # store genes for later use
        current_genes(genes_df)
        
        # check for large number of genes
        gene_count <- nrow(genes_df)
        warning_message <- if (gene_count > 100) {
          shiny::div(
            shiny::p(paste("Create region file with", gene_count, "visible genes (after filtering)")),
            shiny::div(
              style = "color: orange; font-weight: bold; margin: 10px 0;",
              shiny::icon("exclamation-triangle"),
              paste(" Warning: Creating", gene_count, "regions - this may result in a large file!")
            )
          )
        } else {
          shiny::p(paste("Create region file with", gene_count, "visible genes (after filtering)"))
        }
        
        shiny::showModal(shiny::modalDialog(
          title = "Create Regions from Genes",
          shiny::div(
            warning_message,
            shiny::textInput(ns("regions_filename"), "Filename:", 
                           placeholder = "Enter filename (without .txt extension)"),
            shiny::numericInput(ns("window_size_kb"), "Window size (kb):", 
                              value = 20, min = 1, max = 1000, step = 0.1),
            shiny::HTML("<p><small>Window will be centered on gene center</small></p>")
          ),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_create_regions"), "Create")
          ),
          easyClose = TRUE
        ))
      }
      
      # create regions file from genes
      create_regions_from_genes <- function(genes_df, filename, window_size_kb) {
        if (is.null(genes_df) || nrow(genes_df) == 0) {
          return(FALSE)
        }
        
        # ensure regions directory exists
        regions_dir <- get_regions_dir()
        if (!dir.exists(regions_dir)) {
          dir.create(regions_dir, recursive = TRUE)
        }
        
        # check if file exists and confirm overwrite
        file_path <- file.path(regions_dir, paste0(filename, ".txt"))
        if (file.exists(file_path)) {
          return(list(confirm_overwrite = TRUE, file_path = file_path))
        }
        
        # get current state
        current_assembly <- main_state_rv$assembly
        if (is.null(current_assembly)) {
          shiny::showNotification("no assembly selected", type = "error")
          return(FALSE)
        }
        
        # convert window size from kb to bp
        window_size_bp <- window_size_kb * 1000
        half_window <- window_size_bp / 2
        
        
        
        # vectorized calculations (no boundary trimming)
        gene_centers <- (genes_df$start + genes_df$end) / 2
        region_starts <- gene_centers - half_window
        region_ends <- gene_centers + half_window
        
        # create descriptions vectorized
        descriptions <- ifelse(!is.na(genes_df$prot_desc) & 
                              genes_df$prot_desc != "",
                              as.character(genes_df$prot_desc),
                              paste("gene", genes_df$gene))
        
        # create regions data frame (vectorized)
        regions_df <- data.frame(
          id = as.character(genes_df$gene),
          level = 1L,
          description = descriptions,
          assembly = current_assembly,
          contigs = as.character(genes_df$contig),
          zoom_start = as.numeric(region_starts),
          zoom_end = as.numeric(region_ends),
          segment_contig = as.character(genes_df$contig),
          segment_start = as.integer(region_starts),
          segment_end = as.integer(region_ends),
          single_contig = TRUE,
          stringsAsFactors = FALSE
        )
        
        # write to file
        write.table(regions_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        
        shiny::showNotification(paste("created", nrow(regions_df), "regions in", paste0(filename, ".txt")), 
                              type = "message")
        return(TRUE)
      }
      
      # handle create regions confirmation
      observeEvent(input$confirm_create_regions, {
        req(input$regions_filename)
        req(input$window_size_kb)
        
        filename <- trimws(input$regions_filename)
        window_size_kb <- input$window_size_kb
        
        if (filename == "") {
          shiny::showNotification("filename cannot be empty", type = "error")
          return()
        }
        
        if (is.na(window_size_kb) || window_size_kb <= 0) {
          shiny::showNotification("window size must be a positive number", type = "error")
          return()
        }
        
        # use the stored filtered genes
        genes_df <- current_genes()
        
        result <- create_regions_from_genes(genes_df, filename, window_size_kb)
        
        # check if file exists and needs confirmation
        if (is.list(result) && !is.null(result$confirm_overwrite)) {
          shiny::removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "File Already Exists",
            shiny::div(
              shiny::p(paste("The file", paste0(filename, ".txt"), "already exists.")),
              shiny::p("Do you want to overwrite it?")
            ),
            footer = shiny::tagList(
              shiny::modalButton("Cancel"),
              shiny::actionButton(ns("confirm_overwrite"), "Overwrite", class = "btn-warning")
            ),
            easyClose = TRUE
          ))
        } else if (result == TRUE) {
          shiny::removeModal()
        }
      })
      
      # handle overwrite confirmation
      observeEvent(input$confirm_overwrite, {
        req(input$regions_filename)
        req(input$window_size_kb)
        
        filename <- trimws(input$regions_filename)
        window_size_kb <- input$window_size_kb
        
        # use the stored filtered genes and force overwrite
        genes_df <- current_genes()
        
        if (create_regions_from_genes_force(genes_df, filename, window_size_kb)) {
          shiny::removeModal()
        }
      })
      
      # create regions file from genes (force overwrite version)
      create_regions_from_genes_force <- function(genes_df, filename, window_size_kb) {
        if (is.null(genes_df) || nrow(genes_df) == 0) {
          return(FALSE)
        }
        
        # ensure regions directory exists
        regions_dir <- get_regions_dir()
        if (!dir.exists(regions_dir)) {
          dir.create(regions_dir, recursive = TRUE)
        }
        
        file_path <- file.path(regions_dir, paste0(filename, ".txt"))
        
        # get current state
        current_assembly <- main_state_rv$assembly
        if (is.null(current_assembly)) {
          shiny::showNotification("no assembly selected", type = "error")
          return(FALSE)
        }
        
        
        # convert window size from kb to bp
        window_size_bp <- window_size_kb * 1000
        half_window <- window_size_bp / 2
        
        
        
        # vectorized calculations (no boundary trimming)
        gene_centers <- (genes_df$start + genes_df$end) / 2
        region_starts <- gene_centers - half_window
        region_ends <- gene_centers + half_window
        
        # create descriptions vectorized
        descriptions <- ifelse(!is.na(genes_df$prot_desc) & 
                              genes_df$prot_desc != "",
                              as.character(genes_df$prot_desc),
                              paste("gene", genes_df$gene))
        
        # create regions data frame (vectorized)
        regions_df <- data.frame(
          id = as.character(genes_df$gene),
          level = 1L,
          description = descriptions,
          assembly = current_assembly,
          contigs = as.character(genes_df$contig),
          zoom_start = as.numeric(region_starts),
          zoom_end = as.numeric(region_ends),
          segment_contig = as.character(genes_df$contig),
          segment_start = as.integer(region_starts),
          segment_end = as.integer(region_ends),
          single_contig = TRUE,
          stringsAsFactors = FALSE
        )
        
        # write to file (overwriting)
        write.table(regions_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        
        shiny::showNotification(paste("overwritten", nrow(regions_df), "regions in", paste0(filename, ".txt")), 
                              type = "message")
        return(TRUE)
      }
      
      # return function for external use
      return(list(
        show_create_regions_dialog = show_create_regions_dialog
      ))
    }
  )
}
