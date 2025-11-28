# core/regions.r
# Defines the UI and server logic for the region management module.

# Hardcoded fields from the main application region to be tracked
TRACKED_REGION_FIELDS <- c("contigs", "zoom", "assembly")

# ---- UI Definition ----

regions_ui <- function(id = "regions_module") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::div(
          shiny::tags$label("Regions File:", `for` = ns("regions_file_input")),
          shiny::tags$input(id = ns("regions_file_input"), type = "text", class = "form-control", 
                           style = "width: 200px;", readonly = "readonly", value = "")
        ),
        shiny::actionButton(ns("new_region_table"), "New Table"),
        shiny::actionButton(ns("open_region_table"), "Open Table"),
        shiny::actionButton(ns("delete_region_table"), "Delete Table"),
        shiny::actionButton(ns("homologs_btn"), "Homologs", icon = shiny::icon("search")),
        shiny::hr(),
        shiny::actionButton(ns("goto_region_btn"), "Goto"),
        shiny::actionButton(ns("save_region"), "Save"),
        shiny::actionButton(ns("edit_region"), "Edit"),
        shiny::actionButton(ns("delete_region"), "Delete")
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::h4("Saved Regions"),
        DT::DTOutput(ns("region_table_display"))
      )
    )
  )
}

# ---- Server Logic ----

regions_server <- function(id = "regions_module", main_state_rv, session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # --- Helper Functions ---
      create_empty_region_table <- function() {
        data.frame(
          id = character(),
          level = integer(),
          description = character(),
          assembly = character(),
          contigs = character(),
          zoom_start = numeric(),
          zoom_end = numeric(),
          segment_contig = character(),
          segment_start = integer(),
          segment_end = integer(),
          single_contig = logical(),
          stringsAsFactors = FALSE
        )
      }
      
      # migrate region table to add missing columns
      migrate_region_table <- function(rt) {
        required_new_cols <- c("level", "segment_contig", "segment_start", "segment_end", "single_contig")
        missing_cols <- setdiff(required_new_cols, colnames(rt))
        
        if (length(missing_cols) == 0) {
          return(rt)
        }
        
        cat("migrating regions table: adding", paste(missing_cols, collapse = ", "), "columns\n")
        
        # add missing columns with default values
        for (col in missing_cols) {
          if (col == "level") {
            rt[[col]] <- rep(3L, nrow(rt))  # default level is 1
          } else if (col == "segment_contig") {
            rt[[col]] <- character(nrow(rt))
          } else if (col %in% c("segment_start", "segment_end")) {
            rt[[col]] <- integer(nrow(rt))
          } else if (col == "single_contig") {
            rt[[col]] <- logical(nrow(rt))
          }
        }
        
        # ensure level column is in correct position (second column after id)
        if ("level" %in% missing_cols && "id" %in% colnames(rt)) {
          col_order <- colnames(rt)
          id_pos <- which(col_order == "id")
          level_pos <- which(col_order == "level")
          
          # reorder columns to put level right after id
          new_order <- c(col_order[1:id_pos], "level", col_order[setdiff(seq_along(col_order), c(id_pos, level_pos))])
          rt <- rt[, new_order, drop = FALSE]
        }
        
        # compute segment information for each row
        for (i in seq_len(nrow(rt))) {
          segment_info <- compute_segment_info(rt[i, ])
          rt$segment_contig[i] <- segment_info$contig
          rt$segment_start[i] <- segment_info$start
          rt$segment_end[i] <- segment_info$end
          rt$single_contig[i] <- segment_info$single_contig
        }
        
        return(rt)
      }
      
      # compute segment information for a single region row
      compute_segment_info <- function(region_row) {
        # parse contigs
        contigs_list <- if (region_row$contigs == "" || is.na(region_row$contigs)) {
          character(0)
        } else {
          trimws(strsplit(region_row$contigs, ",")[[1]])
        }
        
        # default values for multi-contig or invalid regions
        default_result <- list(
          contig = NA_character_,
          start = 0L,
          end = 0L,
          single_contig = FALSE
        )
        
        if (length(contigs_list) == 0 || 
            is.na(region_row$zoom_start) || is.na(region_row$zoom_end)) {
          return(default_result)
        }
        
        # get segments for the contigs in this region
        all_segments <- get_segments(region_row$assembly)
        if (is.null(all_segments)) {
          return(default_result)
        }
        region_segments <- all_segments[all_segments$contig %in% contigs_list, ]
        if (nrow(region_segments) == 0) {
          return(default_result)
        }
        
        # save current context state
        saved_assembly <- cxt_get_assembly()
        saved_segments <- tryCatch(get_state_segments(), error = function(e) NULL)
        saved_zoom <- tryCatch(cxt_get_xlim(), error = function(e) NULL)
        
        # set up context for this region temporarily
        tryCatch({
          cxt_set_assembly(region_row$assembly)
          cxt_set_view(region_segments)
          cxt_set_zoom(c(region_row$zoom_start, region_row$zoom_end))
          
          # check which contigs the region spans
          cdf <- cxt_get_entire_view()
          zoom_contigs <- cdf[cdf$vstart < region_row$zoom_end & cdf$vend > region_row$zoom_start, ]
          
          result <- if (nrow(zoom_contigs) != 1) {
            # spans multiple contigs or no contigs
            default_result
          } else {
            # single contig case - convert virtual to local coordinates
            contig_name <- zoom_contigs$contig[1]
            vstart <- zoom_contigs$vstart[1]
            local_start_offset <- zoom_contigs$start[1]
            local_start <- round(region_row$zoom_start - vstart) + local_start_offset
            local_end <- round(region_row$zoom_end - vstart) + local_start_offset
            
            list(
              contig = contig_name,
              start = as.integer(local_start),
              end = as.integer(local_end),
              single_contig = TRUE
            )
          }
          
          # restore original context
          if (!is.null(saved_assembly) && !is.null(saved_segments)) {
            cxt_set_assembly(saved_assembly)
            cxt_set_view(saved_segments)
            if (!is.null(saved_zoom)) {
              cxt_set_zoom(saved_zoom)
            }
          }
          
          return(result)
        }, error = function(e) {
          # restore on error
          if (!is.null(saved_assembly) && !is.null(saved_segments)) {
            cxt_set_assembly(saved_assembly)
            cxt_set_view(saved_segments)
            if (!is.null(saved_zoom)) {
              cxt_set_zoom(saved_zoom)
            }
          }
          return(default_result)
        })
        
        return(list(
          contig = contig_name,
          start = as.integer(local_start),
          end = as.integer(local_end),
          single_contig = TRUE
        ))
      }

      format_contigs_for_display <- function(contigs_vector) {
        if (length(contigs_vector) == 0) {
          return("None")
        }
        if (length(contigs_vector) == 1) {
          return(as.character(contigs_vector[1]))
        }
        return(paste(length(contigs_vector), "contigs"))
      }

      format_zoom_for_display <- function(zoom_vector) {
        if (is.null(zoom_vector) || length(zoom_vector) != 2) {
          return("Full range")
        }
        return(paste(round(zoom_vector[1]), "–", round(zoom_vector[2])))
      }

      # --- Reactive Values ---
      region_table <- shiny::reactiveVal(create_empty_region_table())
      undo_stack <- shiny::reactiveVal(list())
      current_regions_file <- shiny::reactiveVal("")
      
      # function to update current file and cache it
      set_current_regions_file <- function(filename) {
        current_regions_file(filename)
        cache_set("regions.current_file", filename)
        session$sendCustomMessage("updateReadonlyInput", list(id = ns("regions_file_input"), value = filename))
      }
      
      # auto-load cached regions file on startup
      shiny::observe({
        if (cache_exists("regions.current_file") && current_regions_file() == "") {
          cached_file <- cache_get("regions.current_file")
          if (cached_file != "" && file.exists(file.path(get_regions_dir(), cached_file))) {
            load_region_table(cached_file)
            # ensure UI updates after session is ready
            session$onFlushed(function() {
              session$sendCustomMessage("updateReadonlyInput", list(id = ns("regions_file_input"), value = cached_file))
            }, once = TRUE)
          }
        }
      })
      
      # update cache when regions table changes
      shiny::observe({
        rt <- region_table()
        segments_data <- convert_regions_to_segments(rt)
        cache_set("segments.current_regions", segments_data)
      })
      
      # convert regions table to segments format
      convert_regions_to_segments <- function(rt) {
        if (is.null(rt) || nrow(rt) == 0) {
          return(data.frame(
            assembly = character(),
            contig = character(),
            start = integer(),
            end = integer(),
            desc = character(),
            id = character(),
            stringsAsFactors = FALSE
          ))
        }
        
        # filter for single-contig segments only
        single_contig_rows <- rt[rt$single_contig %in% TRUE & !is.na(rt$segment_contig), ]
        
        if (nrow(single_contig_rows) == 0) {
          return(data.frame(
            assembly = character(),
            contig = character(),
            start = integer(),
            end = integer(),
            desc = character(),
            id = character(),
            stringsAsFactors = FALSE
          ))
        }
        
        data.frame(
          assembly = single_contig_rows$assembly,
          contig = single_contig_rows$segment_contig,
          start = single_contig_rows$segment_start,
          end = single_contig_rows$segment_end,
          desc = single_contig_rows$description,
          id = single_contig_rows$id,
          stringsAsFactors = FALSE
        )
      }

      # ensure regions directory exists
      ensure_regions_dir <- function() {
        regions_dir <- get_regions_dir()
        if (!dir.exists(regions_dir)) {
          dir.create(regions_dir, recursive = TRUE)
        }
        # ensure versions subdirectory exists
        versions_dir <- file.path(regions_dir, "versions")
        if (!dir.exists(versions_dir)) {
          dir.create(versions_dir)
        }
        return(regions_dir)
      }

      # save region table with versioning
      save_region_table <- function(filename) {
        regions_dir <- ensure_regions_dir()
        file_path <- file.path(regions_dir, filename)
        
        # create backup version if file exists
        if (file.exists(file_path)) {
          versions_dir <- file.path(regions_dir, "versions")
          base_name <- sub("\\.txt$", "", filename)
          existing_versions <- list.files(versions_dir, pattern = paste0("^", base_name, "_v\\d+\\.txt$"))
          version_num <- length(existing_versions) + 1
          version_filename <- paste0(base_name, "_v", version_num, ".txt")
          file.copy(file_path, file.path(versions_dir, version_filename))
        }
        
        write.table(region_table(), file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        shiny::showNotification(paste("regions table saved to", filename), type = "message")
      }

      # load region table
      load_region_table <- function(filename) {
        regions_dir <- get_regions_dir()
        file_path <- file.path(regions_dir, filename)
        
        if (!file.exists(file_path)) {
          stop(paste("region file not found:", file_path))
        }
        
        df <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        
        # validate required columns
        required_cols <- c("id", "description", "assembly", "contigs", "zoom_start", "zoom_end")
        if (!all(required_cols %in% colnames(df))) {
          stop(paste("region file missing required columns:", paste(setdiff(required_cols, colnames(df)), collapse = ", ")))
        }
        
        # migrate table to add segment columns if needed
        df <- migrate_region_table(df)
        
        region_table(df)
        set_current_regions_file(filename)
      }

      # --- Event Handlers ---

      # new table button
      observeEvent(input$new_region_table, {
        shiny::showModal(shiny::modalDialog(
          title = "New Regions Table",
          shiny::textInput(ns("new_table_name"), "Table Name:", placeholder = "Enter table name"),
          shiny::HTML("<p><small>Name will be saved as [name].txt</small></p>"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_new_table"), "Create")
          ),
          easyClose = TRUE
        ))
      })

      observeEvent(input$confirm_new_table, {
        req(input$new_table_name)
        table_name <- trimws(input$new_table_name)
        
        if (table_name == "") {
          shiny::showNotification("table name cannot be empty", type = "error")
          return()
        }
        
        filename <- paste0(table_name, ".txt")
        empty_table <- create_empty_region_table()
        region_table(empty_table)
        
        # immediately create the physical file with headers
        regions_dir <- ensure_regions_dir()
        file_path <- file.path(regions_dir, filename)
        write.table(empty_table, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        
        set_current_regions_file(filename)
        shiny::removeModal()
        shiny::showNotification(paste("created new regions table:", filename), type = "message")
      })

      # open table button
      observeEvent(input$open_region_table, {
        regions_dir <- tryCatch(get_regions_dir(), error = function(e) {
          shiny::showNotification("regions directory not configured", type = "error")
          return(NULL)
        })
        
        if (is.null(regions_dir)) return()
        
        if (!dir.exists(regions_dir)) {
          shiny::showNotification("regions directory does not exist", type = "error")
          return()
        }
        
        files <- list.files(regions_dir, pattern = "\\.txt$", full.names = FALSE)
        if (length(files) == 0) {
          shiny::showModal(shiny::modalDialog(
            title = "Open Regions Table",
            "No region files found in regions directory.",
            footer = shiny::modalButton("OK")
          ))
          return()
        }

        shiny::showModal(shiny::modalDialog(
          title = "Open Regions Table",
          shiny::selectInput(ns("select_region_file"),
            label = "Choose a file:",
            choices = files
          ),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_open_table"), "Open")
          ),
          easyClose = TRUE
        ))
      })

      observeEvent(input$confirm_open_table, {
        req(input$select_region_file)
        tryCatch({
          load_region_table(input$select_region_file)
          shiny::removeModal()
        }, error = function(e) {
          shiny::showNotification(paste("error loading regions table:", e$message), type = "error")
        })
      })

      # delete table button
      observeEvent(input$delete_region_table, {
        if (current_regions_file() == "") {
          shiny::showNotification("no regions table selected", type = "error")
          return()
        }
        
        shiny::showModal(shiny::modalDialog(
          title = "Delete Regions Table",
          paste("are you sure you want to delete the entire table:", current_regions_file(), "?"),
          shiny::HTML("<p><strong>This action cannot be undone!</strong></p>"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_delete_table"), "Delete", class = "btn-danger")
          )
        ))
      })

      observeEvent(input$confirm_delete_table, {
        regions_dir <- get_regions_dir()
        file_path <- file.path(regions_dir, current_regions_file())
        
        if (file.exists(file_path)) {
          file.remove(file_path)
          region_table(create_empty_region_table())
          set_current_regions_file("")
          shiny::showNotification("regions table deleted", type = "message")
        } else {
          shiny::showNotification("file not found", type = "error")
        }
        
        shiny::removeModal()
      })

      # homologs search button
      observeEvent(input$homologs_btn, {
        # get current view information from main state
        current_app_state <- shiny::isolate(shiny::reactiveValuesToList(main_state_rv))
        
        # determine query contig and convert global coordinates to local coordinates
        query_contig <- ""
        query_start <- 1L
        query_end <- 1000L
        
        current_segments <- current_app_state$segments
        if (!is.null(current_segments) && nrow(current_segments) > 0 &&
            !is.null(current_app_state$zoom) && !is.null(current_app_state$assembly)) {
          
          # check which segments the current zoom spans using current context
          cdf <- tryCatch(cxt_get_entire_view(), error = function(e) NULL)
          
          if (!is.null(cdf) && nrow(cdf) > 0) {
            zoom_segments <- cdf[cdf$vstart < current_app_state$zoom[2] & cdf$vend > current_app_state$zoom[1], ]
            
            if (nrow(zoom_segments) >= 1) {
              # use contig of first segment that overlaps with zoom
              query_contig <- zoom_segments$contig[1]
              vstart <- zoom_segments$vstart[1]
              local_start_offset <- zoom_segments$start[1]
              
              # convert virtual coordinates to local coordinates
              query_start <- as.integer(max(1, round(current_app_state$zoom[1] - vstart) + local_start_offset))
              query_end <- as.integer(round(current_app_state$zoom[2] - vstart) + local_start_offset)
              
              # ensure end is at least start
              query_end <- max(query_end, query_start)
            } else {
              # fallback: use first segment's contig
              query_contig <- current_segments$contig[1]
            }
          } else {
            # fallback: use first segment's contig
            query_contig <- current_segments$contig[1]
          }
        }
        
        shiny::showModal(shiny::modalDialog(
          title = "Find Homologous Regions",
          size = "m",
          shiny::fluidRow(
            shiny::column(6,
              shiny::h5("Query Region"),
              shiny::textInput(ns("homologs_contig"), "Contig:", value = query_contig),
              shiny::numericInput(ns("homologs_start"), "Start Position:", value = query_start, min = 1),
              shiny::numericInput(ns("homologs_end"), "End Position:", value = query_end, min = 1)
            ),
            shiny::column(6,
              shiny::h5("Search Parameters"),
              shiny::numericInput(ns("homologs_k"), "Kmer Size:", value = 21, min = 1, max = 50),
              shiny::numericInput(ns("homologs_num_kmers"), "Number of Kmers:", value = 10, min = 1, max = 100),
              shiny::numericInput(ns("homologs_threshold"), "Threshold (%):", value = 80, min = 1, max = 100),
              shiny::numericInput(ns("homologs_margin"), "Margin (bp):", value = 5000, min = 0, step = 1000)
            )
          ),
          shiny::fluidRow(
            shiny::column(12,
              shiny::h5("Output"),
              shiny::textInput(ns("homologs_filename"), "Output Filename:", 
                       value = paste0("homologs_", query_contig, "_", query_start, "_", query_end, ".txt"),
                       placeholder = "homologs_results.txt")
            )
          ),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("homologs_run"), "Search Homologs", class = "btn-primary")
          ),
          easyClose = TRUE
        ))
      })

      # handle homologs search execution
      observeEvent(input$homologs_run, {
        tryCatch({
          # validate inputs
          if (is.null(input$homologs_contig) || input$homologs_contig == "") {
            shiny::showNotification("Please specify a query contig", type = "error")
            return()
          }
          
          if (is.null(input$homologs_filename) || input$homologs_filename == "") {
            shiny::showNotification("Please specify an output filename", type = "error")
            return()
          }
          
          # get assembly sequences
          current_assembly <- shiny::isolate(main_state_rv$assembly)
          assembly_sequences <- tryCatch({
            get_fasta(current_assembly)
          }, error = function(e) {
            shiny::showNotification(paste("Failed to load assembly sequences:", e$message), type = "error")
            return(NULL)
          })
          
          if (is.null(assembly_sequences)) {
            shiny::showNotification("No assembly sequences available", type = "error")
            return()
          }
          
          # validate query contig exists in assembly
          if (!input$homologs_contig %in% names(assembly_sequences)) {
            shiny::showNotification(paste("Contig", input$homologs_contig, "not found in assembly"), type = "error")
            return()
          }
          
          # validate query region
          contig_length <- nchar(assembly_sequences[[input$homologs_contig]])
          if (input$homologs_start > contig_length || input$homologs_end > contig_length) {
            shiny::showNotification(paste("Query region exceeds contig length:", contig_length), type = "error")
            return()
          }
          
          if (input$homologs_start > input$homologs_end) {
            shiny::showNotification("Start position must be <= end position", type = "error")
            return()
          }
          
          # show progress notification
          shiny::showNotification("Running homologs search...", type = "message", duration = 3)
          
          # call alntools homologs_search function (available after init_alntools)
          homologs_results <- tryCatch({
            # check if homologs_search function is available
            if (!exists("homologs_search", mode = "function")) {
              stop("homologs_search function not available - alntools may not be initialized")
            }
            
            # call homologs_search function directly (no package prefix needed)
            homologs_search(
              assembly_sequences = assembly_sequences,
              query_contig = input$homologs_contig,
              query_start = as.integer(input$homologs_start),
              query_end = as.integer(input$homologs_end),
              k = as.integer(input$homologs_k),
              num_kmers = as.integer(input$homologs_num_kmers),
              threshold = as.numeric(input$homologs_threshold),
              num_threads = 0  # use auto-detection for optimal performance
            )
          }, error = function(e) {
            # show error with helpful message
            shiny::showNotification(paste("Homologs search failed:", e$message), 
                                   type = "error", duration = 10)
            return(NULL)
          })
          
          if (is.null(homologs_results) || nrow(homologs_results) == 0) {
            shiny::showNotification("No homologous regions found", type = "warning")
            return()
          }
          
          # save results to regions directory
          regions_dir <- tryCatch(ensure_regions_dir(), error = function(e) {
            shiny::showNotification("Regions directory not configured", type = "error")
            return(NULL)
          })
          
          if (is.null(regions_dir)) return()
          
          # ensure filename has .txt extension
          filename <- input$homologs_filename
          if (!grepl("\\.txt$", filename)) {
            filename <- paste0(filename, ".txt")
          }
          
          file_path <- file.path(regions_dir, filename)
          
          # apply margin to extend regions
          margin <- as.integer(input$homologs_margin %||% 0)
          margin_suffix <- if (margin > 0) paste0("_margin_", format_bp(margin)) else ""
          
          # create better descriptions with length, kmer count, and margin info
          better_descriptions <- paste0("length_", homologs_results$length, 
                                       "_kmers_", homologs_results$kmer_count,
                                       "_coverage_", round(homologs_results$coverage_percent), "pct",
                                       margin_suffix)
          extended_starts <- pmax(1, homologs_results$start - margin)
          extended_ends <- homologs_results$end + margin
          
          # validate margin doesn't extend beyond contig boundaries
          for (i in seq_len(nrow(homologs_results))) {
            contig_name <- homologs_results$contig[i]
            if (contig_name %in% names(assembly_sequences)) {
              contig_length <- nchar(assembly_sequences[[contig_name]])
              extended_ends[i] <- min(extended_ends[i], contig_length)
            }
          }
          
          # convert homologs results to proper region table format
          region_format_results <- data.frame(
            id = homologs_results$id,
            level = rep(1L, nrow(homologs_results)),  # default level
            description = better_descriptions,
            assembly = rep(current_assembly, nrow(homologs_results)),  # use actual assembly name
            contigs = homologs_results$contig,  # single contig per region
            zoom_start = as.numeric(extended_starts),  # use extended coordinates
            zoom_end = as.numeric(extended_ends),      # use extended coordinates
            segment_contig = homologs_results$contig,
            segment_start = as.integer(extended_starts),
            segment_end = as.integer(extended_ends),
            single_contig = rep(TRUE, nrow(homologs_results)),
            stringsAsFactors = FALSE
          )
          
          # compute global zoom coordinates from segment coordinates
          for (i in seq_len(nrow(region_format_results))) {
            # use segment coordinates directly as zoom (they're already in the right space)
            # Note: With segment-based state, zoom coords are relative to the segment set
            region_format_results$zoom_start[i] <- region_format_results$segment_start[i]
            region_format_results$zoom_end[i] <- region_format_results$segment_end[i]
          }
          
          # write results to file
          tryCatch({
            write.table(region_format_results, 
                       file = file_path, 
                       sep = "\t", 
                       row.names = FALSE, 
                       col.names = TRUE, 
                       quote = FALSE)
            
            shiny::showNotification(paste("Found", nrow(homologs_results), "homologous regions. Results saved to:", filename), 
                                   type = "message", duration = 8)
          }, error = function(e) {
            shiny::showNotification(paste("Failed to save results:", e$message), type = "error")
          })
          
          # close dialog
          shiny::removeModal()
          
        }, error = function(e) {
          shiny::showNotification(paste("Homologs search failed:", e$message), type = "error", duration = 10)
        })
      })

      # save current region
      observeEvent(input$save_region, {
        if (current_regions_file() == "") {
          shiny::showNotification("no regions table selected", type = "error")
          return()
        }
        
        # get available region files
        regions_dir <- tryCatch(get_regions_dir(), error = function(e) {
          shiny::showNotification("regions directory not configured", type = "error")
          return(NULL)
        })
        
        if (is.null(regions_dir)) return()
        
        available_files <- list.files(regions_dir, pattern = "\\.txt$", full.names = FALSE)
        if (length(available_files) == 0) {
          available_files <- current_regions_file()
        }
        
        # determine default selection - use last used for adding, fallback to current
        default_selection <- if (cache_exists("regions.last_used_for_add")) {
          last_used <- cache_get("regions.last_used_for_add")
          if (last_used %in% available_files) last_used else current_regions_file()
        } else {
          current_regions_file()
        }
        
        # suggest next ID based on row count
        next_id <- as.character(nrow(region_table()) + 1)
        
        shiny::showModal(shiny::modalDialog(
          title = "Save Current Region",
          shiny::selectInput(ns("save_region_file"), "Save to Region File:", 
                           choices = available_files, 
                           selected = default_selection),
          shiny::textInput(ns("save_region_id"), "Region ID:", value = next_id, placeholder = "Enter region ID"),
          shiny::numericInput(ns("save_region_level"), "Level:", value = 1, min = 1, max = 10, step = 1),
          shiny::textInput(ns("save_region_description"), "Description:", placeholder = "Enter region description"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_save_region"), "Save")
          ),
          easyClose = TRUE
        ))
      })

      # update region ID when file selection changes in add dialog
      observeEvent(input$save_region_file, {
        req(input$save_region_file)
        
        # get next ID for the selected file
        selected_file <- input$save_region_file
        
        target_table <- if (selected_file == current_regions_file()) {
          region_table()
        } else {
          tryCatch({
            regions_dir <- get_regions_dir()
            file_path <- file.path(regions_dir, selected_file)
            if (file.exists(file_path)) {
              read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            } else {
              create_empty_region_table()
            }
          }, error = function(e) {
            create_empty_region_table()
          })
        }
        
        # calculate next ID based on row count
        next_id <- as.character(nrow(target_table) + 1)
        
        # update the ID input
        shiny::updateTextInput(session, "save_region_id", value = next_id)
      })

      observeEvent(input$confirm_save_region, {
        req(input$save_region_description)
        req(input$save_region_id)
        req(input$save_region_level)
        req(input$save_region_file)
        description <- trimws(input$save_region_description)
        custom_id <- trimws(input$save_region_id)
        level <- as.integer(input$save_region_level)
        selected_file <- input$save_region_file
        
        if (description == "") {
          shiny::showNotification("description cannot be empty", type = "error")
          return()
        }
        
        if (custom_id == "") {
          shiny::showNotification("ID cannot be empty", type = "error")
          return()
        }
        
        if (is.na(level) || level < 1 || level > 10) {
          shiny::showNotification("level must be between 1 and 10", type = "error")
          return()
        }
        
        # load target file to check for duplicate IDs
        target_table <- if (selected_file == current_regions_file()) {
          region_table()
        } else {
          tryCatch({
            regions_dir <- get_regions_dir()
            file_path <- file.path(regions_dir, selected_file)
            if (file.exists(file_path)) {
              read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            } else {
              create_empty_region_table()
            }
          }, error = function(e) {
            shiny::showNotification(paste("error loading target file:", e$message), type = "error")
            return()
          })
        }
        
        # check for duplicate ID in target file
        if (custom_id %in% target_table$id) {
          shiny::showNotification(paste("ID", custom_id, "already exists in", selected_file), type = "error")
          return()
        }
        
        current_app_state <- shiny::isolate(shiny::reactiveValuesToList(main_state_rv))
        
        # prepare contigs (comma-separated)
        current_contigs <- unique(current_app_state$segments$contig)
        contigs_str <- if (is.null(current_contigs) || length(current_contigs) == 0) {
          ""
        } else {
          paste(current_contigs, collapse = ",")
        }
        
        # prepare zoom
        zoom_start <- if (is.null(current_app_state$zoom)) NA else current_app_state$zoom[1]
        zoom_end <- if (is.null(current_app_state$zoom)) NA else current_app_state$zoom[2]
        
        new_row <- data.frame(
          id = custom_id,
          level = level,
          description = description,
          assembly = current_app_state$assembly %||% "",
          contigs = contigs_str,
          zoom_start = zoom_start,
          zoom_end = zoom_end,
          segment_contig = NA_character_,
          segment_start = 0L,
          segment_end = 0L,
          single_contig = FALSE,
          stringsAsFactors = FALSE
        )
        
        # compute segment information
        segment_info <- compute_segment_info(new_row)
        new_row$segment_contig <- segment_info$contig
        new_row$segment_start <- segment_info$start
        new_row$segment_end <- segment_info$end
        new_row$single_contig <- segment_info$single_contig
        
        # add region to target table and save
        updated_target_table <- rbind(target_table, new_row)
        
        # save to selected file
        regions_dir <- get_regions_dir()
        file_path <- file.path(regions_dir, selected_file)
        
        # create backup version if file exists
        if (file.exists(file_path)) {
          versions_dir <- file.path(regions_dir, "versions")
          if (!dir.exists(versions_dir)) {
            dir.create(versions_dir, recursive = TRUE)
          }
          base_name <- sub("\\.txt$", "", selected_file)
          existing_versions <- list.files(versions_dir, pattern = paste0("^", base_name, "_v\\d+\\.txt$"))
          version_num <- length(existing_versions) + 1
          version_filename <- paste0(base_name, "_v", version_num, ".txt")
          file.copy(file_path, file.path(versions_dir, version_filename))
        }
        
        write.table(updated_target_table, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        
        # cache the file as last used for adding
        cache_set("regions.last_used_for_add", selected_file)
        
        # update current table if we saved to the current file
        if (selected_file == current_regions_file()) {
          region_table(updated_target_table)
        }
        
        shiny::removeModal()
        shiny::showNotification(paste("region '", description, "' added to", selected_file), type = "message")
      })

      # edit region description
      observeEvent(input$edit_region, {
        selected_rows <- input$region_table_display_rows_selected
        if (length(selected_rows) == 0) {
          shiny::showNotification("please select a region to edit", type = "warning")
          return()
        }
        
        row_data <- region_table()[selected_rows[1], ]
        
        shiny::showModal(shiny::modalDialog(
          title = "Edit Region",
          shiny::textInput(ns("edit_id_input"), "Region ID:", value = row_data$id, placeholder = "Enter region ID"),
          shiny::numericInput(ns("edit_level_input"), "Level:", value = if("level" %in% colnames(row_data)) row_data$level else 1, min = 1, max = 10, step = 1),
          shiny::textInput(ns("edit_description_input"), "Description:", value = row_data$description),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_edit_region"), "Update")
          ),
          easyClose = TRUE
        ))
      })

      observeEvent(input$confirm_edit_region, {
        req(input$edit_description_input)
        req(input$edit_id_input)
        req(input$edit_level_input)
        selected_rows <- input$region_table_display_rows_selected
        if (length(selected_rows) == 0) return()
        
        new_description <- trimws(input$edit_description_input)
        new_id <- trimws(input$edit_id_input)
        new_level <- as.integer(input$edit_level_input)
        
        if (new_description == "") {
          shiny::showNotification("description cannot be empty", type = "error")
          return()
        }
        
        if (new_id == "") {
          shiny::showNotification("ID cannot be empty", type = "error")
          return()
        }
        
        if (is.na(new_level) || new_level < 1 || new_level > 10) {
          shiny::showNotification("level must be between 1 and 10", type = "error")
          return()
        }
        
        current_table <- region_table()
        old_id <- current_table[selected_rows[1], "id"]
        
        # check for duplicate ID (but allow keeping the same ID)
        if (new_id != old_id && new_id %in% current_table$id) {
          shiny::showNotification(paste("ID", new_id, "already exists"), type = "error")
          return()
        }
        
        current_table[selected_rows[1], "id"] <- new_id
        current_table[selected_rows[1], "level"] <- new_level
        current_table[selected_rows[1], "description"] <- new_description
        region_table(current_table)
        save_region_table(current_regions_file())
        
        shiny::removeModal()
        shiny::showNotification("region updated", type = "message")
      })

      # delete region
      observeEvent(input$delete_region, {
        selected_rows <- input$region_table_display_rows_selected
        if (length(selected_rows) == 0) {
          shiny::showNotification("please select a region to delete", type = "warning")
          return()
        }
        
        row_to_delete_id <- region_table()[selected_rows[1], "id"]
        description <- region_table()[selected_rows[1], "description"]
        
        shiny::showModal(shiny::modalDialog(
          title = "Delete Region",
          paste("are you sure you want to delete region:", description, "?"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_delete_region"), "Delete")
          )
        ))
      })

      observeEvent(input$confirm_delete_region, {
        selected_rows <- input$region_table_display_rows_selected
        if (length(selected_rows) == 0) return()
        
        current_table <- region_table()
        region_table(current_table[-selected_rows[1], ])
        save_region_table(current_regions_file())
        
        shiny::removeModal()
        shiny::showNotification("region deleted", type = "message")
      })

      # goto region (with undo functionality)
      goto_region <- function(region_row) {
        # push current region to undo stack
        current_state <- list(
          contigs = main_state_rv$contigs,
          zoom = main_state_rv$zoom,
          assembly = main_state_rv$assembly
        )
        
        current_undo <- undo_stack()
        current_undo[[length(current_undo) + 1]] <- current_state
        undo_stack(current_undo)
        
        # apply region
        main_state_rv$assembly <- region_row$assembly
        
        # parse contigs and get segments
        if (region_row$contigs != "" && !is.na(region_row$contigs)) {
          contig_list <- trimws(strsplit(region_row$contigs, ",")[[1]])
          segments <- get_segments(main_state_rv$assembly)
          selected_segments <- segments[segments$contig %in% contig_list, ]
          main_state_rv$segments <- selected_segments
        } else {
          main_state_rv$segments <- data.frame(segment = character(), contig = character(), 
                                                start = integer(), end = integer(), stringsAsFactors = FALSE)
        }
        
        # apply zoom
        if (!is.na(region_row$zoom_start) && !is.na(region_row$zoom_end)) {
          main_state_rv$zoom <- c(region_row$zoom_start, region_row$zoom_end)
        } else {
          main_state_rv$zoom <- NULL
        }
        
        shiny::showNotification(paste("navigated to region:", region_row$description), type = "message")
      }

      # push current region to undo stack (for external use)
      push_undo_state <- function() {
        current_state <- list(
          contigs = main_state_rv$contigs,
          zoom = main_state_rv$zoom,
          assembly = main_state_rv$assembly
        )
        
        current_undo <- undo_stack()
        current_undo[[length(current_undo) + 1]] <- current_state
        undo_stack(current_undo)
      }

      # undo functionality
      undo_last_action <- function() {
        current_undo <- undo_stack()
        if (length(current_undo) == 0) {
          shiny::showNotification("no actions to undo", type = "warning")
          return()
        }
        
        # restore last region
        last_state <- current_undo[[length(current_undo)]]
        main_state_rv$segments <- last_state$segments
        main_state_rv$zoom <- last_state$zoom
        main_state_rv$assembly <- last_state$assembly
        
        # remove from undo stack
        undo_stack(current_undo[-length(current_undo)])
        shiny::showNotification("undid last action", type = "message")
      }

      # --- Table Display ---
      output$region_table_display <- DT::renderDT({
        rt <- region_table()
        if (nrow(rt) == 0) return(rt)
        
        # create display table with formatted columns
        display_table <- data.frame(
          ID = rt$id,
          Level = if("level" %in% colnames(rt)) rt$level else rep(1L, nrow(rt)),
          Description = rt$description,
          Assembly = rt$assembly,
          Contigs = sapply(rt$contigs, function(x) {
            if (x == "" || is.na(x)) return("None")
            contigs_list <- trimws(strsplit(x, ",")[[1]])
            if (length(contigs_list) == 1) return(contigs_list[1])
            return(paste(length(contigs_list), "contigs"))
          }),
          Window = sapply(seq_len(nrow(rt)), function(i) {
            if (is.na(rt$zoom_start[i]) || is.na(rt$zoom_end[i])) return("Full range")
            window_size <- round(rt$zoom_end[i]) - round(rt$zoom_start[i])
            return(format_bp(window_size))
          }),
          Coordinates = sapply(seq_len(nrow(rt)), function(i) {
            if (is.na(rt$zoom_start[i]) || is.na(rt$zoom_end[i])) return("Full range")
            
            contigs_list <- if (rt$contigs[i] == "" || is.na(rt$contigs[i])) {
              character(0)
            } else {
              trimws(strsplit(rt$contigs[i], ",")[[1]])
            }
            
            if (length(contigs_list) == 0) return("Global")
            
            # For segment-based state, just show zoom coords directly
            # (more complex conversion would require temporarily setting context)
            if (length(contigs_list) == 1) {
              sprintf("%s: %d-%d", contigs_list[1], rt$zoom_start[i], rt$zoom_end[i])
            } else if (length(contigs_list) > 1) {
              "Multiple contigs"
            } else {
              "Global"
            }
          }),
          stringsAsFactors = FALSE
        )
        
        # add goto button column with proper actionButton approach
        display_table$Actions <- sapply(seq_len(nrow(display_table)), function(i) {
          as.character(shiny::actionButton(
            paste0("goto_region_btn_", i), 
            "Goto",
            class = "btn btn-primary btn-sm",
            onclick = sprintf("Shiny.setInputValue('%s', %d, {priority: 'event'})", ns('goto_region_trigger'), i)
          ))
        })
        
        DT::datatable(
          display_table,
          selection = "single",
          escape = FALSE,
          rownames = FALSE,
          options = list(
            pageLength = 20,
            lengthMenu = c(20, 50, 100),
            dom = "ltip",
            columnDefs = list(
              list(targets = ncol(display_table) - 1, orderable = FALSE, width = "80px")
            )
          )
        )
      })

      # goto region button (standalone button)
      observeEvent(input$goto_region_btn, {
        selected_rows <- input$region_table_display_rows_selected
        if (length(selected_rows) == 0) {
          shiny::showNotification("please select a region to go to", type = "warning")
          return()
        }
        
        region_row <- region_table()[selected_rows[1], ]
        goto_region(region_row)
      })

      # handle goto button clicks (table buttons)
      observeEvent(input$goto_region_trigger, {
        req(input$goto_region_trigger)
        row_idx <- as.integer(input$goto_region_trigger)
        rt <- region_table()
        
        # ensure row index is valid
        if (row_idx > 0 && row_idx <= nrow(rt)) {
          region_row <- rt[row_idx, ]
          goto_region(region_row)
        }
      })

      # trigger save region dialog (for keyboard shortcut)
      focus_save_input <- function() {
        if (current_regions_file() == "") {
          shiny::showNotification("no regions table selected", type = "error")
          return()
        }
        
        # get available region files
        regions_dir <- tryCatch(get_regions_dir(), error = function(e) {
          shiny::showNotification("regions directory not configured", type = "error")
          return(NULL)
        })
        
        if (is.null(regions_dir)) return()
        
        available_files <- list.files(regions_dir, pattern = "\\.txt$", full.names = FALSE)
        if (length(available_files) == 0) {
          available_files <- current_regions_file()
        }
        
        # determine default selection - use last used for adding, fallback to current
        default_selection <- if (cache_exists("regions.last_used_for_add")) {
          last_used <- cache_get("regions.last_used_for_add")
          if (last_used %in% available_files) last_used else current_regions_file()
        } else {
          current_regions_file()
        }
        
        # suggest next ID based on row count
        next_id <- as.character(nrow(region_table()) + 1)
        
        shiny::showModal(shiny::modalDialog(
          title = "Save Current Region",
          shiny::selectInput(ns("save_region_file"), "Save to Region File:", 
                           choices = available_files, 
                           selected = default_selection),
          shiny::textInput(ns("save_region_id"), "Region ID:", value = next_id, placeholder = "Enter region ID"),
          shiny::numericInput(ns("save_region_level"), "Level:", value = 1, min = 1, max = 10, step = 1),
          shiny::textInput(ns("save_region_description"), "Description:", placeholder = "Enter region description"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_save_region"), "Save")
          ),
          easyClose = TRUE
        ))
      }

      # get current region table data (for external use)
      get_region_table <- function() {
        return(region_table())
      }

      # return functions for external use
      return(list(
        undo_last_action = undo_last_action,
        push_undo_state = push_undo_state,
        goto_region = goto_region,
        focus_save_input = focus_save_input,
        get_region_table = get_region_table
      ))
    }
  )
}
