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
        shiny::hr(),
        shiny::actionButton(ns("add_region"), "Add"),
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
          id = integer(),
          description = character(),
          assembly = character(),
          contigs = character(),
          zoom_start = numeric(),
          zoom_end = numeric(),
          stringsAsFactors = FALSE
        )
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
        
        region_table(df)
        set_current_regions_file(filename)
        shiny::showNotification(paste("loaded regions from", filename), type = "message")
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
        region_table(create_empty_region_table())
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

      # add current region
      observeEvent(input$add_region, {
        if (current_regions_file() == "") {
          shiny::showNotification("no regions table selected", type = "error")
          return()
        }
        
        shiny::showModal(shiny::modalDialog(
          title = "Add Current Region",
          shiny::textInput(ns("add_region_description"), "Description:", placeholder = "Enter region description"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_add_region"), "Add")
          ),
          easyClose = TRUE
        ))
      })

      observeEvent(input$confirm_add_region, {
        req(input$add_region_description)
        description <- trimws(input$add_region_description)
        
        if (description == "") {
          shiny::showNotification("description cannot be empty", type = "error")
          return()
        }
        
        current_app_state <- shiny::isolate(shiny::reactiveValuesToList(main_state_rv))
        
        # prepare contigs (comma-separated)
        contigs_str <- if (is.null(current_app_state$contigs) || length(current_app_state$contigs) == 0) {
          ""
        } else {
          paste(current_app_state$contigs, collapse = ",")
        }
        
        # prepare zoom
        zoom_start <- if (is.null(current_app_state$zoom)) NA else current_app_state$zoom[1]
        zoom_end <- if (is.null(current_app_state$zoom)) NA else current_app_state$zoom[2]
        
        new_row <- data.frame(
          id = if (nrow(region_table()) == 0) 1 else max(region_table()$id, 0) + 1,
          description = description,
          assembly = current_app_state$assembly %||% "",
          contigs = contigs_str,
          zoom_start = zoom_start,
          zoom_end = zoom_end,
          stringsAsFactors = FALSE
        )
        
        region_table(rbind(region_table(), new_row))
        save_region_table(current_regions_file())
        
        shiny::removeModal()
        shiny::showNotification(paste("region '", description, "' added."), type = "message")
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
        selected_rows <- input$region_table_display_rows_selected
        if (length(selected_rows) == 0) return()
        
        new_description <- trimws(input$edit_description_input)
        if (new_description == "") {
          shiny::showNotification("description cannot be empty", type = "error")
          return()
        }
        
        current_table <- region_table()
        current_table[selected_rows[1], "description"] <- new_description
        region_table(current_table)
        save_region_table(current_regions_file())
        
        shiny::removeModal()
        shiny::showNotification("region description updated", type = "message")
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
        
        # parse contigs
        if (region_row$contigs != "" && !is.na(region_row$contigs)) {
          main_state_rv$contigs <- trimws(strsplit(region_row$contigs, ",")[[1]])
        } else {
          main_state_rv$contigs <- NULL
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
        main_state_rv$contigs <- last_state$contigs
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
            
            cxt <- tryCatch({
              build_context(contigs_list, get_contigs(rt$assembly[i]), 
                          c(rt$zoom_start[i], rt$zoom_end[i]), rt$assembly[i])
            }, error = function(e) NULL)
            
            if (is.null(cxt) || is.null(cxt$mapper)) {
              return("Global")
            }
            
            cdf <- cxt$mapper$cdf
            zoom_contigs <- cdf[cdf$start < rt$zoom_end[i] & cdf$end > rt$zoom_start[i], ]
            
            if (nrow(zoom_contigs) > 1) {
              "Multiple contigs"
            } else if (nrow(zoom_contigs) == 1) {
              contig_name <- zoom_contigs$contig[1]
              local_start <- round(rt$zoom_start[i] - zoom_contigs$start[1]) + 1
              local_end <- round(rt$zoom_end[i] - zoom_contigs$start[1])
              sprintf("%s: %d-%d", contig_name, local_start, local_end)
            } else {
              "Global"
            }
          }),
          stringsAsFactors = FALSE
        )
        
        # add goto button column
        display_table$Actions <- paste0(
          '<button type="button" class="btn btn-primary btn-sm" onclick="Shiny.setInputValue(\'',
          ns('goto_region'), '\', ', rt$id, ', {priority: \'event\'});">Goto</button>'
        )
        
        DT::datatable(
          display_table,
          selection = "single",
          escape = FALSE,
          rownames = FALSE,
          options = list(
            pageLength = 10,
            dom = "tip",
            columnDefs = list(
              list(targets = ncol(display_table) - 1, orderable = FALSE, width = "80px")
            )
          )
        )
      })

      # handle goto button clicks
      observeEvent(input$goto_region, {
        region_id <- input$goto_region
        region_row <- region_table()[region_table()$id == region_id, ]
        if (nrow(region_row) > 0) {
          goto_region(region_row[1, ])
        }
      })

      # trigger add region dialog (for keyboard shortcut)
      focus_add_input <- function() {
        if (current_regions_file() == "") {
          shiny::showNotification("no regions table selected", type = "error")
          return()
        }
        
        shiny::showModal(shiny::modalDialog(
          title = "Add Current Region",
          shiny::textInput(ns("add_region_description"), "Description:", placeholder = "Enter region description"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_add_region"), "Add")
          ),
          easyClose = TRUE
        ))
      }

      # return functions for external use
      return(list(
        undo_last_action = undo_last_action,
        push_undo_state = push_undo_state,
        goto_region = goto_region,
        focus_add_input = focus_add_input
      ))
    }
  )
}
