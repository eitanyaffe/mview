# core/states.R
# Defines the UI and server logic for the state management module.

# Hardcoded fields from the main 'state' reactiveValues to be tracked
TRACKED_STATE_FIELDS <- c("contigs", "zoom", "assembly")

# ---- UI Definition ----

states_default_fn <- paste0(get_project_id() %||% "states1", ".mv")

states_ui <- function(id = "states_module") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::textInput(ns("states_file_input"), "States File:", value = states_default_fn, width = "200px"),
        shiny::actionButton(ns("select_table_trigger"), "Select Table"),
        shiny::actionButton(ns("new_state_table"), "New Table"),
        shiny::actionButton(ns("load_state_table_modal_trigger"), "Load Table"),
        shiny::actionButton(ns("save_table_trigger"), "Save Table"),
        shiny::hr(),
        # UI for adding state with inline text input
        shiny::div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          shiny::textInput(ns("add_state_description_input"),
            label = NULL,
            placeholder = "Enter state description", width = "300px"
          ),
          shiny::actionButton(ns("confirm_add_state_from_input"), "Add Current State", style = "margin-left: 10px;")
        ),
        shiny::actionButton(ns("delete_state"), "Delete Selected State"),
        shiny::actionButton(ns("goto_state"), "Go to Selected State")
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::h4("Saved States"),
        DT::DTOutput(ns("state_table_display"))
      )
    )
  )
}

# ---- Server Logic ----

states_server <- function(id = "states_module", main_state_rv, session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # --- Helper Functions (defined first) ---
      create_empty_state_table <- function() {
        df <- data.frame(
          ID = integer(),
          Description = character(),
          Contigs = character(),
          Zoom = character(),
          Assembly = character(),
          .Contigs_Data = I(list()),
          .Zoom_Data = I(list()),
          .Assembly_Data = I(list()),
          stringsAsFactors = FALSE
        )
        return(df)
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

      format_assembly_for_display <- function(assembly_value) {
        if (is.null(assembly_value) || assembly_value == "") {
          return("None")
        }
        return(as.character(assembly_value))
      }

      # --- Reactive Values for the states module ---
      state_table <- shiny::reactiveVal(create_empty_state_table())
      undo_stack <- shiny::reactiveVal(list())
      current_save_filename <- shiny::reactiveVal("states1.mv")

      # Ensure states directory exists
      ensure_states_dir <- function() {
        states_dir <- "states"
        if (!dir.exists(states_dir)) {
          dir.create(states_dir)
        }
        return(states_dir)
      }

      # Observer for assembly selection changes
      observeEvent(input$assembly_select,
        {
          req(input$assembly_select)
          # Do not trigger on initial NULL or empty string selection from selectInput clearing
          if (is.null(main_state_rv$assembly) || shiny::isolate(main_state_rv$assembly) != input$assembly_select) {
            push_state() # Push current state before changing assembly
            main_state_rv$assembly <- input$assembly_select
            shiny::showNotification(paste("Assembly changed to:", input$assembly_select), type = "message")
          }
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
      ) # ignoreNULL=FALSE to react to deselection (empty string)

      # Update assembly dropdown when main_state_rv$assembly changes (e.g., on Go to State)
      observeEvent(main_state_rv$assembly,
        {
          current_selection <- shiny::isolate(input$assembly_select)
          if (!is.null(main_state_rv$assembly) && main_state_rv$assembly != "" &&
            main_state_rv$assembly != current_selection) {
            shiny::updateSelectInput(session, "assembly_select", selected = main_state_rv$assembly)
          } else if ((is.null(main_state_rv$assembly) || main_state_rv$assembly == "") &&
            !is.null(current_selection) && current_selection != "") {
            # If main state assembly becomes NULL (e.g. new table), clear selection
            shiny::updateSelectInput(session, "assembly_select", selected = "")
          }
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
      )

      # --- Direct Push/Undo Functions ---
      push_state <- function() {
        current_app_state <- shiny::isolate(shiny::reactiveValuesToList(main_state_rv))
        state_to_push <- current_app_state[TRACKED_STATE_FIELDS]
        current_stack <- undo_stack()
        print_state(state_to_push, "Push")
        current_stack <- c(list(state_to_push), current_stack)
        undo_stack(current_stack)
      }

      undo_state <- function() {
        current_stack <- undo_stack()
        if (length(current_stack) > 0) {
          state_to_restore <- current_stack[[1]]
          main_state_rv$contigs <- state_to_restore$contigs
          main_state_rv$zoom <- state_to_restore$zoom
          main_state_rv$assembly <- state_to_restore$assembly
          print_state(state_to_restore, "Undo")
          undo_stack(current_stack[-1])
          shiny::showNotification("Undo successful.", type = "message")
          return(TRUE)
        } else {
          shiny::showNotification("Nothing to undo.", type = "warning")
          return(FALSE)
        }
      }

      # --- Core State Management Functions ---
      observeEvent(input$new_state_table, {
        shiny::showModal(shiny::modalDialog(
          title = "New State Table",
          "Are you sure you want to clear all saved states? This action cannot be undone.",
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_new_state_table"), "Clear Table")
          )
        ))
      })

      observeEvent(input$confirm_new_state_table, {
        state_table(create_empty_state_table())
        main_state_rv$assembly <- NULL # Clear assembly in main state
        shiny::removeModal()
        shiny::showNotification("New state table created.", type = "message")
      })

      # Observer to update the default description for adding a new state
      observe({
        current_table <- state_table()
        default_description <- paste("State", nrow(current_table) + 1)
        shiny::updateTextInput(session, "add_state_description_input", value = default_description)
      })

      # Add current state (using inline text input)
      observeEvent(input$confirm_add_state_from_input, {
        req(input$add_state_description_input)
        description <- input$add_state_description_input
        current_app_state <- shiny::isolate(shiny::reactiveValuesToList(main_state_rv))
        new_row <- data.frame(
          ID = if (nrow(state_table()) == 0) 1 else max(state_table()$ID, 0) + 1,
          Description = description,
          Contigs = format_contigs_for_display(current_app_state$contigs),
          Zoom = format_zoom_for_display(current_app_state$zoom),
          Assembly = format_assembly_for_display(current_app_state$assembly),
          .Contigs_Data = I(list(current_app_state$contigs)),
          .Zoom_Data = I(list(current_app_state$zoom)),
          .Assembly_Data = I(list(current_app_state$assembly)),
          stringsAsFactors = FALSE
        )
        state_table(rbind(state_table(), new_row))
        shiny::showNotification(paste("State '", description, "' added."), type = "message")
      })

      observeEvent(input$delete_state, {
        selected_rows <- input$state_table_display_rows_selected
        if (length(selected_rows) == 0) {
          shiny::showModal(shiny::modalDialog(
            title = "Delete State",
            "Please select a state to delete.", footer = shiny::modalButton("OK")
          ))
          return()
        }
        row_to_delete_id <- state_table()[selected_rows[1], "ID"]
        shiny::showModal(shiny::modalDialog(
          title = "Delete State",
          paste("Are you sure you want to delete state: ", state_table()[
            state_table()$ID == row_to_delete_id,
            "Description"
          ], "?"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_delete_state"), "Delete")
          )
        ))
      })

      observeEvent(input$confirm_delete_state, {
        selected_rows <- input$state_table_display_rows_selected
        req(length(selected_rows) > 0)
        row_to_delete_id <- state_table()[selected_rows[1], "ID"]
        current_table <- state_table()
        current_table <- current_table[current_table$ID != row_to_delete_id, ]
        state_table(current_table)
        shiny::removeModal()
        shiny::showNotification("State deleted.", type = "message")
      })

      observeEvent(input$goto_state, {
        selected_rows <- input$state_table_display_rows_selected
        if (length(selected_rows) == 0) {
          shiny::showModal(shiny::modalDialog(
            title = "Go to State", "Please select a state to go to.",
            footer = shiny::modalButton("OK")
          ))
          return()
        }
        push_state()
        selected_row_data <- state_table()[selected_rows[1], ]
        main_state_rv$contigs <- selected_row_data$.Contigs_Data[[1]]
        main_state_rv$zoom <- selected_row_data$.Zoom_Data[[1]]
        main_state_rv$assembly <- selected_row_data$.Assembly_Data[[1]]
        shiny::showNotification(paste("Went to state: ", selected_row_data$Description), type = "message")
      })

      # --- File Load Operations ---
      observeEvent(input$load_state_table_modal_trigger, {
        states_dir <- ensure_states_dir()
        file_path <- file.path(states_dir, current_save_filename())
        if (file.exists(file_path)) {
          tryCatch(
            {
              loaded_data <- readRDS(file_path)
              if (!is.data.frame(loaded_data) ||
                !all(c(
                  "ID", "Description", "Contigs", "Zoom", "Assembly", ".Contigs_Data",
                  ".Zoom_Data", ".Assembly_Data"
                ) %in% names(loaded_data))) {
                stop("Invalid state file format.")
              }
              state_table(loaded_data)
              shiny::showNotification(paste("State table loaded from", current_save_filename()), type = "message")
            },
            error = function(e) {
              shiny::showModal(shiny::modalDialog(
                title = "Error Loading File",
                paste("Could not load state table:", e$message),
                footer = shiny::modalButton("OK")
              ))
            }
          )
        } else {
          shiny::showModal(shiny::modalDialog(
            title = "File Not Found",
            paste("State table file", current_save_filename(), "not found in states directory."),
            footer = shiny::modalButton("OK")
          ))
        }
      })

      # --- File Save Operations ---
      observeEvent(input$select_table_trigger, {
        states_dir <- ensure_states_dir()
        files <- list.files(states_dir, pattern = "\\.mv$", full.names = FALSE)
        if (length(files) == 0) {
          shiny::showModal(shiny::modalDialog(
            title = "Select Table",
            "No .mv files found in states directory.",
            footer = shiny::modalButton("OK")
          ))
          return()
        }

        shiny::showModal(shiny::modalDialog(
          title = "Select Table",
          shiny::selectInput(ns("select_file_input"),
            label = "Choose a file:",
            choices = files,
            selected = isolate(current_save_filename())
          ),
          shiny::HTML("<p><small>Select a file from the states directory</small></p>"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_select_table"), "Select")
          ),
          easyClose = TRUE
        ))
      })

      observeEvent(input$confirm_select_table, {
        req(input$select_file_input)
        fname <- input$select_file_input
        current_save_filename(fname)
        shiny::updateTextInput(session, "states_file_input", value = fname)
        shiny::removeModal()
        shiny::showNotification(paste("Selected table:", fname), type = "message")
      })

      observeEvent(input$save_table_trigger, {
        states_dir <- ensure_states_dir()
        file_path <- file.path(states_dir, current_save_filename())
        tryCatch(
          {
            saveRDS(state_table(), file_path)
            shiny::showNotification(paste("State table saved to", current_save_filename()), type = "message")
          },
          error = function(e) {
            shiny::showNotification(paste("Error saving state table:", e$message), type = "error")
          }
        )
      })

      # --- Table Display ---
      output$state_table_display <- DT::renderDT({
        st <- state_table()
        # Reorder columns to show Assembly before Contigs
        display_table <- st[, c("ID", "Description", "Assembly", "Contigs", "Zoom"), drop = FALSE]
        DT::datatable(
          display_table,
          rownames = FALSE,
          selection = "single",
          options = list(pageLength = 15, lengthMenu = c(5, 10, 20))
        )
      })

      # Return the direct functions for use outside the module
      return(list(
        push_state = push_state,
        undo_state = undo_state
      ))
    }
  )
}
