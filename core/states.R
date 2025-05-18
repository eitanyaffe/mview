# core/states.R
# Defines the UI and server logic for the state management module.

# Hardcoded fields from the main 'state' reactiveValues to be tracked
TRACKED_STATE_FIELDS <- c("contigs", "zoom")

# ---- UI Definition ----

states_ui <- function(id = "states_module") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::actionButton(ns("new_state_table"), "New Table"),
        shiny::actionButton(ns("load_state_table_modal_trigger"), "Load Table"),
        shiny::actionButton(ns("save_state_modal_trigger"), "Save Table As..."),
        shiny::hr(),
        # UI for adding state with inline text input
        shiny::div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          shiny::textInput(ns("add_state_description_input"), label = NULL, placeholder = "Enter state description", width = "300px"),
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
    # Removed uiOutput(ns("file_input_ui")) as it was a placeholder
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
          .Contigs_Data = I(list()),
          .Zoom_Data = I(list()),
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
        return(paste(zoom_vector[1], "–", zoom_vector[2]))
      }

      # --- Reactive Values for the states module ---
      state_table <- shiny::reactiveVal(create_empty_state_table())
      undo_stack <- shiny::reactiveVal(list())
      current_save_filename <- shiny::reactiveVal("mview_states_default.state")

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
          .Contigs_Data = I(list(current_app_state$contigs)),
          .Zoom_Data = I(list(current_app_state$zoom)),
          stringsAsFactors = FALSE
        )
        state_table(rbind(state_table(), new_row))
        shiny::showNotification(paste("State '", description, "' added."), type = "message")
      })

      observeEvent(input$delete_state, {
        selected_rows <- input$state_table_display_rows_selected
        if (length(selected_rows) == 0) {
          shiny::showModal(shiny::modalDialog(title = "Delete State", "Please select a state to delete.", footer = shiny::modalButton("OK")))
          return()
        }
        row_to_delete_id <- state_table()[selected_rows[1], "ID"]
        shiny::showModal(shiny::modalDialog(
          title = "Delete State",
          paste("Are you sure you want to delete state: ", state_table()[state_table()$ID == row_to_delete_id, "Description"], "?"),
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
          shiny::showModal(shiny::modalDialog(title = "Go to State", "Please select a state to go to.", footer = shiny::modalButton("OK")))
          return()
        }
        push_state()
        selected_row_data <- state_table()[selected_rows[1], ]
        main_state_rv$contigs <- selected_row_data$.Contigs_Data[[1]]
        main_state_rv$zoom <- selected_row_data$.Zoom_Data[[1]]
        shiny::showNotification(paste("Went to state: ", selected_row_data$Description), type = "message")
      })

      # --- File Load Operations ---
      observeEvent(input$load_state_table_modal_trigger, {
        shiny::showModal(shiny::modalDialog(
          title = "Load State Table",
          shiny::fileInput(ns("load_file_input"), "Choose .state File", accept = c(".state")),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_load_state_table"), "Load")
          )
        ))
      })

      observeEvent(input$confirm_load_state_table, {
        req(input$load_file_input)
        file_path <- input$load_file_input$datapath
        tryCatch(
          {
            loaded_data <- readRDS(file_path)
            if (!is.data.frame(loaded_data) ||
              !all(c("ID", "Description", "Contigs", "Zoom", ".Contigs_Data", ".Zoom_Data") %in% names(loaded_data))) {
              stop("Invalid state file format.")
            }
            state_table(loaded_data)
            shiny::removeModal()
            shiny::showNotification(paste("State table loaded from", basename(input$load_file_input$name)), type = "message")
          },
          error = function(e) {
            shiny::removeModal()
            shiny::showModal(shiny::modalDialog(
              title = "Error Loading File",
              paste("Could not load state table:", e$message),
              footer = shiny::modalButton("OK")
            ))
          }
        )
      })

      # --- File Save Operations ---
      observeEvent(input$save_state_modal_trigger, {
        shiny::showModal(shiny::modalDialog(
          title = "Save State Table As...",
          shiny::textInput(ns("save_filename_input"),
            label = "Filename:",
            value = isolate(current_save_filename())
          ),
          shiny::HTML("<p><small>The table will be saved to your browser's default download location. <br>Ensure the name ends with <code>.state</code> (it will be added if missing).</small></p>"),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::downloadButton(ns("execute_save_from_modal"), "Save")
          ),
          easyClose = TRUE
        ))
      })

      output$execute_save_from_modal <- shiny::downloadHandler(
        filename = function() {
          req(input$save_filename_input)
          fname <- input$save_filename_input
          if (is.null(fname) || trimws(fname) == "") {
            fname <- "mview_states_unnamed.state"
          }
          if (!grepl("\\.state$", fname, ignore.case = TRUE)) {
            fname <- paste0(fname, ".state")
          }
          current_save_filename(fname)
          return(basename(fname))
        },
        content = function(file) {
          tryCatch({
            saveRDS(state_table(), file)
            shiny::showNotification("State table saved.", type = "message")
          }, error = function(e) {
            shiny::showNotification(paste("Error saving state table:", e$message), type = "error")
          }, finally = {
            shiny::removeModal()
          })
        }
      )

      # --- Table Display ---
      output$state_table_display <- DT::renderDT({
        st <- state_table()
        display_table <- st[, c("ID", "Description", "Contigs", "Zoom"), drop = FALSE]
        DT::datatable(
          display_table,
          rownames = FALSE,
          selection = "single",
          options = list(pageLength = 5, lengthMenu = c(5, 10, 20))
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
