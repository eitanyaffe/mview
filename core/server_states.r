# Server-side logic for states management

# create UI for state buttons
output$state_buttons_ui <- renderUI({
  # establish reactive dependency to refresh when states are saved
  state_buttons_trigger()
  
  state_buttons <- lapply(1:9, function(i) {
    title <- get_state_title(i)
    title_display <- if (nchar(title) > 0) title else "empty"
    
    div(
      style = "margin-bottom: 5px; display: flex; align-items: center;",
      actionButton(
        inputId = paste0("load_state_", i),
        label = paste0("L", i),
        class = "btn-sm",
        style = "width: 30px; margin-right: 5px;"
      ),
      actionButton(
        inputId = paste0("save_state_", i), 
        label = paste0("S", i),
        class = "btn-sm",
        style = "width: 30px; margin-right: 8px;"
      ),
      span(title_display, style = "font-size: 12px; color: #666;")
    )
  })
  
  div(
    h5("States", style = "margin-bottom: 10px;"),
    do.call(div, state_buttons)
  )
})

# create observers for load buttons
lapply(1:9, function(i) {
  observeEvent(input[[paste0("load_state_", i)]], {
    if (state_exists(i)) {
      success <- load_state(i)
      if (success) {
        showNotification(
          sprintf("loaded state %d", i),
          type = "message",
          duration = 2
        )
        # check if view changed and only refresh manually if it didn't
        current_view_id <- get_current_view_id()
        previous_view_id <- isolate(input$view_id)
        view_changed <- !is.null(current_view_id) && !identical(current_view_id, previous_view_id)
        
        if (!view_changed && exists("refresh_trigger")) {
          # only trigger manual refresh if view didn't change (set_view handles refresh when view changes)
          current_val <- refresh_trigger()
          refresh_trigger(current_val + 1)
          cat("plot refresh triggered after loading state (view unchanged)\n")
        }
        
        # update view dropdown to reflect loaded view (with flag to prevent recursion)
        if (!is.null(current_view_id)) {
          updating_view_programmatically(TRUE)
          updateSelectInput(session, "view_id", selected = current_view_id)
          # reset flag after a short delay
          invalidateLater(100)
          updating_view_programmatically(FALSE)
        }
        # trigger UI refresh to update button titles if needed
        state_buttons_trigger(state_buttons_trigger() + 1)
      } else {
        showNotification(
          sprintf("failed to load state %d", i),
          type = "error",
          duration = 3
        )
      }
    } else {
      showNotification(
        sprintf("state %d does not exist", i),
        type = "warning",
        duration = 2
      )
    }
  })
})

# create observers for save buttons
lapply(1:9, function(i) {
  observeEvent(input[[paste0("save_state_", i)]], {
    # show modal dialog to get title
    showModal(modalDialog(
      title = sprintf("Save State %d", i),
      textInput(
        inputId = "state_title_input",
        label = "Title:",
        value = get_state_title(i),  # pre-fill with existing title if any
        placeholder = "enter a descriptive title for this state"
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_save_state", "Save", class = "btn-primary")
      )
    ))
    
    # store which state number we're saving
    session$userData$saving_state_number <- i
  })
})

# reactive value to trigger state buttons UI refresh
state_buttons_trigger <- reactiveVal(0)

# flag to prevent recursive view updates
updating_view_programmatically <- reactiveVal(FALSE)

# handle save confirmation
observeEvent(input$confirm_save_state, {
  state_number <- session$userData$saving_state_number
  title <- input$state_title_input
  
  if (is.null(title) || nchar(trimws(title)) == 0) {
    showNotification(
      "title cannot be empty",
      type = "error",
      duration = 3
    )
    return()
  }
  
  success <- save_state(state_number, trimws(title))
  
  if (success) {
    showNotification(
      sprintf("saved state %d: '%s'", state_number, trimws(title)),
      type = "message",
      duration = 3
    )
    # trigger UI refresh to show updated title
    state_buttons_trigger(state_buttons_trigger() + 1)
  } else {
    showNotification(
      sprintf("failed to save state %d", state_number),
      type = "error",
      duration = 3
    )
  }
  
  removeModal()
  session$userData$saving_state_number <- NULL
})

# create keyboard shortcuts for states
create_states_keyboard_shortcuts <- function() {
  # create load shortcuts (ctrl+N)
  for (i in 1:9) {
    shortcut_key <- paste0("Ctrl+", i)
    keyboard_shortcuts[[shortcut_key]] <<- list(
      description = paste0("load state ", i),
      type = "state_load",
      state_number = i,
      action = local({
        state_num <- i
        function(state_number, session_obj) {
          if (state_exists(state_number)) {
            # get view before loading to compare
            view_before_load <- get_current_view_id()
            
            success <- load_state(state_number)
            if (success) {
              # check if view changed and only refresh manually if it didn't
              current_view_id <- get_current_view_id()
              view_changed <- !identical(current_view_id, view_before_load)
              
              if (!view_changed && exists("refresh_trigger")) {
                # only trigger manual refresh if view didn't change
                current_val <- refresh_trigger()
                refresh_trigger(current_val + 1)
                cat("plot refresh triggered after keyboard load (view unchanged)\n")
              }
              
              # update view dropdown to reflect loaded view (with flag to prevent recursion)
              if (!is.null(current_view_id)) {
                updating_view_programmatically(TRUE)
                updateSelectInput(session_obj, "view_id", selected = current_view_id)
                # reset flag after a short delay
                invalidateLater(100)
                updating_view_programmatically(FALSE)
              }
              cat(sprintf("keyboard: loaded state %d\n", state_number))
            } else {
              cat(sprintf("keyboard: failed to load state %d\n", state_number))
            }
          } else {
            cat(sprintf("keyboard: state %d does not exist\n", state_number))
          }
        }
      })
    )
  }
  
  # create save shortcuts (ctrl+alt+N)
  for (i in 1:9) {
    shortcut_key <- paste0("Ctrl+Alt+", i)
    keyboard_shortcuts[[shortcut_key]] <<- list(
      description = paste0("save state ", i),
      type = "state_save",
      state_number = i,
      action = local({
        state_num <- i
        function(state_number, session_obj) {
          # trigger save dialog by simulating button click
          session_obj$sendCustomMessage("triggerStateSave", state_number)
        }
      })
    )
  }
}

# call this function to register the shortcuts
create_states_keyboard_shortcuts()
