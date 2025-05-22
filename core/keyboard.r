# keyboard.r
# Handles keyboard shortcuts and events

shifted_numbers <- c("!", "@", "#", "$", "%", "^", "&", "*", "(", ")")

# Define keyboard shortcuts globally
keyboard_shortcuts <- list(
  "Shift+backspace" = list(
    description = "Undo last action",
    type = "general", # Added type for dispatch
    action = function(states_module_output) {
      states_module_output$undo_state()
    }
  ),
  "Shift+X" = list(
    description = "Reset zoom to full range",
    type = "zoom", # Changed from general to zoom type for direct state access
    action = function(states_module_output, main_state_rv) {
      states_module_output$push_state()
      main_state_rv$zoom <- NULL
      cat("zoom reset to NULL by Shift+X\n")
    }
  ),
  "Shift+Z" = list(
    description = "Update zoom from selected range",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      if (!is.null(main_state_rv$current_xlim)) {
        states_module_output$push_state()
        main_state_rv$zoom <- main_state_rv$current_xlim
        cat("Zoom updated by Shift+Z (main_state_rv accessed)\n")
      } else {
        cat("Shift+Z: main_state_rv or main_state_rv$current_xlim is null\n")
      }
    }
  ),
  "Shift+_" = list(
    description = "Zoom out (double the view range)",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        states_module_output$push_state()
        zoom_range <- main_state_rv$zoom[2] - main_state_rv$zoom[1]
        main_state_rv$zoom <- c(
          main_state_rv$zoom[1] - zoom_range * 0.5,
          main_state_rv$zoom[2] + zoom_range * 0.5
        )
        cat("zoomed out by Shift+-\n")
      } else {
        cat("Shift+-: zoom is null, cannot zoom out\n")
      }
    }
  )
)

# Function to generate and add tab-switching shortcuts
create_and_add_tab_shortcuts <- function() {
  # Fixed vector of tab titles as they appear in the UI
  tabs <- c(
    "States", "Contigs", "Genomes", "Contig Map", "Selected Contigs",
    "Options", "Parameters"
  )

  for (i in seq_along(tabs)) {
    shortcut_key <- paste0("Shift+", shifted_numbers[i])
    tab <- tabs[i]

    keyboard_shortcuts[[shortcut_key]] <<- list(
      description = paste0("Switch to ", tab, " tab"),
      type = "tab_switch",
      action = local({
        tab_value <- tab
        function(session_obj) {
          updateTabsetPanel(session_obj, "mainTabs", selected = tab_value)
        }
      })
    )
  }
}

# Call this function once to populate the shortcuts when the script is sourced
create_and_add_tab_shortcuts()

# Initialize keyboard shortcuts
keyboard_initialize <- function() {
  keyboard_js <- HTML("
    document.addEventListener('keydown', function(e) {
      let keyCombo = [];
      if (e.metaKey) keyCombo.push('Command');
      if (e.ctrlKey) keyCombo.push('Ctrl');
      if (e.altKey) keyCombo.push('Alt');
      if (e.shiftKey) keyCombo.push('Shift');

      let key = e.key;
      if (key === 'Meta') return; // skip raw 'Meta' keypress
      keyCombo.push(key);

      Shiny.setInputValue('key_pressed', keyCombo.join('+'), {priority: 'event'});
    });
  ")

  return(tags$script(keyboard_js))
}

# Server-side handler for keyboard shortcuts
keyboard_server <- function(input, output, session, main_state_rv, states_module_output) {
  # Create a reactive value to track keyboard events
  keyboard_events <- reactiveVal(0)

  # Handle keyboard events
  observeEvent(input$key_pressed, {
    key_combo <- input$key_pressed

    # Check if the key combo matches any of our shortcuts
    if (!is.null(key_combo) && key_combo != "") {
      current_shortcut_names <- names(keyboard_shortcuts) # Use a local copy of names
      for (shortcut_name in current_shortcut_names) {
        if (tolower(key_combo) == tolower(shortcut_name)) {
          shortcut_details <- keyboard_shortcuts[[shortcut_name]]

          # Dispatch action based on type
          if (!is.null(shortcut_details$type) && shortcut_details$type == "tab_switch") {
            shortcut_details$action(session) # Pass session for tab switching
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "general") {
            # Default for "general" type or if type is not specified
            shortcut_details$action(states_module_output)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "zoom") {
            shortcut_details$action(states_module_output, main_state_rv)
          } else {
            cat("Unknown shortcut type:", shortcut_details$type, "\n")
          }

          keyboard_events(keyboard_events() + 1)
          break
        }
      }
    }
  })

  # Help button shows keyboard shortcuts
  observeEvent(input$helpBtn, {
    showModal(modalDialog(
      title = "Keyboard Shortcuts",
      keyboard_summary(),
      easyClose = TRUE
    ))
  })

  # Output for the last key pressed
  output$last_key_output <- renderText({
    paste("Pressed:", input$key_pressed)
  })

  return(list(
    keyboard_events = keyboard_events,
    last_key_pressed = reactive({
      input$key_pressed
    })
  ))
}

# Generate keyboard shortcuts documentation for help
keyboard_summary <- function() {
  html_content <- "<h5>Actions</h5><ul>"
  for (i in seq_along(keyboard_shortcuts)) {
    shortcut <- names(keyboard_shortcuts)[i]
    type <- keyboard_shortcuts[[i]]$type
    if (type == "tab_switch") {
      next
    }
    html_content <- paste0(
      html_content,
      "<li><strong>", shortcut, ":</strong> ",
      keyboard_shortcuts[[shortcut]]$description, "</li>"
    )
  }
  html_content <- paste0(html_content, "</ul>")

  return(HTML(html_content))
}
