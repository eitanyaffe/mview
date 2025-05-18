# keyboard.r
# Handles keyboard shortcuts and events

# Define keyboard shortcuts globally
keyboard_shortcuts <- list(
  "Ctrl+z" = list(
    description = "Undo last action",
    action = function(states_module_output) {
      states_module_output$undo_state()
    }
  ),
  "Ctrl+s" = list(
    description = "Zoom to current region",
    action = function(states_module_output) {
      # Code to zoom to current region
      # This will be implemented later
      print("Zoom to current region action triggered")
    }
  )
)

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
keyboard_server <- function(input, session, states_module_output) {
  # Create a reactive value to track keyboard events
  keyboard_events <- reactiveVal(0)

  # Handle keyboard events
  observeEvent(input$key_pressed, {
    key_combo <- input$key_pressed

    # Check if the key combo matches any of our shortcuts
    if (!is.null(key_combo) && key_combo != "") {
      for (shortcut in names(keyboard_shortcuts)) {
        if (tolower(key_combo) == tolower(shortcut)) {
          keyboard_shortcuts[[shortcut]]$action(states_module_output)
          keyboard_events(keyboard_events() + 1)
          break
        }
      }
    }
  })

  # Help button shows keyboard shortcuts
  observeEvent(input$helpBtn, {
    showModal(modalDialog(
      title = "Help",
      keyboard_summary(),
      easyClose = TRUE
    ))
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
  html_content <- "<h4>Keyboard Shortcuts</h4><ul>"
  for (shortcut in names(keyboard_shortcuts)) {
    html_content <- paste0(
      html_content,
      "<li><strong>", shortcut, ":</strong> ",
      keyboard_shortcuts[[shortcut]]$description, "</li>"
    )
  }
  html_content <- paste0(html_content, "</ul>")

  return(HTML(html_content))
}
