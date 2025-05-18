# keyboard.r
# Handles keyboard shortcuts and events

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

  # Define keyboard shortcuts
  keyboard_shortcuts <- list(
    "Shift+Backspace" = list(
      description = "Undo last action",
      action = function() {
        states_module_output$undo_state()
      }
    ),
    "Z" = list(
      description = "Zoom to current region",
      action = function() {
        # Code to zoom to current region
        # This will be implemented later
        print("Zoom to current region action triggered")
      }
    )
  )

  # Handle keyboard events
  observeEvent(input$key_pressed, {
    key_combo <- input$key_pressed

    # Log the key press for debugging
    # print(paste("Key pressed:", key_combo)) # Already done in server logic

    # Check if the key combo matches any of our shortcuts
    if (!is.null(key_combo) && key_combo != "") {
      for (shortcut in names(keyboard_shortcuts)) {
        if (tolower(key_combo) == tolower(shortcut)) {
          keyboard_shortcuts[[shortcut]]$action()
          keyboard_events(keyboard_events() + 1)
          break
        }
      }
    }
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
  shortcuts <- list(
    "Shift+Backspace" = "Undo last action",
    "Z" = "Zoom to current region"
  )

  html_content <- "<h4>Keyboard Shortcuts</h4><ul>"
  for (shortcut in names(shortcuts)) {
    html_content <- paste0(
      html_content,
      "<li><strong>", shortcut, ":</strong> ",
      shortcuts[[shortcut]], "</li>"
    )
  }
  html_content <- paste0(html_content, "</ul>")

  return(HTML(html_content))
}
