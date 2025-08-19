# keyboard.r
# Handles keyboard shortcuts and events

shifted_numbers <- c("!", "@", "#", "$", "%", "^", "&", "*", "(", ")")

# helper function to zoom to 10^N basepairs around current center
zoom_to_power_of_ten <- function(n, states_module_output, main_state_rv) {
  # calculate target range (10^N basepairs)
  target_range <- 10^n
  half_range <- target_range / 2
  
  # determine current center
  current_center <- if (!is.null(main_state_rv$zoom)) {
    # if zoomed, use center of current zoom
    (main_state_rv$zoom[1] + main_state_rv$zoom[2]) / 2
  } else {
    # if not zoomed, need to get center of all selected contigs
    # this requires getting the contig data to find the center
    NULL
  }
  
  if (is.null(current_center)) {
    # fallback: try to get center from contig data if available
    if (length(main_state_rv$contigs) > 0) {
      # attempt to get center from first contig - this is a fallback
      cat(sprintf("alt+%d: no current zoom, cannot determine center\n", n))
      return()
    } else {
      cat(sprintf("alt+%d: no contigs selected, cannot determine center\n", n))
      return()
    }
  }
  
  # save current state before changing zoom
  states_module_output$push_state()
  
  # set new zoom range centered on current center
  new_start <- current_center - half_range
  new_end <- current_center + half_range
  main_state_rv$zoom <- c(new_start, new_end)
  
  cat(sprintf("zoomed to %s basepairs around center %.0f\n", 
              format(target_range, big.mark = ","), current_center))
}

# Define keyboard shortcuts globally
keyboard_shortcuts <- list(
  "Alt+backspace" = list(
    description = "Undo last action",
    type = "general", # Added type for dispatch
    action = function(states_module_output) {
      states_module_output$undo_state()
    }
  ),
  "Alt+0" = list(
    description = "Reset zoom to full range",
    type = "zoom", # Changed from general to zoom type for direct state access
    action = function(states_module_output, main_state_rv) {
      states_module_output$push_state()
      main_state_rv$zoom <- NULL
    }
  ),
  "Alt+Z" = list(
    description = "Update zoom from selected range",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      if (!is.null(main_state_rv$current_xlim)) {
        states_module_output$push_state()
        main_state_rv$zoom <- main_state_rv$current_xlim
      }
    }
  ),
  "Alt+_" = list(
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
      }
    }
  ),
  "Alt+=" = list(
    description = "Zoom in (halve the view range)",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        states_module_output$push_state()
        zoom_range <- main_state_rv$zoom[2] - main_state_rv$zoom[1]
        # halve the range while keeping the same center
        main_state_rv$zoom <- c(
          main_state_rv$zoom[1] + zoom_range * 0.25,
          main_state_rv$zoom[2] - zoom_range * 0.25
        )
      }
    }
  ),
  "Alt+>" = list(
    description = "Move zoom area one window to the right",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        states_module_output$push_state()
        zoom_range <- main_state_rv$zoom[2] - main_state_rv$zoom[1]
        main_state_rv$zoom <- c(
          main_state_rv$zoom[1] + zoom_range,
          main_state_rv$zoom[2] + zoom_range
        )
      }
    }
  ),
  "Alt+<" = list(
    description = "Move zoom area one window to the left",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        states_module_output$push_state()
        zoom_range <- main_state_rv$zoom[2] - main_state_rv$zoom[1]
        main_state_rv$zoom <- c(
          main_state_rv$zoom[1] - zoom_range,
          main_state_rv$zoom[2] - zoom_range
        )
        cat("moved zoom area left by Alt+<\n")
      } else {
        cat("Alt+<: zoom is null, cannot move zoom area\n")
      }
    }
  ),
  "Alt+C" = list(
    description = "Clear selected contigs",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      states_module_output$push_state()
      main_state_rv$contigs <- character()
      cat("cleared selected contigs by Alt+C\n")
    }
  ),
  "Alt+2" = list(
    description = "Zoom to 100 basepairs around current center",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      zoom_to_power_of_ten(2, states_module_output, main_state_rv)
    }
  ),
  "Alt+3" = list(
    description = "Zoom to 1,000 basepairs around current center",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      zoom_to_power_of_ten(3, states_module_output, main_state_rv)
    }
  ),
  "Alt+4" = list(
    description = "Zoom to 10,000 basepairs around current center",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      zoom_to_power_of_ten(4, states_module_output, main_state_rv)
    }
  ),
  "Alt+5" = list(
    description = "Zoom to 100,000 basepairs around current center",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      zoom_to_power_of_ten(5, states_module_output, main_state_rv)
    }
  ),
  "Alt+6" = list(
    description = "Zoom to 1,000,000 basepairs around current center",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      zoom_to_power_of_ten(6, states_module_output, main_state_rv)
    }
  ),
  "Alt+7" = list(
    description = "Zoom to 10,000,000 basepairs around current center",
    type = "zoom",
    action = function(states_module_output, main_state_rv) {
      zoom_to_power_of_ten(7, states_module_output, main_state_rv)
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
    shortcut_key <- paste0("Alt+", shifted_numbers[i])
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
      // prevent default behavior for Alt key combinations first
      if (e.altKey && !e.metaKey && !e.ctrlKey) {
        e.preventDefault();
        e.stopPropagation();
      }

      let keyCombo = [];
      if (e.metaKey) keyCombo.push('Command');
      if (e.ctrlKey) keyCombo.push('Ctrl');
      if (e.altKey) keyCombo.push('Alt');
      if (e.shiftKey) keyCombo.push('Shift');

      // use e.code for more reliable key detection
      let key = e.code;
      
      // convert common codes to readable format
      if (key.startsWith('Key')) {
        key = key.replace('Key', '');
      } else if (key.startsWith('Digit')) {
        key = key.replace('Digit', '');
      } else if (key === 'Backspace') {
        key = 'backspace';
      } else if (key === 'Comma') {
        key = '<';
      } else if (key === 'Period') {
        key = '>';
      } else if (key === 'Minus') {
        key = '_';
      } else if (key === 'Equal') {
        key = '=';
      }
      
      if (key === 'Meta' || key === 'Alt' || key === 'Control' || key === 'Shift') return;
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

  # About button shows version and license info
  observeEvent(input$aboutBtn, {
    # gather system statistics
    num_assemblies <- length(get_assemblies())
    num_views <- length(get_view_ids())
    num_profiles <- length(profiles_get_all())
    num_tabs <- length(get_registered_tabs())
    
    # build statistics text
    stats_text <- paste(
      sprintf("<p><strong>Registered objects:</strong></p>"),
      sprintf("<ul>"),
      sprintf("<li>Assemblies: %d</li>", num_assemblies),
      sprintf("<li>Views: %d</li>", num_views),
      sprintf("<li>Profiles: %d</li>", num_profiles),
      sprintf("<li>Dynamic tabs: %d</li>", num_tabs),
      sprintf("</ul>"),
      sep = ""
    )
    
    showModal(modalDialog(
      title = "About mview",
      HTML(paste(
        "<p><strong>mview v1.01</strong></p>",
        "<p>Developed by Eitan Yaffe</p>",
        "<p>August 2025</p>",
        "<p>This software is provided under an open source license.</p>",
        "<hr>",
        stats_text,
        sep = ""
      )),
      easyClose = TRUE
    ))
  })

  # Output for the last key pressed
  output$last_key_output <- renderText({
    paste(input$key_pressed, " ")
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
  # separate alt+N zoom shortcuts for cleaner display
  alt_zoom_shortcuts <- grep("^Alt\\+[2-7]$", names(keyboard_shortcuts), value = TRUE)
  other_shortcuts <- setdiff(names(keyboard_shortcuts), alt_zoom_shortcuts)
  
  html_content <- "<h5>Actions</h5><ul>"
  
  # add non-tab switch, non-alt+N shortcuts first
  for (shortcut_name in other_shortcuts) {
    type <- keyboard_shortcuts[[shortcut_name]]$type
    if (type == "tab_switch") {
      next
    }
    html_content <- paste0(
      html_content,
      "<li><strong>", shortcut_name, ":</strong> ",
      keyboard_shortcuts[[shortcut_name]]$description, "</li>"
    )
  }
  
  # add alt+N shortcuts as a group
  if (length(alt_zoom_shortcuts) > 0) {
    html_content <- paste0(html_content, 
      "<li><strong>Alt+N (N=2-7):</strong> zoom to 10^N basepairs around current center</li>")
  }
  
  html_content <- paste0(html_content, "</ul>")

  return(HTML(html_content))
}
