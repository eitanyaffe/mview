# keyboard.r
# Handles keyboard shortcuts and events

shifted_numbers <- c("!", "@", "#", "$", "%", "^", "&", "*", "(", ")")

# helper function to zoom to 10^N basepairs around current center
zoom_to_power_of_ten <- function(n, regions_module_output, main_state_rv) {
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
  
  # save current region before changing zoom
  regions_module_output$push_undo_state()
  
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
    action = function(regions_module_output) {
      regions_module_output$undo_last_action()
    }
  ),
  "Alt+0" = list(
    description = "Reset zoom to full range",
    type = "zoom", # Changed from general to zoom type for direct region access
    action = function(regions_module_output, main_state_rv) {
      regions_module_output$push_undo_state()
      main_state_rv$zoom <- NULL
    }
  ),
  "Alt+Z" = list(
    description = "Update zoom from selected range",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$current_xlim)) {
        regions_module_output$push_undo_state()
        main_state_rv$zoom <- main_state_rv$current_xlim
      }
    }
  ),
  "Alt+_" = list(
    description = "Zoom out (double the view range)",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        regions_module_output$push_undo_state()
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
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        regions_module_output$push_undo_state()
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
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        regions_module_output$push_undo_state()
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
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        regions_module_output$push_undo_state()
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
action = function(regions_module_output, main_state_rv) {
      regions_module_output$push_undo_state()
      main_state_rv$contigs <- character()
      cat("cleared selected contigs by Alt+C\n")
    }
  ),
  "Alt+2" = list(
    description = "Zoom to 100 basepairs around current center",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      zoom_to_power_of_ten(2, regions_module_output, main_state_rv)
    }
  ),
  "Alt+3" = list(
    description = "Zoom to 1,000 basepairs around current center",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      zoom_to_power_of_ten(3, regions_module_output, main_state_rv)
    }
  ),
  "Alt+4" = list(
    description = "Zoom to 10,000 basepairs around current center",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      zoom_to_power_of_ten(4, regions_module_output, main_state_rv)
    }
  ),
  "Alt+5" = list(
    description = "Zoom to 100,000 basepairs around current center",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      zoom_to_power_of_ten(5, regions_module_output, main_state_rv)
    }
  ),
  "Alt+6" = list(
    description = "Zoom to 1,000,000 basepairs around current center",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      zoom_to_power_of_ten(6, regions_module_output, main_state_rv)
    }
  ),
  "Alt+7" = list(
    description = "Zoom to 10,000,000 basepairs around current center",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      zoom_to_power_of_ten(7, regions_module_output, main_state_rv)
    }
  ),
  "Alt+S" = list(
    description = "Add current region",
    type = "region",
    action = function(regions_module_output) {
      # trigger the add region input focus
      regions_module_output$focus_add_input()
    }
  ),
  "Alt+ArrowDown" = list(
    description = "Navigate to next contig in contig table",
    type = "navigation",
    action = function(regions_module_output, main_state_rv) {
      # get contig table for current assembly
      contigs_data <- get_contigs(main_state_rv$assembly)
      if (is.null(contigs_data) || nrow(contigs_data) == 0) {
        cat("alt+arrowdown: no contigs available\n")
        return()
      }
      
      # check if any contigs are currently selected
      if (length(main_state_rv$contigs) == 0) {
        # no contigs selected, select the first contig in the table
        regions_module_output$push_undo_state()
        main_state_rv$contigs <- contigs_data$contig[1]
        main_state_rv$zoom <- NULL
        cat(sprintf("selected first contig: %s\n", contigs_data$contig[1]))
        return()
      }
      
      # find the position of the last currently selected contig in the table
      last_contig <- tail(main_state_rv$contigs, 1)
      current_index <- match(last_contig, contigs_data$contig)
      
      if (is.na(current_index)) {
        # current contig not found in table, select first contig
        regions_module_output$push_undo_state()
        main_state_rv$contigs <- contigs_data$contig[1]
        main_state_rv$zoom <- NULL
        cat(sprintf("current contig not found, selected first: %s\n", contigs_data$contig[1]))
        return()
      }
      
      # find next contig
      if (current_index >= nrow(contigs_data)) {
        # already at last contig, wrap to first
        next_index <- 1
      } else {
        next_index <- current_index + 1
      }
      
      next_contig <- contigs_data$contig[next_index]
      
      # push undo state before changing
      regions_module_output$push_undo_state()
      
      # set to next contig
      main_state_rv$contigs <- next_contig
      main_state_rv$zoom <- NULL
      
      cat(sprintf("navigated to next contig: %s (index %d)\n", next_contig, next_index))
    }
  ),
  "Alt+ArrowUp" = list(
    description = "Navigate to previous contig in contig table",
    type = "navigation",
    action = function(regions_module_output, main_state_rv) {
      # get contig table for current assembly
      contigs_data <- get_contigs(main_state_rv$assembly)
      if (is.null(contigs_data) || nrow(contigs_data) == 0) {
        cat("alt+arrowup: no contigs available\n")
        return()
      }
      
      # check if any contigs are currently selected
      if (length(main_state_rv$contigs) == 0) {
        # no contigs selected, select the last contig in the table
        regions_module_output$push_undo_state()
        main_state_rv$contigs <- contigs_data$contig[nrow(contigs_data)]
        main_state_rv$zoom <- NULL
        cat(sprintf("selected last contig: %s\n", contigs_data$contig[nrow(contigs_data)]))
        return()
      }
      
      # find the position of the last currently selected contig in the table
      last_contig <- tail(main_state_rv$contigs, 1)
      current_index <- match(last_contig, contigs_data$contig)
      
      if (is.na(current_index)) {
        # current contig not found in table, select last contig
        regions_module_output$push_undo_state()
        main_state_rv$contigs <- contigs_data$contig[nrow(contigs_data)]
        main_state_rv$zoom <- NULL
        cat(sprintf("current contig not found, selected last: %s\n", contigs_data$contig[nrow(contigs_data)]))
        return()
      }
      
      # find previous contig
      if (current_index <= 1) {
        # already at first contig, wrap to last
        prev_index <- nrow(contigs_data)
      } else {
        prev_index <- current_index - 1
      }
      
      prev_contig <- contigs_data$contig[prev_index]
      
      # push undo state before changing
      regions_module_output$push_undo_state()
      
      # set to previous contig
      main_state_rv$contigs <- prev_contig
      main_state_rv$zoom <- NULL
      
      cat(sprintf("navigated to previous contig: %s (index %d)\n", prev_contig, prev_index))
    }
  )
)

# Function to generate and add tab-switching shortcuts
create_and_add_tab_shortcuts <- function() {
  # Fixed vector of tab titles as they appear in the UI
  tabs <- c(
    "Regions", "Contigs", "Genomes", "Contig Map", "Selected Contigs",
    "Legends", "Options", "Parameters"
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
      } else if (key === 'ArrowDown') {
        key = 'ArrowDown';
      } else if (key === 'ArrowUp') {
        key = 'ArrowUp';
      }
      
      if (key === 'Meta' || key === 'Alt' || key === 'Control' || key === 'Shift') return;
      keyCombo.push(key);

      Shiny.setInputValue('key_pressed', keyCombo.join('+'), {priority: 'event'});
    });
  ")

  return(tags$script(keyboard_js))
}

# Server-side handler for keyboard shortcuts
keyboard_server <- function(input, output, session, main_state_rv, regions_module_output) {
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
            shortcut_details$action(regions_module_output)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "zoom") {
            shortcut_details$action(regions_module_output, main_state_rv)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "region") {
            shortcut_details$action(regions_module_output)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "navigation") {
            shortcut_details$action(regions_module_output, main_state_rv)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "state_load") {
            shortcut_details$action(shortcut_details$state_number, session)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "state_save") {
            shortcut_details$action(shortcut_details$state_number, session)
          } else {
            cat("unknown shortcut type:", shortcut_details$type, "\n")
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
  # separate different types of shortcuts for cleaner display
  alt_zoom_shortcuts <- grep("^Alt\\+[2-7]$", names(keyboard_shortcuts), value = TRUE)
  state_load_shortcuts <- grep("^Ctrl\\+[1-9]$", names(keyboard_shortcuts), value = TRUE)
  state_save_shortcuts <- grep("^Ctrl\\+Alt\\+[1-9]$", names(keyboard_shortcuts), value = TRUE)
  tab_shortcuts <- names(keyboard_shortcuts)[sapply(keyboard_shortcuts, function(x) x$type == "tab_switch")]
  other_shortcuts <- setdiff(names(keyboard_shortcuts), c(alt_zoom_shortcuts, state_load_shortcuts, state_save_shortcuts, tab_shortcuts))
  
  html_content <- "<h5>Actions</h5><ul>"
  
  # add main action shortcuts first
  for (shortcut_name in other_shortcuts) {
    html_content <- paste0(
      html_content,
      "<li><strong>", shortcut_name, ":</strong> ",
      keyboard_shortcuts[[shortcut_name]]$description, "</li>"
    )
  }
  
  # add grouped shortcuts
  if (length(alt_zoom_shortcuts) > 0) {
    html_content <- paste0(html_content, 
      "<li><strong>Alt+N (N=2-7):</strong> zoom to 10^N basepairs around current center</li>")
  }
  
  if (length(state_load_shortcuts) > 0) {
    html_content <- paste0(html_content, 
      "<li><strong>Ctrl+N (N=1-9):</strong> load state N</li>")
  }
  
  if (length(state_save_shortcuts) > 0) {
    html_content <- paste0(html_content, 
      "<li><strong>Ctrl+Alt+N (N=1-9):</strong> save state N</li>")
  }
  
  html_content <- paste0(html_content, "</ul>")

  return(HTML(html_content))
}
