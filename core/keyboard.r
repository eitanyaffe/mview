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
# Helper function for contig navigation
change_contig_action <- function(regions_module_output, main_state_rv, input, direction) {
  # check if DataTable is initialized
  if (is.null(input$contigTable_rows_all)) {
    cat(sprintf("alt+arrow%s: contig table not initialized\n", direction))
    return()
  }
  
  # get current table data (filtered/sorted as displayed)
  contigs_data <- get_contigs(main_state_rv$assembly)
  if (is.null(contigs_data) || nrow(contigs_data) == 0) {
    cat(sprintf("alt+arrow%s: no contigs available\n", direction))
    return()
  }
  
  # apply same filtering as the table
  if (!is.null(contigs_data) && "length" %in% names(contigs_data)) {
    min_length_kb <- if (is.null(input$min_contig_length)) 0 else input$min_contig_length
    max_length_kb <- input$max_contig_length
    min_length_bp <- min_length_kb * 1000
    max_length_bp <- if (!is.null(max_length_kb) && !is.na(max_length_kb)) max_length_kb * 1000 else NULL
    contigs_data <- contigs_data[contigs_data$length >= min_length_bp, ]
    if (!is.null(max_length_bp)) {
      contigs_data <- contigs_data[contigs_data$length <= max_length_bp, ]
    }
  }
  
  # get current row selection
  selected_rows <- input$contigTable_rows_selected
  
  if (is.null(selected_rows) || length(selected_rows) == 0) {
    # no selection, select first row
    regions_module_output$push_undo_state()
    contig_proxy <- DT::dataTableProxy("contigTable")
    DT::selectRows(contig_proxy, 1)
    main_state_rv$contigs <- contigs_data$contig[input$contigTable_rows_all[1]]
    main_state_rv$zoom <- NULL
    cat("selected first contig in table\n")
    return()
  }
  
  if (length(selected_rows) != 1) {
    # multiple rows selected, do nothing
    cat(sprintf("alt+arrow%s: multiple rows selected, doing nothing\n", direction))
    return()
  }
  
  # find current position in sorted table
  current_table_index <- which(input$contigTable_rows_all == selected_rows[1])
  
  # find next/previous position based on direction
  if (direction == "down") {
    # find next position (wrap to first if at end)
    if (current_table_index >= length(input$contigTable_rows_all)) {
      new_table_index <- 1
    } else {
      new_table_index <- current_table_index + 1
    }
  } else {
    # find previous position (wrap to last if at beginning)
    if (current_table_index <= 1) {
      new_table_index <- length(input$contigTable_rows_all)
    } else {
      new_table_index <- current_table_index - 1
    }
  }
  
  # get the actual row index for the new position
  new_row_index <- input$contigTable_rows_all[new_table_index]
  new_contig <- contigs_data$contig[new_row_index]
  
  # push undo state before changing
  regions_module_output$push_undo_state()
  
  # update selection and contig
  contig_proxy <- DT::dataTableProxy("contigTable")
  DT::selectRows(contig_proxy, new_row_index)
  main_state_rv$contigs <- new_contig
  main_state_rv$zoom <- NULL
  
  cat(sprintf("navigated to %s contig: %s\n", 
              if (direction == "down") "next" else "previous", 
              new_contig))
}

# Helper function for region navigation
change_region_action <- function(regions_module_output, main_state_rv, direction) {
  # get current regions table
  region_table_data <- regions_module_output$get_region_table()
  if (is.null(region_table_data) || nrow(region_table_data) == 0) {
    cat(sprintf("ctrl+alt+arrow%s: no regions available\n", direction))
    return()
  }
  
  # determine current region by matching current state
  current_assembly <- main_state_rv$assembly
  current_contigs <- main_state_rv$contigs
  current_zoom <- main_state_rv$zoom
  
  # find matching region or get first one if none match
  current_index <- -1
  
  # try to find current region by matching state
  for (i in seq_len(nrow(region_table_data))) {
    region_row <- region_table_data[i, ]
    
    # check if this region matches current state
    region_contigs <- if (region_row$contigs == "" || is.na(region_row$contigs)) {
      character(0)
    } else {
      trimws(strsplit(region_row$contigs, ",")[[1]])
    }
    
    region_zoom <- if (is.na(region_row$zoom_start) || is.na(region_row$zoom_end)) {
      NULL
    } else {
      c(region_row$zoom_start, region_row$zoom_end)
    }
    
    # check for match
    assembly_match <- identical(current_assembly, region_row$assembly)
    contigs_match <- identical(sort(current_contigs), sort(region_contigs))
    zoom_match <- identical(current_zoom, region_zoom)
    
    if (assembly_match && contigs_match && zoom_match) {
      current_index <- i
      break
    }
  }
  
  # determine next/previous region index based on direction
  if (direction == "down") {
    # find next region
    if (current_index == -1) {
      # no current region found, go to first region
      new_index <- 1
    } else if (current_index >= nrow(region_table_data)) {
      # at last region, wrap to first
      new_index <- 1
    } else {
      new_index <- current_index + 1
    }
  } else {
    # find previous region
    if (current_index == -1) {
      # no current region found, go to last region
      new_index <- nrow(region_table_data)
    } else if (current_index <= 1) {
      # at first region, wrap to last
      new_index <- nrow(region_table_data)
    } else {
      new_index <- current_index - 1
    }
  }
  
  # navigate to new region
  region_row <- region_table_data[new_index, ]
  regions_module_output$goto_region(region_row)
  
  # update the regions table selection
  regions_proxy <- DT::dataTableProxy("regions_module-region_table_display")
  DT::selectRows(regions_proxy, new_index)
  
  cat(sprintf("navigated to %s region: %s (index %d)\n", 
              if (direction == "down") "next" else "previous", 
              region_row$description, new_index))
}

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
  "Alt+ArrowRight" = list(
    description = "Move zoom area 3/4 window to the right",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        regions_module_output$push_undo_state()
        zoom_range <- main_state_rv$zoom[2] - main_state_rv$zoom[1]
        move_distance <- zoom_range * 0.75
        main_state_rv$zoom <- c(
          main_state_rv$zoom[1] + move_distance,
          main_state_rv$zoom[2] + move_distance
        )
      }
    }
  ),
  "Alt+ArrowLeft" = list(
    description = "Move zoom area 3/4 window to the left",
    type = "zoom",
    action = function(regions_module_output, main_state_rv) {
      if (!is.null(main_state_rv$zoom)) {
        regions_module_output$push_undo_state()
        zoom_range <- main_state_rv$zoom[2] - main_state_rv$zoom[1]
        move_distance <- zoom_range * 0.75
        main_state_rv$zoom <- c(
          main_state_rv$zoom[1] - move_distance,
          main_state_rv$zoom[2] - move_distance
        )
        cat("moved zoom area left by Alt+ArrowLeft\n")
      } else {
        cat("Alt+ArrowLeft: zoom is null, cannot move zoom area\n")
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
    action = function(regions_module_output, main_state_rv, input) {
      change_contig_action(regions_module_output, main_state_rv, input, "down")
    }
  ),
  "Alt+ArrowUp" = list(
    description = "Navigate to previous contig in contig table",
    type = "navigation",
    action = function(regions_module_output, main_state_rv, input) {
      change_contig_action(regions_module_output, main_state_rv, input, "up")
    }
  ),
  "Ctrl+Alt+ArrowUp" = list(
    description = "Navigate to previous region",
    type = "region_navigation",
    action = function(regions_module_output, main_state_rv) {
      change_region_action(regions_module_output, main_state_rv, "up")
    }
  ),
  "Ctrl+Alt+ArrowDown" = list(
    description = "Navigate to next region",
    type = "region_navigation",
    action = function(regions_module_output, main_state_rv) {
      change_region_action(regions_module_output, main_state_rv, "down")
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
      } else if (key === 'ArrowLeft') {
        key = 'ArrowLeft';
      } else if (key === 'ArrowRight') {
        key = 'ArrowRight';
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
            shortcut_details$action(regions_module_output, main_state_rv, input)
          } else if (!is.null(shortcut_details$type) && shortcut_details$type == "region_navigation") {
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
