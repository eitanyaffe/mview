# Server-side logic for parameter UI and observation

# Create UI for parameter tabs in right panel
output$parameter_tabs_ui <- renderUI({
  param_ui_trigger() # Establish reactive dependency on the counter

  all_params <- get_registered_parameters()

  if (length(all_params) == 0) {
    return(NULL) # Display nothing if no parameters
  }

  # Group parameters for UI organization
  grouped_params <- list()
  for (param_key in names(all_params)) {
    param <- all_params[[param_key]]
    group <- param$group

    if (!is.element(group, names(grouped_params))) {
      grouped_params[[group]] <- list()
    }
    grouped_params[[group]][[param$id]] <- param
  }

  # Create tab panels for each group
  tab_panels <- lapply(names(grouped_params), function(group_name) {
    group_params <- grouped_params[[group_name]]

    # Create inputs for each parameter in this group
    param_inputs <- lapply(names(group_params), function(param_id) {
      param <- group_params[[param_id]]
      input_id <- paste0("param_input_", param$group, "_", param$id)
      current_val <- param$reactive_accessor()

      # Create appropriate Shiny input based on param type
      input_element <- switch(param$type,
        "string" = textInput(
          inputId = input_id,
          label = param$id,
          value = current_val
        ),
        "integer" = numericInput(
          inputId = input_id,
          label = param$id,
          value = current_val,
          step = 1
        ),
        "double" = numericInput(
          inputId = input_id,
          label = param$id,
          value = current_val
        ),
        "boolean" = checkboxInput(
          inputId = input_id,
          label = param$id,
          value = current_val
        ),
        "select" = radioButtons(
          inputId = input_id,
          label = param$id,
          choices = param$choices,
          selected = current_val
        ),
        # Default case
        shiny::p(paste("Unknown parameter type:", param$type))
      )

      return(input_element)
    })

    # Create a tab panel for this group
    shiny::tabPanel(
      title = group_name,
      shiny::div(
        class = "parameter-group-content",
        param_inputs
      )
    )
  })

  # Create tabset panel with all parameter groups
  do.call(shiny::tabsetPanel, c(list(id = "parameterTabs"), tab_panels))
})



# Observe parameter inputs and update values
observe({
  param_ui_trigger() # Establish reactive dependency on the counter

  # Get all registered parameters
  all_params <- get_registered_parameters()
  # For each parameter, create an observer for its UI input
  lapply(names(all_params), function(param_key) {
    param <- all_params[[param_key]]
    input_id <- paste0("param_input_", param$group, "_", param$id)

    # Create an observer for this specific input
    observeEvent(input[[input_id]],
      {
        if (!is.null(input[[input_id]])) {
          update_parameter_from_ui(param$group, param$id, input[[input_id]])
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )
  })
})

# Server logic for collapsible parameter panel
parameter_panel_collapsed <- reactiveVal(FALSE)

observeEvent(input$toggleParameterPanel, {
  current_state <- parameter_panel_collapsed()
  parameter_panel_collapsed(!current_state)
  
  # toggle CSS classes and update button icon
  if (!current_state) {
    # collapsing
    shinyjs::addClass(id = "parameter-panel", class = "collapsed")
    shinyjs::addClass(id = "parameter-panel-column", class = "collapsed")
    shinyjs::addClass(id = "profile-plots-column", class = "expanded")
    shinyjs::html(id = "toggleParameterPanel", 
                  html = '<i class="fa fa-chevron-right"></i>')
  } else {
    # expanding
    shinyjs::removeClass(id = "parameter-panel", class = "collapsed")
    shinyjs::removeClass(id = "parameter-panel-column", class = "collapsed")
    shinyjs::removeClass(id = "profile-plots-column", class = "expanded")
    shinyjs::html(id = "toggleParameterPanel", 
                  html = '<i class="fa fa-chevron-left"></i>')
  }
})
