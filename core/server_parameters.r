# Server-side logic for parameter UI and observation

# Create UI tab for parameters
output$parameters_ui <- renderUI({
  param_ui_trigger() # Establish reactive dependency on the counter

  all_params <- get_registered_parameters()

  if (length(all_params) == 0) {
    return(shiny::tagList(shiny::p("No parameters registered.")))
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

  # Create UI elements for each group
  ui_elements <- lapply(names(grouped_params), function(group_name) {
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
        # Default case
        shiny::p(paste("Unknown parameter type:", param$type))
      )

      return(input_element)
    })

    # Create a section for this group
    shiny::div(
      shiny::h4(group_name),
      param_inputs,
      shiny::hr()
    )
  })

  # Return all UI elements
  return(shiny::tagList(ui_elements))
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
