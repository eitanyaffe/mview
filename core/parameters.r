# Parameter management system
# Handles registration and access to application parameters

# 1. Initialize Parameter Storage
# param_registry: A named list to store parameter definitions and their reactive accessors.
# Each element will be a list: list(group, id, type, default_value, reactive_accessor)
param_registry <- list()

# Trigger for UI and observer updates
param_ui_trigger <- shiny::reactiveVal(0)

# 2. register_param(group, id, type, default_value)
#    - group (character): The group name for UI organization.
#    - id (character): The unique ID for the parameter within its group.
#    - type (character): "string", "integer", "double", "boolean".
#    - default_value: The default value, matching the specified type.
register_param <- function(group, id, type, default_value) {
  cat(sprintf("register_param: %s_%s\n", group, id))

  param_key_registry <- paste0(group, "_", id) # Key for the param_registry
  param_key_cache <- paste0("param_", group, "_", id) # Key for the cache module

  if (is.element(param_key_registry, names(param_registry))) {
    stop(sprintf("parameter %s in group %s already registered.", id, group))
  }

  # Validate type and default_value
  if (!type %in% c("string", "integer", "double", "boolean")) {
    stop(sprintf("invalid parameter type: %s for %s_%s", type, group, id))
  }

  # Create the reactive value for this parameter
  current_rv <- shiny::reactiveVal()

  # Initialize the reactive value: Cache > Default
  initial_value <- default_value
  if (cache_exists(param_key_cache)) {
    cached_val <- cache_get(param_key_cache)
    initial_value <- cached_val
  }
  current_rv(initial_value)

  # Store in registry
  param_registry[[param_key_registry]] <<- list(
    group = group,
    id = id,
    type = type,
    default_value = default_value,
    reactive_accessor = current_rv,
    cache_key = param_key_cache
  )

  return(current_rv)
}

# 3. get_param(group, id)
#    - Returns the reactive_accessor function for the parameter.
get_param <- function(group, id) {
  param_key_registry <- paste0(group, "_", id)
  if (!is.element(param_key_registry, names(param_registry))) {
    stop(sprintf("parameter %s in group %s not found.", id, group))
  }
  return(param_registry[[param_key_registry]]$reactive_accessor)
}

# 4. finalize_parameter_registration_and_trigger_ui()
#    - Call this function once after all parameters have been registered.
#    - It signals that the UI for parameters can now be built/updated.
param_registration_done <- function() {
  param_ui_trigger(param_ui_trigger() + 1)
}

# 5. Helper to get all registered parameters (for UI generation)
#    - This function will be called by the server-side UI rendering logic.
get_registered_parameters <- function() {
  # Return all registered parameters
  return(param_registry)
}

# 6. update_parameter_from_ui(group, id, new_value)
#    - This function will be called by observeEvent handlers in server_parameters.r
#      when a UI input changes.
#    - It updates the reactiveVal and the cache.
update_parameter_from_ui <- function(group, id, new_value) {
  param_key_registry <- paste0(group, "_", id)
  if (!is.element(param_key_registry, names(param_registry))) {
    warning(sprintf("attempted to update non-existent parameter %s in group %s from UI.", id, group))
    return(invisible(NULL))
  }

  param_definition <- param_registry[[param_key_registry]]

  # Type coercion based on param_definition$type
  coerced_value <- new_value
  tryCatch(
    {
      if (param_definition$type == "integer") {
        if (!is.integer(new_value)) {
          temp_val <- as.integer(new_value)
          if (is.na(temp_val) || temp_val != new_value) {
            stop("value cannot be safely coerced to integer.")
          }
          coerced_value <- temp_val
        }
      } else if (param_definition$type == "double") {
        coerced_value <- as.numeric(new_value)
      } else if (param_definition$type == "string") {
        coerced_value <- as.character(new_value)
      } else if (param_definition$type == "boolean") {
        coerced_value <- as.logical(new_value)
      }

      if (is.na(coerced_value) && !is.na(new_value) &&
        !(param_definition$type == "string" && new_value == "NA")) {
        stop("type coercion failed, resulted in NA.")
      }
    },
    error = function(e) {
      warning(sprintf(
        "type coercion failed for parameter %s_%s. Value '%s' to type '%s'. Error: %s.",
        group, id, new_value, param_definition$type, e$message
      ))
      return(invisible(NULL))
    }
  )

  # Update reactive value
  param_definition$reactive_accessor(coerced_value)

  # Update cache
  cache_set(param_definition$cache_key, coerced_value)
}

# Clear all registered parameters
clear_parameters <- function(clear_cache = FALSE) {
  # If clear_cache is TRUE, remove all parameter values from cache
  if (clear_cache) {
    for (param_key in names(param_registry)) {
      param <- param_registry[[param_key]]
      if (cache_exists(param$cache_key)) {
        cache_unset(param$cache_key)
      }
    }
  }

  # Clear the registry
  param_registry <<- list()

  # Trigger UI update by changing the counter
  param_ui_trigger(param_ui_trigger() + 1)
}
