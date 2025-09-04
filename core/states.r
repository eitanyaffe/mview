# States management system
# Handles saving and loading of view configurations and parameter values

library(jsonlite)

# ensure .states directory exists
ensure_states_directory <- function() {
  if (!dir.exists(".states")) {
    dir.create(".states")
    cat("created .states directory\n")
  }
}

# get all current parameter values organized by group
get_all_parameter_values <- function() {
  all_params <- get_registered_parameters()
  param_values <- list()
  
  for (param_key in names(all_params)) {
    param <- all_params[[param_key]]
    group <- param$group
    id <- param$id
    current_value <- param$reactive_accessor()  # get current value from reactiveVal
    
    if (!is.element(group, names(param_values))) {
      param_values[[group]] <- list()
    }
    param_values[[group]][[id]] <- current_value
  }
  
  return(param_values)
}

# set all parameter values from loaded state
set_all_parameter_values <- function(param_values) {
  if (is.null(param_values) || length(param_values) == 0) {
    cat("no parameter values to set\n")
    return()
  }
  
  for (group_name in names(param_values)) {
    group_params <- param_values[[group_name]]
    for (param_id in names(group_params)) {
      tryCatch({
        update_parameter_from_ui(group_name, param_id, group_params[[param_id]])
        # cat(sprintf("set parameter %s_%s = %s\n", group_name, param_id, group_params[[param_id]]))
      }, error = function(e) {
        warning(sprintf("failed to set parameter %s_%s: %s", group_name, param_id, e$message))
      })
    }
  }
}

# get current view id
get_current_view_id <- function() {
  current_view <- .view_env$current_view
  if (is.null(current_view)) return(NULL)
  return(current_view$id)
}

# save current state to file
save_state <- function(state_number, title) {
  if (!is.numeric(state_number) || state_number < 1 || state_number > 9) {
    stop("state_number must be between 1 and 9")
  }
  
  if (is.null(title) || nchar(title) == 0) {
    stop("title cannot be empty")
  }
  
  ensure_states_directory()
  
  # collect current state
  state_data <- list(
    metadata = list(
      title = title,
      timestamp = Sys.time(),
      mview_version = "1.01",
      state_number = state_number
    ),
    view = get_current_view_id(),
    parameters = get_all_parameter_values()
  )
  
  # save to file
  filename <- file.path(".states", sprintf("state_%d.json", state_number))
  
  tryCatch({
    json_str <- toJSON(state_data, pretty = TRUE, auto_unbox = TRUE)
    writeLines(json_str, filename)
    cat(sprintf("saved state %d: '%s' to %s\n", state_number, title, filename))
    return(TRUE)
  }, error = function(e) {
    warning(sprintf("failed to save state %d: %s", state_number, e$message))
    return(FALSE)
  })
}

# load state from file
load_state <- function(state_number) {
  if (!is.numeric(state_number) || state_number < 1 || state_number > 9) {
    stop("state_number must be between 1 and 9")
  }
  
  filename <- file.path(".states", sprintf("state_%d.json", state_number))
  
  if (!file.exists(filename)) {
    warning(sprintf("state file %s does not exist", filename))
    return(FALSE)
  }
  
  tryCatch({
    # load and parse json
    json_str <- readLines(filename, warn = FALSE)
    state_data <- fromJSON(paste(json_str, collapse = "\n"))
    
    # validate basic structure
    if (is.null(state_data$metadata) || is.null(state_data$metadata$title)) {
      warning(sprintf("invalid state file format in %s", filename))
      return(FALSE)
    }
    
    cat(sprintf("loading state %d: '%s'\n", state_number, state_data$metadata$title))
    
    # set view if specified and exists
    if (!is.null(state_data$view)) {
      current_views <- get_view_ids()
      if (state_data$view %in% current_views) {
        set_view(state_data$view)
        cat(sprintf("set view to: %s\n", state_data$view))
      } else {
        warning(sprintf("view '%s' no longer exists, skipping view change", state_data$view))
      }
    }
    
    # set parameters
    if (!is.null(state_data$parameters)) {
      set_all_parameter_values(state_data$parameters)
    }
    
    
    cat(sprintf("loaded state %d successfully\n", state_number))
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("failed to load state %d: %s", state_number, e$message))
    return(FALSE)
  })
}

# get state title from file (for UI display)
get_state_title <- function(state_number) {
  if (!is.numeric(state_number) || state_number < 1 || state_number > 9) {
    return("")
  }
  
  filename <- file.path(".states", sprintf("state_%d.json", state_number))
  
  if (!file.exists(filename)) {
    return("")
  }
  
  tryCatch({
    json_str <- readLines(filename, warn = FALSE)
    state_data <- fromJSON(paste(json_str, collapse = "\n"))
    
    if (!is.null(state_data$metadata) && !is.null(state_data$metadata$title)) {
      return(state_data$metadata$title)
    } else {
      return("")
    }
  }, error = function(e) {
    return("")
  })
}

# check if state file exists
state_exists <- function(state_number) {
  if (!is.numeric(state_number) || state_number < 1 || state_number > 9) {
    return(FALSE)
  }
  
  filename <- file.path(".states", sprintf("state_%d.json", state_number))
  return(file.exists(filename))
}
