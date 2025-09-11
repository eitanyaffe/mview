# core/profile_manager.r

# ---- Profile Environment ----

# Global environment to store registered profiles
.profile_env <- new.env(parent = emptyenv())
.profile_env$registered_profiles <- list()

# ---- View Management ----

# global list of registered views
.view_env <- new.env(parent = emptyenv())
.view_env$registered_views <- list()
.view_env$current_view <- NULL

# register a view with id and filename
view_register <- function(id, filename, ...) {
  stopifnot(is.character(id) && length(id) == 1 && nzchar(id))
  stopifnot(is.character(filename) && length(filename) == 1 && nzchar(filename))
  if (!file.exists(filename)) {
    stop(sprintf("view_register: file not found: %s", filename))
  }
  params <- list(...)
  .view_env$registered_views[[id]] <- list(id = id, filename = filename, params = params)
  cat(sprintf("registering view: %s, filename: %s\n", id, filename))
  invisible(NULL)
}

# set the current view by id, clearing profiles and sourcing the view file
set_view <- function(id) {
  views <- .view_env$registered_views
  if (!(id %in% names(views))) {
    stop(sprintf("set_view: view id '%s' not found", id))
  }
  cat(sprintf("Setting view to %s\n", id))
  profiles_clear_all()
  clear_parameters(clear_cache = FALSE)

  .view_env$current_view <- views[[id]]
  view_file <- views[[id]]$filename
  cat(sprintf("Sourcing view file: %s\n", view_file))
  # source the view file to register its profiles
  source(view_file, local = TRUE)

  # Trigger UI update to display parameters
  param_registration_done()

  invisible(NULL)
}

get_view_ids <- function() {
  names(.view_env$registered_views)
}

# return the current view parameters list
get_current_view_parameters <- function() {
  if (is.null(.view_env$current_view)) {
    return(list())
  }
  params <- .view_env$current_view$params
  if (is.null(params)) {
    return(list())
  }
  params
}

# return a specific parameter for the current view
get_current_view_parameter <- function(id) {
  stopifnot(is.character(id) && length(id) == 1 && nzchar(id))
  params <- get_current_view_parameters()
  if (!(id %in% names(params))) {
    return(NULL)
  }
  params[[id]]
}

# ---- Profile Creation ----

#' Create a new profile object
#'
#' @param id Character, unique identifier for the profile.
#' @param name Character, display name for the profile.
#' @param type Character, type of profile (e.g., "points", "lines").
#' @param height Numeric, weight for vertical space allocation.
#' @param attr List, additional profile attributes.
#' @param params List, optionalUI parameters for the profile.
#' @param plot_f Function, adds layers to ggplot object.
#' @param auto_register Logical, if TRUE, registers profile immediately.
#' @param ... Additional arguments passed to profile_register.
#' @return A profile object (list).
profile_create <- function(
    id, name, type, height,
    attr = list(), params = list(), plot_f = NULL, get_annotations = NULL, auto_register = FALSE, ...) {
  # Basic validation
  stopifnot(
    is.character(id) && length(id) == 1 && nzchar(id),
    is.character(name) && length(name) == 1 && nzchar(name),
    is.character(type) && length(type) == 1 && nzchar(type),
    is.numeric(height) && length(height) == 1 && height > 0,
    is.list(attr),
    is.function(plot_f)
  )
  
  # Validate height is specified in pixels (positive integer)
  if (height <= 0) {
    stop(sprintf("profile_create: height must be > 0 pixels, got: %s", height))
  }

  cat(sprintf("registering profile: %s\n", id))
  profile <- list(
    id = id,
    name = name,
    type = type,
    height = height,
    attr = attr,
    plot_f = plot_f,
    get_annotations = get_annotations,
    params = params,
    ...
  )

  # register UI parameters
  if (length(params) > 0) {
    for (i in seq_along(params)) {
      param_id <- names(params)[i]
      param <- params[[i]]
      group_id <- if (is.null(param$group_id)) id else param$group_id
      register_param(
        group = group_id,
        id = param_id,
        type = param$type,
        default_value = param$default,
        choices = param$choices
      )
    }
  }

  # Set defaults for common attributes
  if (is.null(profile$attr$title)) {
    profile$attr$title <- name
  }
  if (is.null(profile$attr$should_plot_grid)) {
    profile$attr$should_plot_grid <- TRUE
  }

  if (auto_register) {
    profile_register(profile)
  }

  return(profile)
}

# ---- Profile Registration and Retrieval ----

#' Register a profile
#' @param profile A profile object.
#' @return The profile object, invisibly.
profile_register <- function(profile) {
  stopifnot(is.list(profile) && !is.null(profile$id) && nzchar(profile$id))
  if (profile$id %in% names(.profile_env$registered_profiles)) {
    warning(paste0("profile_register: Profile with id '", profile$id, "' already exists. Overwriting."))
  }
  .profile_env$registered_profiles[[profile$id]] <- profile
  invisible(profile)
}

#' Get a specific profile by its ID
#' @param id Character, the unique ID of the profile to retrieve.
#' @return The profile object if found, otherwise NULL.
profile_get <- function(id) {
  stopifnot(is.character(id) && length(id) == 1)
  .profile_env$registered_profiles[[id]]
}

#' Get all registered profiles
#' @return A list of all registered profile objects.
profiles_get_all <- function() {
  .profile_env$registered_profiles
}

#' Get a parameter value from profiles of a specific type
#' Searches through all registered profiles for the first one matching the type
#' and returns the specified parameter value, or the default if not found.
#' @param param_name Character, name of the parameter to retrieve.
#' @param profile_type Character, type of profile to search (optional).
#' @param default_value Any, value to return if parameter not found.
#' @return The parameter value if found, otherwise default_value.
profile_get_param <- function(param_name, profile_type = NULL, default_value = NULL) {
  stopifnot(is.character(param_name) && length(param_name) == 1 && nzchar(param_name))
  
  profiles <- .profile_env$registered_profiles
  if (length(profiles) == 0) {
    return(default_value)
  }
  
  for (profile in profiles) {
    # check if profile type matches (if specified)
    if (!is.null(profile_type) && profile$type != profile_type) {
      next
    }
    
    # check if parameter exists in this profile
    if (!is.null(profile[[param_name]])) {
      return(profile[[param_name]])
    }
  }
  
  return(default_value)
}

#' Clear all registered profiles
#' Removes all profiles from the registry.
profiles_clear_all <- function() {
  .profile_env$registered_profiles <- list()
  invisible(NULL)
}
