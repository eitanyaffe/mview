# core/profile_manager.r

# ---- Profile Environment ----

# Global environment to store registered profiles
.profile_env <- new.env(parent = emptyenv())
.profile_env$registered_profiles <- list()

# ---- View Management ----

# global list of registered views
.view_env <- new.env(parent = emptyenv())
.view_env$registered_views <- list()

# register a view with id and filename
view_register <- function(id, filename) {
  stopifnot(is.character(id) && length(id) == 1 && nzchar(id))
  stopifnot(is.character(filename) && length(filename) == 1 && nzchar(filename))
  .view_env$registered_views[[id]] <- list(id = id, filename = filename)
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

# ---- Profile Creation ----

#' Create a new profile object
#'
#' @param id Character, unique identifier for the profile.
#' @param name Character, display name for the profile.
#' @param type Character, type of profile (e.g., "points", "lines").
#' @param height Numeric, weight for vertical space allocation.
#' @param attr List, additional profile attributes.
#' @param params List, optionalUI parameters for the profile.
#' @param data_f Function, returns data for plotting.
#' @param plot_f Function, adds layers to ggplot object.
#' @param auto_register Logical, if TRUE, registers profile immediately.
#' @return A profile object (list).
profile_create <- function(
    id, name, type, height,
    attr = list(), params = list(), data_f = NULL, plot_f = NULL, auto_register = FALSE) {
  # Basic validation
  stopifnot(
    is.character(id) && length(id) == 1 && nzchar(id),
    is.character(name) && length(name) == 1 && nzchar(name),
    is.character(type) && length(type) == 1 && nzchar(type),
    is.numeric(height) && length(height) == 1 && height > 0,
    is.list(attr),
    is.function(data_f),
    is.function(plot_f)
  )
  cat(sprintf("registering profile: %s\n", id))
  profile <- list(
    id = id,
    name = name,
    type = type,
    height = height,
    attr = attr,
    params = params,
    data_f = data_f,
    plot_f = plot_f
  )

  # register UI parameters
  if (length(params) > 0) {
    for (i in seq_along(params)) {
      param_id <- names(params)[i]
      param <- params[[i]]
      register_param(
        group = id,
        id = param_id,
        type = param$type,
        default_value = param$default
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

#' Clear all registered profiles
#' Removes all profiles from the registry.
profiles_clear_all <- function() {
  .profile_env$registered_profiles <- list()
  invisible(NULL)
}
