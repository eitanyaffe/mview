# core/tabs.r
# Tab management system for dynamic tab registration and loading

# Private storage for registered tabs and current loading state
.tabs_env <- new.env(parent = emptyenv())
.tabs_env$registered_tabs <- list()
.tabs_env$current_loading_tab <- NULL
.tabs_env$export_functions <- list()

# Register a new tab with its configuration
register_tab <- function(tab_id, tab_label, tab_code, ...) {
  stopifnot(is.character(tab_id) && length(tab_id) == 1 && nzchar(tab_id))
  stopifnot(is.character(tab_label) && length(tab_label) == 1 && nzchar(tab_label))
  stopifnot(is.character(tab_code) && length(tab_code) == 1 && nzchar(tab_code))
  
  if (!file.exists(tab_code)) {
    stop(sprintf("register_tab: tab code file not found: %s", tab_code))
  }
  
  # Store tab configuration
  tab_config <- list(
    tab_id = tab_id,
    tab_label = tab_label,
    tab_code = tab_code,
    tab_panel_f = NULL,  # Will be set by the tab code itself
    ...
  )
  
  .tabs_env$registered_tabs[[tab_id]] <- tab_config
  cat(sprintf("registering tab: %s, label: %s, code: %s\n", tab_id, tab_label, tab_code))
  invisible(NULL)
}

# Set the tab panel function for the currently loading tab
set_tab_panel_f <- function(tab_panel_f) {
  current_tab <- .tabs_env$current_loading_tab
  if (is.null(current_tab)) {
    stop("set_tab_panel_f can only be called from within a loading tab file")
  }
  
  # Update the registered tab with the panel function
  .tabs_env$registered_tabs[[current_tab$tab_id]]$tab_panel_f <- tab_panel_f
  cat(sprintf("set tab panel for: %s\n", current_tab$tab_id))
  invisible(NULL)
}

# Get all registered tabs
get_registered_tabs <- function() {
  .tabs_env$registered_tabs
}

# Get a specific tab by its ID
get_tab_by_id <- function(tab_id) {
  stopifnot(is.character(tab_id) && length(tab_id) == 1 && nzchar(tab_id))
  
  tabs <- .tabs_env$registered_tabs
  if (tab_id %in% names(tabs)) {
    return(tabs[[tab_id]])
  }
  
  NULL
}

# Load all registered tab code files
load_tabs <- function() {
  tabs <- .tabs_env$registered_tabs
  if (length(tabs) == 0) {
    cat("no tabs registered\n")
    return(invisible(NULL))
  }
  
  for (tab_id in names(tabs)) {
    tab <- tabs[[tab_id]]
    cat(sprintf("loading tab: %s from %s\n", tab_id, tab$tab_code))
    
    # Set current loading tab so get_loaded_tab() works
    .tabs_env$current_loading_tab <- tab
    
    # Source the tab code file into the calling environment (server)
    tryCatch({
      source(tab$tab_code, local = parent.frame())
    }, error = function(e) {
      cat(sprintf("error loading tab %s: %s\n", tab_id, e$message))
    })
    
    # Clear current loading tab
    .tabs_env$current_loading_tab <- NULL
  }
  
  cat(sprintf("loaded %d tabs\n", length(tabs)))
  invisible(NULL)
}

# Get tab panels for UI integration
get_tab_panels <- function() {
  tabs <- .tabs_env$registered_tabs
  panels <- list()
  
  for (tab_id in names(tabs)) {
    tab <- tabs[[tab_id]]
    if (!is.null(tab$tab_panel_f)) {
      # Call the function to create the tabPanel and add as unnamed element
      panels <- c(panels, list(tab$tab_panel_f()))
    }
  }
  
  panels
}

# Register export function for a tab
register_tab_export_function <- function(tab_id, export_type, export_function) {
  stopifnot(is.character(tab_id) && length(tab_id) == 1 && nzchar(tab_id))
  stopifnot(is.character(export_type) && length(export_type) == 1 && nzchar(export_type))
  stopifnot(is.function(export_function))
  stopifnot(export_type %in% c("pdf", "table"))
  
  # Check that the tab exists and supports export
  tab <- get_tab_by_id(tab_id)
  if (is.null(tab)) {
    stop(sprintf("register_tab_export_function: tab '%s' not found", tab_id))
  }
  
  if (!isTRUE(tab$supports_export)) {
    warning(sprintf("register_tab_export_function: tab '%s' does not have supports_export = TRUE", tab_id))
  }
  
  # Initialize tab entry if it doesn't exist
  if (is.null(.tabs_env$export_functions[[tab_id]])) {
    .tabs_env$export_functions[[tab_id]] <- list()
  }
  
  .tabs_env$export_functions[[tab_id]][[export_type]] <- export_function
  cat(sprintf("registered %s export function for tab: %s\n", export_type, tab_id))
  invisible(NULL)
}

# Get export function for a tab
get_tab_export_function <- function(tab_id, export_type) {
  stopifnot(is.character(tab_id) && length(tab_id) == 1 && nzchar(tab_id))
  stopifnot(is.character(export_type) && length(export_type) == 1 && nzchar(export_type))
  stopifnot(export_type %in% c("pdf", "table"))
  
  export_functions <- .tabs_env$export_functions
  if (tab_id %in% names(export_functions) && 
      export_type %in% names(export_functions[[tab_id]])) {
    return(export_functions[[tab_id]][[export_type]])
  }
  
  NULL
}

# Get all tabs that have export functions registered
get_exportable_tabs <- function() {
  export_functions <- .tabs_env$export_functions
  registered_tabs <- .tabs_env$registered_tabs
  
  exportable_tabs <- list()
  for (tab_id in names(export_functions)) {
    if (tab_id %in% names(registered_tabs)) {
      tab <- registered_tabs[[tab_id]]
      if (isTRUE(tab$supports_export)) {
        exportable_tabs[[tab_id]] <- tab
      }
    }
  }
  
  exportable_tabs
}

# Get available export types for a tab
get_tab_export_types <- function(tab_id) {
  stopifnot(is.character(tab_id) && length(tab_id) == 1 && nzchar(tab_id))
  
  export_functions <- .tabs_env$export_functions
  if (tab_id %in% names(export_functions)) {
    return(names(export_functions[[tab_id]]))
  }
  
  character(0)
}