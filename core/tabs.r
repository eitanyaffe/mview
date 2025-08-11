# core/tabs.r
# Tab management system for dynamic tab registration and loading

# Private storage for registered tabs and current loading state
.tabs_env <- new.env(parent = emptyenv())
.tabs_env$registered_tabs <- list()
.tabs_env$current_loading_tab <- NULL

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