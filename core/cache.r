# Cache module for storing computation results in memory
# Primarily designed for caching large data frames

# Cache initialization
cache_init <- function(project_id = "default", clean = FALSE) {
  if (!exists("g.cache") || clean) {
    g.cache <<- list()
  }
  g.project_id <<- project_id
}

get_project_id <- function() {
  if (exists("g.project_id")) {
    return(g.project_id)
  } else {
    return("default")
  }
}

# Clean the cache (remove all entries)
cache_clear <- function() {
  cache_init(project_id = g.project_id, clean = TRUE)
}

# Get full key with project prefix
make_full_key <- function(key) {
  paste0(g.project_id, "_", key)
}

# Check if a key exists in the cache
cache_exists <- function(key) {
  is.element(make_full_key(key), names(g.cache))
}

# Get a value from the cache
cache_get <- function(key) {
  if (!cache_exists(key)) {
    stop(sprintf("key %s not found in cache\n", key))
  }
  g.cache[[make_full_key(key)]]
}

# Get a value from the cache if it exists, otherwise return default
cache_get_if_exists <- function(key, default) {
  if (cache_exists(key)) {
    cache_get(key)
  } else {
    default
  }
}

# Set a value in the cache
cache_set <- function(key, value) {
  g.cache[[make_full_key(key)]] <<- value
}

# Remove a value from the cache
cache_unset <- function(key) {
  g.cache[[make_full_key(key)]] <<- NULL
}

# Main cache function: evaluates expression or retrieves cached result
cache <- function(key, expr, force = FALSE) {
  if (cache_exists(key) && !force) {
    return(cache_get(key))
  }

  # Evaluate the expression
  result <- eval(substitute(expr), parent.frame())

  # Store in cache
  cache_set(key, result)

  return(result)
}

# Get information about the cache contents
cache_info <- function() {
  cache_size <- length(g.cache)

  if (cache_size == 0) {
    cat("cache is empty\n")
    return(invisible(NULL))
  }

  # Calculate total size in GB
  total_bytes <- 0
  for (item in g.cache) {
    total_bytes <- total_bytes + as.numeric(object.size(item))
  }

  total_gb <- total_bytes / (1024^3)

  cat("cache contains", cache_size, "elements\n")
  cat("total cache size:", sprintf("%.2f GB", total_gb), "\n")
  cat("current project:", g.project_id, "\n")

  invisible(NULL)
}

# format size to kb/mb/gb
format_size <- function(size) {
  if (size < 1024) {
    return(sprintf("%.2f B", size))
  } else if (size < 1024^2) {
    return(sprintf("%.2f KB", size / 1024))
  } else if (size < 1024^3) {
    return(sprintf("%.2f MB", size / (1024^2)))
  } else {
    return(sprintf("%.2f GB", size / (1024^3)))
  }
}

# cache list, print items, their type and sizes
cache_list <- function(prefix = "") {
  project_prefix <- paste0(g.project_id, "_")
  for (key in names(g.cache)) {
    if (startsWith(key, project_prefix) && startsWith(key, paste0(project_prefix, prefix))) {
      item <- g.cache[[key]]
      size <- object.size(item)
      size_str <- format_size(size)
      # Remove project prefix from display
      display_key <- sub(project_prefix, "", key)
      cat(sprintf("%s: %s (%s)\n", display_key, typeof(item), size_str))
    }
  }
}

# cache test
cache_test <- function() {
  cache_init(project_id = "test")
  cache_set("test", 1:10)
  cache_get("test")
  cache_list()
}
