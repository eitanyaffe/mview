# Cache module for storing computation results in memory
# Primarily designed for caching large data frames

# Cache initialization
cache_init <- function(clean = FALSE) {
  if (!exists("g.cache") || clean) {
    g.cache <<- list()
  }
}

# Clean the cache (remove all entries)
cache_clean <- function() {
  cache_init(clean = TRUE)
}

# Check if a key exists in the cache
cache_exists <- function(key) {
  is.element(key, names(g.cache))
}

# Get a value from the cache
cache_get <- function(key) {
  if (!cache_exists(key)) {
    stop(sprintf("key %s not found in cache\n", key))
  }
  g.cache[[key]]
}

# Set a value in the cache
cache_set <- function(key, value) {
  cache_init()
  g.cache[[key]] <<- value
}

# Remove a value from the cache
cache_unset <- function(key) {
  cache_init()
  g.cache[[key]] <<- NULL
}

# Main cache function: evaluates expression or retrieves cached result
cache <- function(key, expr) {
  cache_init()
  if (cache_exists(key)) {
    return(cache_get(key))
  }

  # Evaluate the expression
  result <- eval(expr)

  # Store in cache
  cache_set(key, result)

  return(result)
}

# Get information about the cache contents
cache_info <- function() {
  cache_init()
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
cache_list <- function() {
  for (key in names(g.cache)) {
    item <- g.cache[[key]]
    size <- object.size(item)
    size_str <- format_size(size)
    cat(sprintf("%s: %s (%s)\n", key, typeof(item), size_str))
  }
}

# cache test
cache_test <- function() {
  cache_init()
  cache_set("test", 1:10)
  cache_get("test")
  cache_list()
}

cache_init()
