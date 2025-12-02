# segment color registry
# pure registry functions, no Shiny dependencies

# private storage
.seg_colors_env <- new.env(parent = emptyenv())
.seg_colors_env$color_scheme <- "order"
.seg_colors_env$user_schemes <- list()
.seg_colors_env$color_mappings <- list()

# dynamic schemes (always available)
.seg_colors_env$dynamic_schemes <- c("order", "degree")

# register user-defined color schemes
# schemes: named list, e.g. list(bin = "bin_color", genome = "genome_color")
register_segment_colors <- function(schemes) {
  if (!is.list(schemes)) {
    stop("schemes must be a named list")
  }
  for (name in names(schemes)) {
    cat(sprintf("registering segment color scheme: %s -> %s\n", name, schemes[[name]]))
    .seg_colors_env$user_schemes[[name]] <- schemes[[name]]
  }
}

# register a color mapping function for a style
# map_f: function(ids) -> colors (vectorized)
register_color_mapping_f <- function(style, map_f) {
  .seg_colors_env$color_mappings[[style]] <- map_f
}

# get color mapping function for a style
get_color_mapping_f <- function(style) {
  .seg_colors_env$color_mappings[[style]]
}

# get the current color mapping function (uses current style)
# returns function(ids) -> colors (vectorized)
get_current_color_map <- function() {
  style <- get_segment_color_scheme()
  map_f <- get_color_mapping_f(style)
  if (is.null(map_f)) {
    return(function(ids) rep("#E0E0E0", length(ids)))
  }
  map_f
}

# get all available color schemes
get_segment_color_schemes <- function() {
  c(.seg_colors_env$dynamic_schemes, names(.seg_colors_env$user_schemes))
}

# set current color scheme
set_segment_color_scheme <- function(name) {
  available <- get_segment_color_schemes()
  if (!name %in% available) {
    warning(sprintf("unknown color scheme: %s, available: %s", name, paste(available, collapse = ", ")))
    return(FALSE)
  }
  .seg_colors_env$color_scheme <- name
  TRUE
}

# get current color scheme
get_segment_color_scheme <- function() {
  .seg_colors_env$color_scheme
}

# get user schemes (for server to look up field names)
get_user_color_schemes <- function() {
  .seg_colors_env$user_schemes
}

# check if scheme is dynamic
is_dynamic_scheme <- function(name) {
  name %in% .seg_colors_env$dynamic_schemes
}

