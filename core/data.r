# Data module for loading and managing lookups, assemblies, contigs, genomes, etc.

# Private storage for data and registered functions
.data_env <- new.env(parent = emptyenv())
.data_env$lookups <- NULL
.data_env$register_contigs_f <- NULL
.data_env$register_genomes_f <- NULL
.data_env$register_contig_map_f <- NULL
.data_env$register_fasta_f <- NULL
.data_env$assemblies <- NULL
.data_env$regions_dir <- NULL

# Set lookup tables from files
set_lookup <- function(lookup_files) {
  if (length(lookup_files) == 0) {
    stop("no lookup files provided")
  }

  # Load and combine all lookup files
  lookup_tables <- list()
  for (file in lookup_files) {
    if (!file.exists(file)) {
      stop(sprintf("lookup file not found: %s", file))
    }
    table <- read.delim(file, stringsAsFactors = FALSE)
    if (!all(c("id", "path") %in% colnames(table))) {
      stop(sprintf("lookup file missing required columns (id, path): %s", file))
    }
    lookup_tables[[length(lookup_tables) + 1]] <- table
  }

  # Combine all tables
  combined_lookup <- do.call(rbind, lookup_tables)

  # print number of rows
  cat(sprintf("number of lookup files: %d\n", nrow(combined_lookup)))

  # Store in the module environment
  .data_env$lookups <- combined_lookup
}

list_lookup <- function() {
  # print lookup ids
  cat(sprintf("lookup ids:\n%s\n", paste(.data_env$lookups$id, collapse = "\n")))
}

# Get data by ID with optional tag and custom read function
get_data <- function(id, tag = "", read_f = read.delim, null.on.missing = FALSE) {
  if (!exists("lookups", envir = .data_env) || is.null(.data_env$lookups)) {
    stop("lookup table not set, call set_lookup first")
  }

  id <- if (tag != "") paste(id, tag, sep = ":") else id

  # Check if ID exists in lookup table
  if (!id %in% .data_env$lookups$id) {
    if (null.on.missing) {
      return(NULL)
    } else {
      stop(sprintf("data id not found: %s", id))
    }
  }

  # Get path from lookup table
  path <- .data_env$lookups$path[.data_env$lookups$id == id]

  # Use the cache module to load or retrieve from cache
  cache(path, {
    # Check if file exists
    if (!file.exists(path)) {
      stop(sprintf("data file not found: %s", path))
    }
    cat(sprintf("reading %s data from %s\n", id, path))
    read_f(path)
  })
}

make_get_data_f <- function(id, tag = "", read_f = read.delim) {
  function(cxt) {
    get_data(id, tag, read_f)
  }
}

# Clear data from cache
clear_data <- function(id, tag = "") {
  id <- if (tag != "") paste(id, tag, sep = ":") else id
  cache_key <- paste0("DATA_", id)
  if (cache_exists(cache_key)) {
    cache_unset(cache_key)
    return(TRUE)
  } else {
    warning(sprintf("no cached data found for id: %s", cache_key))
    return(FALSE)
  }
}

# Register functions for contigs, genomes, and contig mapping
register_contigs_f <- function(f) {
  cat(sprintf("registering contigs function\n"))
  .data_env$register_contigs_f <- f
}

register_genomes_f <- function(f) {
  cat(sprintf("registering genomes function\n"))
  .data_env$register_genomes_f <- f
}

register_contig_map_f <- function(f) {
  cat(sprintf("registering contig map function\n"))
  .data_env$register_contig_map_f <- f
}

register_fasta_f <- function(f) {
  cat(sprintf("registering fasta function\n"))
  .data_env$register_fasta_f <- f
}

# Get contigs using the registered function
get_contigs <- function(assembly = NULL) {
  if (is.null(.data_env$register_contigs_f)) {
    stop("contigs function not registered, call register_contigs_f first")
  }
  key <- sprintf("global.contigs.%s", if (is.null(assembly)) "default" else assembly)
  # cat(sprintf("getting contigs for assembly: %s\n", key))
  rr <- cache(key, {
    df <- .data_env$register_contigs_f(assembly)
    df <- df[order(df$length, decreasing = TRUE), ]
  })
  rr
}

# Get genomes using the registered function
get_genomes <- function(assembly = NULL) {
  if (is.null(.data_env$register_genomes_f)) {
    stop("genomes function not registered, call register_genomes_f first")
  }
  key <- sprintf("global.genomes.%s", if (is.null(assembly)) "default" else assembly)
  cache(key, {
    .data_env$register_genomes_f(assembly)
  })
}

# Get contig map using the registered function
get_contig_map <- function(assembly = NULL) {
  if (is.null(.data_env$register_contig_map_f)) {
    stop("contig map function not registered, call register_contig_map_f first")
  }
  key <- sprintf("global.contig_map.%s", if (is.null(assembly)) "default" else assembly)
  cache(key, {
    .data_env$register_contig_map_f(assembly)
  })
}

# Set and get assemblies, which is a vector of assembly string IDs
set_assemblies <- function(assemblies) {
  .data_env$assemblies <- assemblies
}

get_assemblies <- function() {
  .data_env$assemblies
}

# Set and get regions directory
set_regions_dir <- function(regions_dir) {
  .data_env$regions_dir <- regions_dir
}

get_regions_dir <- function() {
  if (is.null(.data_env$regions_dir)) {
    stop("regions directory not set, call set_regions_dir first")
  }
  .data_env$regions_dir
}

# get fasta using the registered function
get_fasta <- function(assembly = NULL) {
  if (is.null(.data_env$register_fasta_f)) {
    stop("fasta function not registered, call register_fasta_f first")
  }
  key <- sprintf("global.fasta.%s", if (is.null(assembly)) "default" else assembly)
  cache(key, {
    .data_env$register_fasta_f(assembly)
  })
}
