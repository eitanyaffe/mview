########################################################
# Register lookup files
########################################################

# search for files under MAKESHIFT_ROOT/export/long/pb4/
dir <- paste(Sys.getenv("MAKESHIFT_ROOT"), "export", "long", "pb4", "default", sep = "/")
fns <- list.files(dir, full.names = TRUE, pattern = "*.txt")
set_lookup(fns)

########################################################
# Register functions
########################################################

aids <- get_data("ASSEMBLY_TABLE")$ASSEMBLY_ID

# Set assemblies (as a vector of string IDs)
set_assemblies(aids)

# Register contigs function
register_contigs_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(contig = df$contig, length = df$length)
})

# Register genomes function
register_genomes_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_GENOME_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  df$gid <- df$contig
  data.frame(gid = df$gid, length = df$length)
})

# Register contig map function
register_contig_map_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(contig = df$contig, gid = df$contig)
})

# Register views
view_register("long_view", "views/long_view.r")
