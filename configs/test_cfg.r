# register here proejct specific paths and decisions, including choice of views

########################################################
# Create test data
########################################################

# define contig table with assemby, cid, length
test_contigs_f <- function(N, M) {
  rr <- NULL
  for (i in 1:N) {
    assembly <- paste0("a", i)
    rr <- rbind(rr, data.frame(
      assembly = assembly,
      cid = paste0(assembly, "_c", 1:M),
      length = sample(100:500, M)
    ))
  }
  rr$gid <- paste0("g", 1:N)
  return(rr)
}

# make genome table from contig table
test_genomes_f <- function(contigs) {
  ss <- split(contigs, contigs$assembly)
  rr <- NULL
  for (i in 1:length(ss)) {
    rr <- rbind(rr, data.frame(
      gid = paste0("g", i),
      length = sum(ss[[i]]$length)
    ))
  }
  return(rr)
}

test_contigs <- test_contigs_f(3, 10)
test_genomes <- test_genomes_f(test_contigs)

########################################################
# Register functions
########################################################

# Set assemblies (as a vector of string IDs)
set_assemblies(unique(test_contigs$assembly))

# Register contigs function
register_contigs_f(function(assembly = NULL) {
  test_contigs[test_contigs$assembly == assembly, ]
})

# Register genomes function
register_genomes_f(function(assembly = NULL) {
  test_contigs[test_contigs$assembly == assembly, ]
})

# Register contig map function
register_contig_map_f(function(assembly = NULL) {
  test_contigs[test_contigs$assembly == assembly, ]
})

# Register views
view_register("test_view1", "views/test_view1.r")
view_register("test_view2", "views/test_view2.r")
