########################################################
# helper functions
########################################################

get_contigs_f <- function(cxt) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = cxt$assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(contig = df$contig, coord = df$length / 2, value = df$length, desc = paste0(df$contig, " (", format_bp(df$length), ")"))
}

########################################################
# Align code
########################################################

library(Rcpp)

alntools_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/alntools/cpp")
cdir <- ".Rcpp_dir"
rebuild <- F
sourceCpp(paste0(alntools_dir, "/aln_R.cpp"), verbose = T, cacheDir = cdir, rebuild = rebuild)

########################################################
# Register profiles
########################################################


get_map_tag <- function(assembly) {
  switch(assembly,
    "BAM" = "L1",
    "EBP" = "L2",
    "EBI" = "L9",
    "EBC" = "L11"
  )
}

align_profile(
  id = "align",
  name = "Align",
  aln_f = function(cxt) {
    tag <- paste0(cxt$assembly, "_", get_map_tag(cxt$assembly))
    get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
  }
)

axis_profile()
