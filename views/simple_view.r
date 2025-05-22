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

align_profile(
  id = "aln_L1",
  name = "Align L1",
  aln_f = make_get_data_f(
    id = "MINIMAP_LIB_ALN", tag = "BAM_L1",
    read_f = aln_load
  )
)

axis_profile()
