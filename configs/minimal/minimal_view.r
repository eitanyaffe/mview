########################################################
# Minimal view: genes, alignments and axis
########################################################

########################################################
# alignment setup
########################################################

library(Rcpp)

alntools_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/alntools/cpp")
cdir <- ".Rcpp_dir"
rebuild <- F

# Configure OpenMP for threading support in bin queries
if (Sys.info()["sysname"] == "Darwin" && file.exists("/opt/homebrew/opt/libomp/include/omp.h")) {
  # macOS with Homebrew libomp - enables parallel processing of large alignment datasets
  Sys.setenv(PKG_CPPFLAGS = "-I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp -D_OPENMP=200805")
  Sys.setenv(PKG_LIBS = "-L/opt/homebrew/opt/libomp/lib -lomp")
}

sourceCpp(paste0(alntools_dir, "/aln_R.cpp"), verbose = T, cacheDir = cdir, rebuild = rebuild)

########################################################
# alignment profiles
########################################################

# single alignment profile that adapts to current assembly context
align_profile(
  id = "alignments",
  name = "Alignments",
  aln_f = function(cxt) {
    get_aln_f(cxt$assembly)
  },
  height = 400,
  params = default_alignment_params
)

########################################################
# basic gene profile
########################################################

gene_params <- list(
  color_field = list(
    group_id = "gene",
    type = "select",
    choices = c("tax", "mge"),
    default = "tax"
  )
)

gene_profile(
  id = "genes",
  name = "Genes",
  height = 50,
  gene_f = get_genes_f,
  color_field = "tax",
  label_field = "label",
  params = gene_params
)

########################################################
# axis profile
########################################################

axis_profile()
