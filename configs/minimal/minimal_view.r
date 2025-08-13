########################################################
# Example main view - genes, alignments and axis
########################################################

########################################################
# alignment setup
########################################################

library(Rcpp)

alntools_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/alntools/cpp")
cdir <- ".Rcpp_dir"
rebuild <- F
sourceCpp(paste0(alntools_dir, "/aln_R.cpp"), verbose = T, cacheDir = cdir, rebuild = rebuild)

########################################################
# alignment profiles
########################################################

# single alignment profile that adapts to current assembly context
align_profile(
  id = "alignments",
  name = "Alignments",
  aln_f = function(cxt) {
    # get current assembly from context
    assembly <- cxt$assembly
    aln_path <- file.path("examples", "aln", paste0(assembly, ".aln"))
    
    if (!file.exists(aln_path)) {
      warning(sprintf("alignment file not found for assembly %s: %s", assembly, aln_path))
      return(NULL)
    }
    
    cat(sprintf("loading alignment for assembly: %s\n", assembly))
    read_cached(paste0("aln_", assembly), aln_path, read_f = aln_load)
  },
  height = 400,
  params = default_alignment_params
)

########################################################
# basic gene profile
########################################################

gene_params_std <- list(
  color_style = list(
    group_id = "gene",
    type = "select",
    choices = c("by_taxonomy", "by_regex"),
    default = "by_taxonomy"
  )
)

gene_profile(
  id = "genes",
  name = "Genes",
  height = 50,
  gene_f = get_genes_f,
  params = gene_params_std
)

########################################################
# axis profile
########################################################

axis_profile()

