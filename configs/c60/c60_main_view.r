########################################################
# Align code
########################################################

library(Rcpp)

alntools_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/alntools/cpp")
cdir <- ".Rcpp_dir"
rebuild <- F
sourceCpp(paste0(alntools_dir, "/aln_R.cpp"), verbose = T, cacheDir = cdir, rebuild = rebuild)

########################################################
# Alignment profiles
########################################################

lib.table <- get_data("LONG_LIB_TABLE")
get_map_tag <- function(assembly, timepoint) {
  ix = lib.table$ASSEMBLY_ID == assembly & lib.table$SAMPLE_TYPE == timepoint
  if (sum(ix) == 0 || sum(ix) > 1) {
    return (NULL)
  }
  paste0(assembly, "_", lib.table$LIB_ID[ix])
}

make_align_profile <- function(timepoint) {
  align_profile(
    id = paste0("align_", timepoint),
    name = timepoint,
    aln_f = function(cxt) {
      tag <- get_map_tag(cxt$assembly, timepoint)
      get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
    },
    height = 0.4,
    params = default_alignment_params
  )
}

timepoints <- get_current_view_parameter("timepoints")
for (timepoint in timepoints) {
  make_align_profile(timepoint)
}

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
  height = 0.05,
  gene_f = get_genes_f,
  params = gene_params_std
)

########################################################
# MGE profile
########################################################

gene_params_mge <- list(
  color_style = list(
    group_id = "mge_genes",
    type = "select",
    choices = c("by_taxonomy", "by_regex"),
    default = "by_regex"
  )
)

mge_gene_groups <- list(
  plasmid = c("plasmid", "conjugation", "conjugative", "mobC"),
  phage = c("tail", "head", "phage", "capsid"),
  mobile = c("mobilization", "transposase", "integrase", "toxin"),
  abx = c("mepA", "efflux", "MATE", "multidrug")
)

mge_gene_colors <- list(
  plasmid = "#29e111", 
  phage = "#ffb300", 
  mobile = "#b1c5ec",
  abx = "red"
)

gene_profile(
  id = "mge_genes",
  name = "MGE",
  height = 0.05,
  gene_f = get_genes_f,
  select_groups = mge_gene_groups,
  select_colors = mge_gene_colors,
  params = gene_params_mge
)

########################################################
# axis profile
########################################################

axis_profile()
