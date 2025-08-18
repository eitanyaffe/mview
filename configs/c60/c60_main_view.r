# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Initialize alntools with GPU support
init_alntools(enable_gpu = TRUE, verbose = FALSE)

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

gene_params <- list(
  color_field = list(
    group_id = "gene",
    type = "select",
    choices = c("tax", "mge"),
    default = "tax"
  ),
  height = list(
    group_id = "gene",
    type = "integer",
    default = 50
  )
)

gene_profile(
  id = "genes",
  name = "Genes",
  height = gene_params$height$default,
  gene_f = get_genes_f,
  color_field = "tax",
  label_field = "label",
  params = gene_params
)

########################################################
# axis profile
########################################################

axis_profile()
