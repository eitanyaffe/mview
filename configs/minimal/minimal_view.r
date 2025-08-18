########################################################
# Minimal view: genes, alignments and axis
########################################################

########################################################
# alignment setup
########################################################

# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Initialize alntools with GPU support if available
init_alntools(enable_gpu = TRUE, verbose = FALSE)

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
