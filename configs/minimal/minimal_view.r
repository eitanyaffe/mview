########################################################
# Minimal view: genes, alignments and axis
########################################################

########################################################
# alignment setup
########################################################

# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Load segments profile
source("profiles/segments_profile.r")

# Initialize alntools
init_alntools(verbose = FALSE)

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

gene_profile(
  id = "genes",
  name = "Genes",
  height = default_gene_params$height$default,
  gene_f = get_genes_f,
  color_field = "tax",
  label_field = "label",
  params = default_gene_params
)

########################################################
# regions segments profile
########################################################

segments_profile(
  id = "regions",
  name = "Regions",
  segments_data = "segments.current_regions"
)

########################################################
# axis profile
########################################################

axis_profile()
