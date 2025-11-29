########################################################
# Minimal view: genes, alignments and axis
########################################################

########################################################
# alignment setup
########################################################

# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Load interval profile
source("profiles/interval_profile.r")

# Initialize alntools
init_alntools(verbose = FALSE)

########################################################
# alignment profiles
########################################################

# single alignment profile that adapts to current assembly context
align_profile(
  id = "alignments",
  name = "Alignments",
  aln_f = function() {
    get_aln_f(cxt_get_assembly())
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
  height = 30,
  gene_f = get_genes_f,
  color_field = "tax",
  label_field = "label",
  params = default_gene_params
)

########################################################
# regions interval profile
########################################################

interval_profile(
  id = "regions",
  name = "Regions",
  intervals_f = "segments.current_regions",
  merge_adjacent = TRUE
)

########################################################
# axis profile
########################################################

axis_profile()
