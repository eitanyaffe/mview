# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Load segments profile
source("profiles/segments_profile.r")

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

# self-align profile (assembly only, no timepoint)
if (get_current_view_parameter("show_self_align")) {
  align_profile(
    id = "align_self",
    name = "Self Align",
    aln_f = function(cxt) {
    get_data("MINIMAP_SELF_ALN", tag = cxt$assembly, read_f = aln_load)
    },
    params = default_alignment_params
  )
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

get_map_tag_other <- function(assembly, assembly_other, timepoint) {
  ix = lib.table$ASSEMBLY_ID == assembly_other & lib.table$SAMPLE_TYPE == timepoint
  if (sum(ix) == 0 || sum(ix) > 1) {
    return (NULL)
  }
  paste0(assembly, "_", lib.table$LIB_ID[ix])
}

make_align_profile_other <- function(assembly_other, timepoint) {
  align_profile(
    id = paste0("align_", assembly_other, "_", timepoint),
    name = paste0(assembly_other, " ", timepoint),
    aln_f = function(cxt) {
      tag <- get_map_tag_other(cxt$assembly, assembly_other, timepoint)
      ifn = "/Users/eitany/work/makeshift-dev/export/long/pb-b20/default/minimap/v1.08/map-hifi/lib/EBC/EBU_-21/align.aln"
      aln_load(ifn)
      # get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
    },
    params = default_alignment_params
  )
}

subject_ids <- get_current_view_parameter("other_subject_ids")
for (assembly_other in subject_ids) {
  make_align_profile_other(assembly_other, timepoint)
}

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
