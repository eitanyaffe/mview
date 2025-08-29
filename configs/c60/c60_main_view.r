# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Load segments profile
source("profiles/segments_profile.r")

# Load synteny profile
source("profiles/synteny/synteny_profile.r")
source("profiles/synteny/synteny_profile_detail.r") 
source("profiles/synteny/synteny_profile_summary.r")

# Initialize alntools
init_alntools(verbose = FALSE)

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
      ifn = "/Users/eitany/work/makeshift-dev/export/long/pb-b20/default/minimap/v1.08/map-hifi/lib/EBC/EBN_-1/align.aln"
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
# synteny profiles
########################################################

# synteny data function
get_synteny_f <- function(cxt, binsize, hide_self = TRUE) {
  # load sequenced bp data
  sequenced_bp_data <- get_data("MINIMAP_SYNTENY_ASSEMBLY_TABLE", 
                               tag = paste0(cxt$assembly, "_sequenced_bp_", binsize))
  
  # load mutation density data  
  mutation_data <- get_data("MINIMAP_SYNTENY_ASSEMBLY_TABLE", 
                           tag = paste0(cxt$assembly, "_median_mutation_density_", binsize))
  
  # filter out self-assembly libraries if hide_self is TRUE
  if (hide_self && !is.null(sequenced_bp_data) && !is.null(mutation_data)) {
    # get library columns (exclude contig, start, end)
    coord_cols <- c("contig", "start", "end")
    lib_cols <- setdiff(colnames(sequenced_bp_data), coord_cols)
    
    # filter out libraries that start with current assembly ID
    assembly_prefix <- paste0(cxt$assembly, "_")
    self_libs <- lib_cols[startsWith(lib_cols, assembly_prefix)]
    keep_libs <- setdiff(lib_cols, self_libs)
    
    if (length(self_libs) > 0) {
      cat(sprintf("hiding %d self-assembly libraries: %s\n", 
                  length(self_libs), 
                  paste(self_libs[seq_len(min(3, length(self_libs)))], collapse = ", ")))
    }
    
    # keep only non-self libraries plus coordinate columns
    keep_cols <- c(coord_cols, keep_libs)
    sequenced_bp_data <- sequenced_bp_data[, keep_cols, drop = FALSE]
    mutation_data <- mutation_data[, keep_cols, drop = FALSE]
  }
  
  return(list(
    sequenced_bp = sequenced_bp_data,
    mutations = mutation_data
  ))
}

# summary synteny profile
synteny_profile(
  id = "synteny",
  name = "synteny",
  synteny_f = get_synteny_f,
  params = default_synteny_params
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
