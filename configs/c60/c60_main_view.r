# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Load RNA profile API (short_R.cpp bridge) and profile definitions
source("profiles/rna/rna_profile_api.r")
source("profiles/rna/rna_cov_profile.r")
source("profiles/rna/rna_gene_profile.r")

# Load interval profile
source("profiles/interval_profile.r")

# Load variants profile
source("profiles/variants_profile.r")

# Load rearrangements profile
source("profiles/rearrangements_profile.r")

# Load rRNA profile
source("profiles/rrna_profile.r")

# Load tRNA profile
source("profiles/trna_profile.r")

# Load synteny profile
source("profiles/synteny/synteny_profile.r")
source("profiles/synteny/synteny_profile_detail.r") 
source("profiles/synteny/synteny_profile_summary.r")

# Load CME predicted sites profile
source("profiles/sites_profile.r")

# Initialize alntools
init_alntools(verbose = FALSE)

# Initialize short (srt_load / srt_query for RNA coverage profiles)
init_short(verbose = FALSE)

########################################################
# Alignment profiles
########################################################

# get_map_tag is now defined in c60_cfg.r

# self-align profile (assembly only, no timepoint)
if (get_current_view_parameter("show_self_align")) {
  align_profile(
    id = "align_self",
    name = "Self Align",
    aln_f = function() {
    get_data("MINIMAP_SELF_ALN", tag = cxt_get_assembly(), read_f = aln_load)
    },
    params = default_alignment_params
  )
}

make_align_profile <- function(timepoint) {
  align_profile(
    id = paste0("align_", timepoint),
    name = timepoint,
    aln_f = function() {
      tag <- get_map_tag(cxt_get_assembly(), timepoint)
      if (is.null(tag)) {
        return(NULL)
      }
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
    aln_f = function() {
      tag <- get_map_tag_other(cxt_get_assembly(), assembly_other, timepoint)
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

# synteny data function - simple user function
get_synteny_f <- function(assembly, field, binsize, hide_self = TRUE) {
  synteny_tag <- paste0(assembly, "_", field, "_", binsize)
  data <- get_data("SYNTENY_ASSEMBLY_TABLE", 
                   tag = synteny_tag,
                   null.on.missing = TRUE)
  if (is.null(data)) {
    return(NULL)
  }

  required_cols <- c("contig", "start", "end")
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    path <- get_path("SYNTENY_ASSEMBLY_TABLE", tag = synteny_tag)
    stop(sprintf("synteny table missing required columns (%s) in %s",
                 paste(missing_cols, collapse = ", "),
                 path))
  }
  
  # filter out self-assembly libraries if hide_self is TRUE and data exists
  if (hide_self) {
    keep_cols <- !grepl(paste0("^", assembly, "_"), colnames(data))
    keep_cols[match(required_cols, colnames(data))] <- TRUE
    data <- data[, keep_cols, drop = FALSE]
  }
  
  return(data)
}

# consensus data function - loads the merged consensus RDS file
get_consensus_f <- function(assembly) {
  data <- get_data("SYNTENY_CONSENSUS_MERGED", 
                   tag = assembly,
                   read_f = readRDS)
  return(data)
}

# summary synteny profile
## synteny_profile(
##   id = "synteny",
##   name = "synteny",
##   synteny_f = get_synteny_f,
##   consensus_f = get_consensus_f,
##   params = default_synteny_params
## )

########################################################
# genes profile
########################################################

fancy_gene_profile(
  id     = "genes",
  name   = "Genes",
  gene_f = get_genes_f
)

########################################################
# rRNA gene profile
########################################################

rrna_profile(
  id = "rrna",
  name = "rRNA",
  get_gff_f = function() {
    get_data("BARRNAP_TABLE", tag = cxt_get_assembly(), null.on.missing = TRUE, read_f = read_gff_f)
  }
)

########################################################
# tRNA gene profile
########################################################

trna_profile(
  id = "trna",
  name = "tRNA",
  get_trna_f = function() {
    get_data("TRNA_TABLE", tag = cxt_get_assembly(), null.on.missing = TRUE)
  }
)

########################################################
# assembly interval profile
########################################################

## interval_profile(
##   id = "assembly_segments",
##   name = "Assembly Segments",
##   intervals_f = function(assembly) {
##     data <- get_data("CAV_REFINE_SEGMENT_TABLE", 
##                      tag = assembly,
##                      null.on.missing = TRUE)
##     data$desc <- data$segment        
##     data$id <- ""
##     data
##   }
## )

########################################################
# regions interval profile
########################################################

interval_profile(
  id = "regions",
  name = "Regions",
  intervals_f = "segments.current_regions",
  merge_adjacent = TRUE,
  height = 60
)

########################################################
# bin segments interval profile
########################################################

interval_profile(
  id = "bin_segments",
  name = "Segments",
  intervals_f = get_bin_segments_f,
  color_f = function(ids) { get_current_color_map()(ids) },
  merge_adjacent = FALSE,
  height = 60
)

########################################################
# malign hotspots interval profile (with window param)
########################################################

local({
  params <- list(
    window = list(
      group_id = "hotspots",
      type     = "integer",
      default  = 100
    )
  )

  plot_f <- function(profile, gg) {
    assembly <- cxt_get_assembly()
    hs <- get_malign_hotspots_f(assembly)
    if (is.null(hs) || nrow(hs) == 0)
      return(list(plot = gg, legends = list()))

    w  <- as.integer(if (!is.null(profile$window)) profile$window else min(hs$window))
    hs <- hs[hs$window == w, ]
    if (nrow(hs) == 0)
      return(list(plot = gg, legends = list()))

    filtered <- cxt_filter_intervals(hs, merge_adjacent = FALSE)
    if (is.null(filtered) || nrow(filtered) == 0)
      return(list(plot = gg, legends = list()))

    xlim       <- cxt_get_xlim()
    hover_text <- paste0(filtered$contig, ": ", filtered$start, "-", filtered$end,
                         "\n", filtered$desc)

    gg <- gg +
      ggplot2::geom_rect(
        data = filtered,
        ggplot2::aes(xmin = gstart, xmax = gend, ymin = -0.15, ymax = 0.15, text = hover_text),
        fill = "#CC3333", color = "black", size = 0.5
      ) +
      ggplot2::geom_text(
        data = filtered,
        ggplot2::aes(x = gstart + diff(xlim) * 0.01, y = 0.25, label = id, text = hover_text),
        color = "black", size = 2, hjust = 1, vjust = 0.5
      ) +
      ggplot2::ylim(-0.5, 0.6)

    list(plot = gg, legends = list())
  }

  profile_create(
    id            = "hotspots",
    name          = "Hotspots",
    type          = "intervals",
    height        = 60,
    is_fixed      = TRUE,
    attr          = list(hide_y_ticks = TRUE),
    params        = params,
    plot_f        = plot_f,
    auto_register = TRUE
  )
})

########################################################
# variants profile
########################################################

variants_profile(
  id = "variants",
  name = "Variants"
)

########################################################
# rearrangements profile
########################################################

rearrangements_profile(
  id = "rearrangements",
  name = "Rearrangements"
)

########################################################
# CME predicted sites profile (must register before axis)
########################################################

sites_cme_id_df <- get_data("SITES_MOTIF_CME_ID_TABLE", null.on.missing = TRUE)

if (!is.null(sites_cme_id_df) && "CME_ID" %in% names(sites_cme_id_df)) {
  cme_choices <- sort(unique(as.character(sites_cme_id_df$CME_ID)))
  default_cme <- if ("c4" %in% cme_choices) "c4" else if (length(cme_choices) > 0) cme_choices[1] else ""

  rv_cme_id <- register_param(
    "sites", "cme_id", "select",
    default = default_cme,
    choices = cme_choices
  )

  sites_profile(
    id   = "sites_predict",
    name = "CME Sites",
    sites_f = function(assembly) get_data("SITES_PREDICT_ANNOTATE_SITES", tag = assembly, null.on.missing = TRUE),
    rv_cme_id = rv_cme_id
  )
}

########################################################
# RNA profiles
# RNA_LIB_TABLE columns: ASSEMBLY_ID, LIB_ID, SAMPLE_TYPE, PAIRED_DNA_LIB_ID
# SHORT_SRT_ALL:<LIB_ID>   → srt_all.srt for any RNA-flow lib (RNA or paired DNA)
# RNA_ASS_GENE_TABLE:<AID> → per-subject gene expression table
# RNA_ASS_LIB_STATS:<AID>  → per-subject lib stats (total_mapped per lib)
########################################################

local({
    rna_lib_table <- get_data("RNA_LIB_TABLE", null.on.missing = TRUE)

    # map style name → export variable key
    .srt_key <- function(style) switch(style,
        "unique_strict"     = "SHORT_SRT_ALL_RNA_UNIQUE_STRICT",
        "non_unique_strict" = "SHORT_SRT_ALL_RNA_NON_UNIQUE_STRICT",
                              "SHORT_SRT_ALL_RNA_NON_UNIQUE_LOOSE")

    # look up RNA srt for (assembly, sample_type, map_style)
    get_rna_srt_f <- function(assembly, sample_type, map_style = "non_unique_loose") {
        if (is.null(rna_lib_table)) return(NULL)
        ix <- rna_lib_table$ASSEMBLY_ID == assembly &
              rna_lib_table$SAMPLE_TYPE  == sample_type
        if (sum(ix) != 1) return(NULL)
        lib_id <- rna_lib_table$LIB_ID[ix]
        get_data(.srt_key(map_style), tag = lib_id, read_f = srt_load,
                 null.on.missing = TRUE)
    }

    # look up paired DNA srt for (assembly, sample_type, map_style)
    get_dna_srt_f <- function(assembly, sample_type, map_style = "non_unique_loose") {
        if (is.null(rna_lib_table)) return(NULL)
        ix <- rna_lib_table$ASSEMBLY_ID == assembly &
              rna_lib_table$SAMPLE_TYPE  == sample_type
        if (sum(ix) != 1) return(NULL)
        dna_lib_id <- rna_lib_table$PAIRED_DNA_LIB_ID[ix]
        get_data(.srt_key(map_style), tag = dna_lib_id, read_f = srt_load,
                 null.on.missing = TRUE)
    }

    get_rna_gene_table_f <- function(assembly) {
        # check.names=FALSE preserves column names like D+7_rna_count (+ would become . otherwise)
        get_data("RNA_ASS_GENE_TABLE", tag = assembly, null.on.missing = TRUE,
                 read_f = function(f) read.delim(f, check.names = FALSE))
    }

    get_rna_lib_stats_f <- function(assembly) {
        get_data("RNA_ASS_LIB_STATS", tag = assembly, null.on.missing = TRUE)
    }

    rna_sample_types <- c("D-1", "D+4", "D+7", "D+18")

    rna_cov_profile(
        id              = "rna_cov",
        name            = "RNA coverage",
        get_rna_srt_f   = get_rna_srt_f,
        get_dna_srt_f   = get_dna_srt_f,
        get_lib_stats_f = get_rna_lib_stats_f,
        sample_types    = rna_sample_types,
        height          = 150,
        shared_id       = "rna_params"
    )

    rna_cov_profile(
        id              = "dna_cov",
        name            = "DNA coverage",
        get_rna_srt_f   = get_rna_srt_f,
        get_dna_srt_f   = get_dna_srt_f,
        get_lib_stats_f = get_rna_lib_stats_f,
        sample_types    = rna_sample_types,
        height          = 150,
        shared_id       = "rna_params",
        mode_fixed      = "dna"
    )

    rna_gene_profile(
        id               = "rna_genes",
        name             = "RNA genes",
        get_gene_table_f = get_rna_gene_table_f,
        get_lib_stats_f  = get_rna_lib_stats_f,
        sample_types     = rna_sample_types,
        height           = 100,
        shared_id        = "rna_params"
    )
})

########################################################
# axis profile (must be last — adds coord_cartesian xlim for bottom panel)
########################################################

axis_profile(height = 80)
