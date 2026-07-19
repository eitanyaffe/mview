# Load shared alignment API
source("profiles/align/align_profile_api.r")

# Load interval profile
source("profiles/interval_profile.r")

# Load variants profile
source("profiles/variants_profile.r")

# Load rearrangements profile
source("profiles/rearrangements_profile.r")

# Load rRNA profile
source("profiles/rrna_profile.r")

# Load synteny profile
source("profiles/synteny/synteny_profile.r")
source("profiles/synteny/synteny_profile_detail.r") 
source("profiles/synteny/synteny_profile_summary.r")

# Load differing variants profile
source("profiles/differing_variants_profile.r")

# Load CME predicted sites profile
source("profiles/sites_profile.r")

# Initialize alntools
init_alntools(verbose = FALSE)

########################################################
# Sample pair parameters (always registered, used by differing variants in all views)
########################################################

sample_choices <- sort(unique(sample.map$sample))

rv_s1 <- register_param("aln_samples", "sample_1", "select", sample_choices[1], choices = sample_choices)
rv_s2 <- register_param("aln_samples", "sample_2", "select", sample_choices[2], choices = sample_choices)

########################################################
# Alignment profiles
########################################################

# self-align profile
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

view_type <- get_current_view_parameter("view_type")

if (identical(view_type, "any_pair")) {
  make_pair_profile <- function(rv_s, pid) {
    align_profile(
      id   = paste0("align_", pid),
      name = rv_s(),
      attr = list(title_f = function() rv_s()),
      aln_f = function() {
        row <- sample.map[sample.map$sample == rv_s(), ]
        if (nrow(row) == 0) return(NULL)
        tag <- paste0("ALL_", row$lib[1])
        get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
      },
      params = default_alignment_params
    )
  }

  make_pair_profile(rv_s1, "sample_1")
  make_pair_profile(rv_s2, "sample_2")

} else {
  # per-sample alignment profiles for the selected mouse
  make_align_profile <- function(name, lib) {
    tag <- paste0("ALL_", lib)
    align_profile(
      id = paste0("align_", lib),
      name = name,
      aln_f = function() {
        get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
      },
      params = default_alignment_params
    )
  }

  timepoint_samples <- get_current_view_parameter("timepoint_samples")
  if (!is.null(timepoint_samples)) {
    condition_samples <- get_current_view_parameter("condition_samples")
    rows <- sample.map[sample.map$timepoint == timepoint_samples, , drop = FALSE]
    if (!is.null(condition_samples)) {
      rows <- rows[rows$condition == condition_samples, , drop = FALSE]
    }
    for (i in seq_len(nrow(rows))) {
      make_align_profile(as.character(rows$individual[i]), as.character(rows$lib[i]))
    }
  } else {
    individual <- get_current_view_parameter("individual")
    if (!is.null(individual)) {
      mouse_samples <- sample.map[sample.map$individual == individual, ]
      mouse_samples$tp_rank <- match(mouse_samples$timepoint, tp_order)
      mouse_samples <- mouse_samples[order(mouse_samples$tp_rank), ]

      for (i in seq_len(nrow(mouse_samples))) {
        make_align_profile(mouse_samples$timepoint[i], mouse_samples$lib[i])
      }
    }
  }
}

########################################################
# basic gene profile
########################################################

gene_profile(
  id = "genes",
  name = "Genes",
  gene_f = get_genes_f,
  color_field = "tax",
  label_field = "label",
  params = default_gene_params
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
# differing variants profile (sample pair from aln_samples params)
########################################################

get_differing_variants_f <- function() {
  s1 <- rv_s1()
  s2 <- rv_s2()

  df <- get_data("UNI_DIST_MERGE_VARIANTS_RAW", null.on.missing = TRUE)
  if (is.null(df)) {
    return(NULL)
  }

  mask <- (df$sample1 == s1 & df$sample2 == s2) | (df$sample1 == s2 & df$sample2 == s1)
  df <- df[mask, ]
  if (nrow(df) == 0) return(NULL)

  # normalize reversed rows so freq1/freq2 always match s1/s2
  reversed <- df$sample1 == s2 & df$sample2 == s1
  if (any(reversed)) {
    tmp_freq     <- df$freq1[reversed];     df$freq1[reversed]     <- df$freq2[reversed];     df$freq2[reversed]     <- tmp_freq
    tmp_support  <- df$support1[reversed];  df$support1[reversed]  <- df$support2[reversed];  df$support2[reversed]  <- tmp_support
    tmp_coverage <- df$coverage1[reversed]; df$coverage1[reversed] <- df$coverage2[reversed]; df$coverage2[reversed] <- tmp_coverage
  }

  df
}

differing_variants_profile(
  id = "differing_variants",
  name = "Diff Variants",
  get_variants_f = get_differing_variants_f
)

########################################################
# axis profile
########################################################

axis_profile()

########################################################
# CME predicted sites profile
########################################################

sites_df <- get_data("SITES_PREDICT_TABLE", null.on.missing = TRUE)

if (!is.null(sites_df)) {
  cme_choices <- sort(unique(as.character(sites_df$cme_id)))
  default_cme <- if (length(cme_choices) > 0) cme_choices[1] else ""

  rv_cme_id <- register_param(
    "sites", "cme_id", "select",
    default = default_cme,
    choices = cme_choices
  )

  sites_profile(
    id   = "sites_predict",
    name = "CME Sites",
    sites_f = function(assembly) sites_df,
    rv_cme_id = rv_cme_id
  )
}
