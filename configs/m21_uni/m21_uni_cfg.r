########################################################
# Set regions directory
########################################################

set_regions_dir("configs/m21_uni/regions")

########################################################
# Register lookup files
########################################################

project_name <- "m21_uni"

dir <- paste(Sys.getenv("MAKESHIFT_ROOT"), "export", "long", "m21", "default", sep = "/")
fns <- list.files(dir, full.names = TRUE, pattern = "*.txt")
set_lookup(fns)

########################################################
# set default navigation mode
########################################################

navigate_up_down_type <- "genome"
if (!navigate_up_down_type %in% c("genome", "contig")) {
  stop(sprintf("invalid navigate_up_down_type: %s", navigate_up_down_type))
}
cache_set("navigate_up_down_type", navigate_up_down_type)

########################################################
# Register assemblies, genomes and contigs
########################################################

aids <- c("ALL")

set_assemblies(aids)

# Register contigs function
register_contigs_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(contig = df$contig, length = df$length, coverage = df$coverage, circular = df$circular)
})

# Register segments function
register_segments_f(function(assembly = NULL) {
  df <- get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(
    segment = df$segment,
    contig = df$contig,
    start = df$start,
    end = df$end,
    length = df$length,
    stringsAsFactors = FALSE
  )
})

# Register genomes function
register_genomes_f(function(assembly = NULL) {
  df <- get_data("UNI_GENOME_TABLE_FULL")
  if (is.null(df)) return(NULL)

  rr <- data.frame(gid = df$bin, df[, -match(c("bin", "domain"), names(df))])
  rr <- rr[order(-rr$mean_abundance), ]
  return(rr)
})

# Register segment map function: combines b-prefix bins (BINNING_BIN_SEGMENT_TABLE)
# and m-prefix bins (UNI_BIN_SEGMENTS) to cover all unified genome types
register_segment_map_f(function(assembly = NULL) {
    binning_df <- get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly, null.on.missing = TRUE)
    uni_df     <- get_data("UNI_BIN_SEGMENTS", null.on.missing = TRUE)

    parts <- list()
    if (!is.null(binning_df)) {
        parts[["binning"]] <- data.frame(segment = binning_df$segment, gid = binning_df$bin, stringsAsFactors = FALSE)
    }
    if (!is.null(uni_df)) {
        parts[["uni"]] <- data.frame(segment = uni_df$segment, gid = uni_df$bin, stringsAsFactors = FALSE)
    }

    if (length(parts) == 0) return(NULL)
    unique(do.call(rbind, parts))
})

read_fasta_f <- function(path) {
  seqinr::read.fasta(file = path, seqtype = "DNA", as.string = TRUE)
}

get_fasta_f <- function(assembly = NULL) {
  get_data("ASSEMBLY_CONTIG_FILE", tag = assembly, read_f = read_fasta_f)
}

# Register fasta function
register_fasta_f(get_fasta_f)

# Register seg_bins function (segment to bin mapping with coordinates and colors)
register_seg_bins_f(function(assembly) {
  seg_table <- get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly, null.on.missing = TRUE)
  if (is.null(seg_table)) return(NULL)
  
  host_table <- get_data("BINNING_HOST_TABLE", tag = assembly, null.on.missing = TRUE)
  
  if (is.null(host_table)) {
    unique_bins <- sort(unique(seg_table$bin))
    bin_colors <- setNames(rainbow(length(unique_bins), s = 0.5, v = 0.9), unique_bins)
    seg_table$bin_color <- bin_colors[seg_table$bin]
    return(seg_table)
  }
  
  host_bins <- host_table$bin
  
  tax_cols <- c("phylum", "class", "order", "family", "genus", "species")
  available_cols <- tax_cols[tax_cols %in% names(host_table)]
  if (length(available_cols) > 0) {
    host_table <- host_table[do.call(order, host_table[available_cols]), ]
    host_bins <- host_table$bin
  }
  
  host_colors <- setNames(rainbow(length(host_bins), s = 0.6, v = 0.85), host_bins)
  
  all_bins <- unique(seg_table$bin)
  non_host_bins <- all_bins[!all_bins %in% host_bins]
  non_host_colors <- setNames(rep("#CCCCCC", length(non_host_bins)), non_host_bins)
  
  bin_colors <- c(host_colors, non_host_colors)
  seg_table$bin_color <- bin_colors[seg_table$bin]
  
  return(seg_table)
})

# Register segment color schemes
register_segment_colors(list(bin = "bin_color"))

# Register seg_adj function
register_seg_adj_f(function(assembly) {
  count_mat <- get_data("BINNING_SEG_ADJ_count", tag = assembly, null.on.missing = TRUE)
  total_mat <- get_data("BINNING_SEG_ADJ_total_read_count", tag = assembly, null.on.missing = TRUE)
  assoc_mat <- get_data("BINNING_SEG_ADJ_associated_read_count", tag = assembly, null.on.missing = TRUE)
  
  if (!is.null(count_mat)) {
    names(count_mat)[names(count_mat) == "seg_src"] <- "src"
    names(count_mat)[names(count_mat) == "seg_tgt"] <- "tgt"
  }
  if (!is.null(total_mat)) {
    names(total_mat)[names(total_mat) == "seg_src"] <- "src"
    names(total_mat)[names(total_mat) == "seg_tgt"] <- "tgt"
  }
  if (!is.null(assoc_mat)) {
    names(assoc_mat)[names(assoc_mat) == "seg_src"] <- "src"
    names(assoc_mat)[names(assoc_mat) == "seg_tgt"] <- "tgt"
  }
  
  list(count = count_mat, total = total_mat, associated = assoc_mat)
})

# register csegments function
register_csegments_f(function(assembly = NULL) {
  df <- get_data("BINNING_CSEG_TABLE", null.on.missing = TRUE)
  map <- get_data("BINNING_CSEG_MAP", null.on.missing = TRUE)
  map  = map[map$assembly == assembly, ]
  df = df[is.element(df$cseg, map$csegment), ]
  
  if (is.null(df)) return(NULL)
  data.frame(
    csegment = df$cseg,
    seg_cluster = df$seg_cluster,
    count = df$count,
    length = df$length,
    stringsAsFactors = FALSE
  )
})

# register cluster mapping function
register_cluster_mapping_f(function(assembly = NULL) {
  get_data("BINNING_ASS_CLUSTER_MAPPING", tag = assembly, null.on.missing = TRUE)
})

# register cseg adjacency function
register_cseg_adj_f(function(assembly) {
  count_mat <- get_data("BINNING_CSEG_ADJ_count", tag = assembly, null.on.missing = TRUE)
  total_mat <- get_data("BINNING_CSEG_ADJ_total_read_count", tag = assembly, null.on.missing = TRUE)
  assoc_mat <- get_data("BINNING_CSEG_ADJ_associated_read_count", tag = assembly, null.on.missing = TRUE)
  
  if (!is.null(count_mat)) {
    names(count_mat)[names(count_mat) == "cluster_src"] <- "src"
    names(count_mat)[names(count_mat) == "cluster_tgt"] <- "tgt"
  }
  if (!is.null(total_mat)) {
    names(total_mat)[names(total_mat) == "cluster_src"] <- "src"
    names(total_mat)[names(total_mat) == "cluster_tgt"] <- "tgt"
  }
  if (!is.null(assoc_mat)) {
    names(assoc_mat)[names(assoc_mat) == "cluster_src"] <- "src"
    names(assoc_mat)[names(assoc_mat) == "cluster_tgt"] <- "tgt"
  }
  
  list(count = count_mat, total = total_mat, associated = assoc_mat)
})

########################################################
# register views
########################################################

sample.map <- get_data("UNI_SAMPLE_MAP")

tp_order <- c("early", "pre", "post", "post1", "post2", "late")

view_file <- "configs/m21_uni/m21_uni_main_view.r"

view_register("any_pair", view_file, view_type = "any_pair", show_self_align = FALSE)

tps <- unique(as.character(sample.map$timepoint))
tps <- c(intersect(tp_order, tps), setdiff(tps, tp_order))
for (tp in tps) {
  rows <- sample.map[sample.map$timepoint == tp, , drop = FALSE]
  if (nrow(rows) == 0) {
    next
  }
  view_register(sprintf("%s (samples)", tp), view_file, timepoint_samples = tp, show_self_align = FALSE)
}

conds <- sort(unique(as.character(sample.map$condition)))
for (cond in conds) {
  for (tp in tps) {
    rows <- sample.map[sample.map$timepoint == tp & sample.map$condition == cond, , drop = FALSE]
    if (nrow(rows) == 0) {
      next
    }
    view_register(paste0(cond, "_", tp), view_file,
                  timepoint_samples = tp, condition_samples = cond, show_self_align = FALSE)
  }
}

individuals <- unique(sample.map[, c("individual", "condition")])
individuals <- individuals[order(individuals$individual), ]

for (i in seq_len(nrow(individuals))) {
  mouse <- individuals$individual[i]
  condition <- individuals$condition[i]
  view_name <- paste0(condition, "_", mouse)
  view_register(view_name, view_file, individual = mouse, show_self_align = FALSE)
}

view_register("self-align", view_file, individual = NULL, show_self_align = TRUE)

########################################################
# genes
########################################################

mge_groups <- list(
  plasmid = c("plasmid", "conjugation", "conjugative", "mobC"),
  phage = c("tail", "head", "phage", "capsid"),
  mobile = c("mobilization", "transposase", "integrase", "toxin"),
  abx = c("mepA", "efflux", "MATE", "multidrug")
)

mge_colors <- list(
  plasmid = "#29e111", 
  phage = "#ffb300", 
  mobile = "#b1c5ec",
  abx = "red"
)

get_tax_color <- function(tax_values) {
  unique_tax <- sort(unique(tax_values[!is.na(tax_values) & tax_values != "none"]))
  if (length(unique_tax) > 0) {
    tax_colors <- rainbow(length(unique_tax))
    names(tax_colors) <- unique_tax
    return(tax_colors)
  }
  return(NULL)
}

get_mge_color <- function(prot_desc_vector, mge_groups, mge_colors) {
  result <- rep("#E8E8E8", length(prot_desc_vector))
  
  for (group_name in names(mge_groups)) {
    patterns <- mge_groups[[group_name]]
    group_color <- mge_colors[[group_name]]
    
    for (pattern in patterns) {
      matches <- grepl(pattern, prot_desc_vector, ignore.case = TRUE)
      result[matches] <- group_color
    }
  }
  
  return(result)
}

get_genes_f <- function(assembly) {
  cache(paste0(assembly, "_genes"), {
    genes <- get_data("PRODIGAL_GENE_TABLE", tag = assembly)
    uniref <- get_data("UNIREF_GENE_TAX_TABLE", tag = assembly)
    ix <- match(genes$gene, uniref$gene)
    fields <- c("uniref", "identity", "coverage", "evalue", "bitscore", "prot_desc", "tax", "uniref_count")
    for (field in fields) {
      if (is.numeric(uniref[[field]])) {
        genes[[field]] <- ifelse(is.na(ix), 0, uniref[[field]][ix])
      } else {
        genes[[field]] <- ifelse(is.na(ix), "none", uniref[[field]][ix])
      }
    }
    
    tax_color_map <- get_tax_color(genes$tax)
    genes$tax_color <- tax_color_map[genes$tax]
    
    genes$mge_color <- get_mge_color(genes$prot_desc, mge_groups, mge_colors)
    
    genes$label <- paste0(
      "Gene: ", genes$gene, "\n",
      if (!is.null(genes$prot_desc)) paste0("Description: ", genes$prot_desc, "\n") else "",
      if (!is.null(genes$tax)) paste0("Taxonomy: ", genes$tax) else ""
    )
    
    genes
  })
}

register_tab(
  tab_id = "genes",
  tab_label = "Genes", 
  tab_code = "tabs/gene_tab.r",
  get_genes_f = get_genes_f
)

########################################################
# register alignment tab
########################################################

register_tab(
  tab_id = "alignments",
  tab_label = "Alignments",
  tab_code = "tabs/alignment_tab.r"
)

########################################################
# register variants tab
########################################################

get_map_tag <- function(assembly, lib_id) {
  paste0(assembly, "_", lib_id)
}

get_aln_f <- function(assembly, library_id) {
  tag <- get_map_tag(assembly, library_id)
  get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
}

get_variants_gene_table_f <- function(assembly) {
  genes_data <- get_genes_f(assembly)

  data.frame(
    gene = genes_data$gene,
    contig = genes_data$contig,
    start = genes_data$start,
    end = genes_data$end,
    strand = genes_data$strand,
    desc = genes_data$prot_desc,
    stringsAsFactors = FALSE
  )
}

# build library_id_map from sample map for variant trajectory display
build_library_id_map <- function(smap, tp_order, aid) {
  smap$tp_rank <- match(smap$timepoint, tp_order)
  smap$tp_rank[is.na(smap$tp_rank)] <- length(tp_order) + 1

  sets <- list()

  # all-samples set: sorted by condition, then individual, then timepoint
  all_rows <- smap[order(smap$condition, smap$individual, smap$tp_rank), ]
  all_entry <- data.frame(aid = aid, set_name = "all", stringsAsFactors = FALSE)
  all_entry$set_ids <- list(all_rows$sample)
  sets[[length(sets) + 1]] <- all_entry

  # per-individual sets (temporal trajectory per mouse)
  for (ind in sort(unique(smap$individual))) {
    rows <- smap[smap$individual == ind, ]
    rows <- rows[order(rows$tp_rank), ]
    cond <- rows$condition[1]
    entry <- data.frame(aid = aid, set_name = paste0(cond, "_", ind), stringsAsFactors = FALSE)
    entry$set_ids <- list(rows$sample)
    sets[[length(sets) + 1]] <- entry
  }

  # per-condition sets (all samples for a treatment)
  for (cond in sort(unique(smap$condition))) {
    rows <- smap[smap$condition == cond, ]
    rows <- rows[order(rows$individual, rows$tp_rank), ]
    entry <- data.frame(aid = aid, set_name = cond, stringsAsFactors = FALSE)
    entry$set_ids <- list(rows$sample)
    sets[[length(sets) + 1]] <- entry
  }

  do.call(rbind, sets)
}

variant_library_id_map <- build_library_id_map(sample.map, tp_order, aids[1])

register_tab(
  tab_id = "variants",
  tab_label = "Variants",
  tab_code = "tabs/variants/variants_tab.r",
  is.dynamic = FALSE,
  library_id_map = variant_library_id_map,
  sample_map = sample.map[, c("sample", "individual")],
  get_variants_table_f = function(assembly) {
    vars_table <- get_data("UNI_VARIANTS_TABLE")
    vars_genic <- get_data("UNI_VARIANTS_GENIC", null.on.missing = TRUE)
    if (is.null(vars_table)) {
      return(NULL)
    }
    if (is.null(vars_genic)) {
      df <- vars_table
    } else {
      # left join by variant_id: one row per variants row; first genic hit per id (match pattern)
      df <- vars_table
      genic_cols <- setdiff(names(vars_genic), "variant_id")
      genic_cols <- genic_cols[!genic_cols %in% names(df)]
      ix <- match(df$variant_id, vars_genic$variant_id)
      for (col in genic_cols) {
        df[[col]] <- vars_genic[[col]][ix]
      }
    }
    df$row_id <- df$variant_id
    df
  },
  get_variants_support_f = function(assembly) get_data("UNI_VARIANTS_SUPPORT"),
  get_variants_coverage_f = function(assembly) get_data("UNI_VARIANTS_COVERAGE"),
  use_genes = TRUE,
  supports_export = TRUE
)

########################################################
# register rearrangements tab
########################################################

register_tab(
  tab_id = "rearrangements",
  tab_label = "Rearrangements",
  tab_code = "tabs/rearrangements/rearrangements_tab.r",
  is.dynamic = FALSE,
  library_id_map = variant_library_id_map,
  sample_map = sample.map[, c("sample", "individual")],
  get_rearrange_events_f = function(assembly) get_data("UNI_REARRANGE_TABLE"),
  get_rearrange_support_f = function(assembly) get_data("UNI_REARRANGE_SUPPORT"),
  get_rearrange_coverage_f = function(assembly) get_data("UNI_REARRANGE_COVERAGE"),
  use_genes = TRUE,
  supports_export = TRUE
)

########################################################
# register associations tab
########################################################

get_associations_library_ids <- function(assembly) {
  lib_table <- get_data("BINNING_ASSOCIATIONS_LIBRARY_TABLE", tag = assembly, null.on.missing = TRUE)
  if (is.null(lib_table)) {
    return(c("early", "pre", "post", "late"))
  }
  lib_table <- lib_table[lib_table$ASSEMBLY_ID == assembly, ]
  if (nrow(lib_table) == 0) {
    return(c("early", "pre", "post", "late"))
  }
  return(lib_table$LIB_ID)
}

register_tab(
  tab_id = "associations",
  tab_label = "Associations",
  tab_code = "tabs/associations/associations_tab.r",
  get_abundance_f = function(assembly) get_data("BINNING_ABUNDANCE_LR_MAT", tag = assembly, null.on.missing = TRUE),
  get_abundance_summary_f = function(assembly) get_data("BINNING_ABUNDANCE_LR_SUMMARY", tag = assembly, null.on.missing = TRUE),
  get_coverage_f = function(assembly) get_data("BINNING_COV_LR_MAT", tag = assembly, null.on.missing = TRUE),
  get_library_ids_f = get_associations_library_ids,
  get_bin_adj_total_f = function(assembly) get_data("BINNING_BIN_ADJ_mean_total_read_count", tag = assembly, null.on.missing = TRUE),
  get_bin_adj_associated_f = function(assembly) get_data("BINNING_BIN_ADJ_mean_associated_read_count", tag = assembly, null.on.missing = TRUE),
  get_bin_adj_support_f = function(assembly) get_data("BINNING_BIN_ADJ_mean_support_read_count", tag = assembly, null.on.missing = TRUE),
  get_seg_adj_total_f = function(assembly) get_data("BINNING_SEG_ADJ_total_read_count", tag = assembly, null.on.missing = TRUE),
  get_seg_adj_associated_f = function(assembly) get_data("BINNING_SEG_ADJ_associated_read_count", tag = assembly, null.on.missing = TRUE),
  get_seg_adj_count_f = function(assembly) get_data("BINNING_SEG_ADJ_count", tag = assembly, null.on.missing = TRUE),
  get_bin_segment_table_f = function(assembly) get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly, null.on.missing = TRUE),
  get_host_table_f = function(assembly) get_data("BINNING_HOST_TABLE", tag = assembly, null.on.missing = TRUE)
)

########################################################
# register poly tab
########################################################

#register_tab(
#  tab_id = "poly",
#  tab_label = "Poly",
#  tab_code = "tabs/poly/poly_tab.r",
#  library_ids = c("early", "pre", "post", "late"),
#  get_unify_table_f = function(assembly) get_data("EVO_UNIFY_TABLE_ASSEMBLY", tag = assembly, null.on.missing = TRUE),
#  get_unify_support_f = function(assembly) get_data("EVO_UNIFY_SUPPORT_ASSEMBLY", tag = assembly, null.on.missing = TRUE),
#  get_unify_coverage_f = function(assembly) get_data("EVO_UNIFY_COVERAGE_ASSEMBLY", tag = assembly, null.on.missing = TRUE),
#  get_abundance_f = function(assembly) get_data("BINNING_ABUNDANCE_LR_MAT", tag = assembly, null.on.missing = TRUE)
#)

########################################################
# bin segments profile function
########################################################

get_bin_segments_f <- function(assembly) {
  seg_table <- get_seg_bins(assembly)
  if (is.null(seg_table)) {
    return(NULL)
  }
  
  seg_table$assembly <- assembly
  seg_table$desc <- paste0("Segment: ", seg_table$segment, "\nBin: ", seg_table$bin)
  seg_table$id <- seg_table$segment
  
  return(seg_table)
}
