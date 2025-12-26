########################################################
# Set regions directory
########################################################

set_regions_dir("configs/c60/regions")

########################################################
# Register lookup files
########################################################

# project_name <- "pb4"
# project_name <- "pb-b20"
# project_name <- "pb-pilot2"
project_name <- "h40"

# search for lookup files
dir <- paste(Sys.getenv("MAKESHIFT_ROOT"), "export", "long", project_name, "default", sep = "/")
fns <- list.files(dir, full.names = TRUE, pattern = "*.txt")
set_lookup(fns)

########################################################
# set default navigation mode
########################################################

navigate_up_down_type <- "genome"
# validate navigation mode
if (!navigate_up_down_type %in% c("genome", "contig")) {
  stop(sprintf("invalid navigate_up_down_type: %s", navigate_up_down_type))
}
# cache navigation mode
cache_set("navigate_up_down_type", navigate_up_down_type)

########################################################
# Register assemblies, genomes and contigs
########################################################

aids <- get_data("ASSEMBLY_TABLE")$ASSEMBLY_ID
aids <- c("AAK", "EAP", "BAA", "BAM", "EBC", "BAH", "EAL", "EAA", "EAV", "EBQ")

# Sort assembly IDs to prioritize EBC
aids <- aids[order(aids != "EBC")]

# Set assemblies (as a vector of string IDs)
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
  # Return required columns: segment, contig, start, end, length
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
    # df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
    df <- get_data("BINNING_HOST_TABLE", tag = assembly)
    
    if (is.null(df)) {
        return(NULL)
    }
    
    rr = data.frame(gid = df$bin, df[,-match(c("bin","domain"), names(df))])
    
    # load sweeps table (global, no tag)
    sweeps_df <- get_data("EVO_HOST_SWEEPS_CLASSIFIED", null.on.missing = TRUE)
    if (!is.null(sweeps_df)) {
        # filter by assembly and join on host_bin = gid
        sweeps_assembly <- sweeps_df[sweeps_df$aid == assembly, ]
        if (nrow(sweeps_assembly) > 0) {
            # merge on gid = host_bin
            rr <- merge(rr, sweeps_assembly, by.x = "gid", by.y = "host_bin", all.x = TRUE, sort = FALSE)
            # rename and select fields
            if ("var_sweep_class_label" %in% names(rr)) {
                rr$sweep_class <- rr$var_sweep_class_label
            }
            if ("num_rearrange" %in% names(rr)) {
                rr$rearrange <- rr$num_rearrange
            }
            if ("num_sweeping_rearrange" %in% names(rr)) {
                rr$sweeping_rearrange <- rr$num_sweeping_rearrange
            }
            if ("num_sweeping_locals" %in% names(rr)) {
                rr$sweeping_locals <- rr$num_sweeping_locals
            }
            if ("num_sweeping_elements" %in% names(rr)) {
                rr$sweeping_elements <- rr$num_sweeping_elements
            }
            if ("num_local_variants" %in% names(rr)) {
                rr$local_variants <- rr$num_local_variants
            }
            if ("num_associated_elements" %in% names(rr)) {
                rr$associated_elements <- rr$num_associated_elements
            }
            # remove original columns we don't need
            cols_to_remove <- c("aid", "host_bin", "num_rearrange", "num_sweeping_rearrange", 
                              "num_sweeping_locals", "num_sweeping_elements", "num_local_variants",
                              "num_associated_elements", "var_sweep_class_label", "var_sweep_class_id")
            cols_to_remove <- cols_to_remove[cols_to_remove %in% names(rr)]
            if (length(cols_to_remove) > 0) {
                rr <- rr[, !names(rr) %in% cols_to_remove]
            }
        }
    }
    
    rr = rr[order(-rr$length),]
    return(rr)
})

# Register segment map function
register_segment_map_f(function(assembly = NULL) {
    df <- get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly)
    if (is.null(df)) {
        return(NULL)
    }
    rr = data.frame(segment = df$segment, gid = df$bin)
    unique(rr)
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
  
  # load host table to determine which bins are hosts
  host_table <- get_data("BINNING_HOST_TABLE", tag = assembly, null.on.missing = TRUE)
  
  if (is.null(host_table)) {
    # fallback: simple rainbow over all bins
    unique_bins <- sort(unique(seg_table$bin))
    bin_colors <- setNames(rainbow(length(unique_bins), s = 0.5, v = 0.9), unique_bins)
    seg_table$bin_color <- bin_colors[seg_table$bin]
    return(seg_table)
  }
  
  # hosts: bins in host_table, sorted by taxonomy
  host_bins <- host_table$bin
  
  # sort hosts by taxonomy: phylum, class, order, family, genus, species
  tax_cols <- c("phylum", "class", "order", "family", "genus", "species")
  available_cols <- tax_cols[tax_cols %in% names(host_table)]
  if (length(available_cols) > 0) {
    host_table <- host_table[do.call(order, host_table[available_cols]), ]
    host_bins <- host_table$bin
  }
  
  # assign rainbow colors to hosts
  host_colors <- setNames(rainbow(length(host_bins), s = 0.6, v = 0.85), host_bins)
  
  # non-hosts get light gray
  all_bins <- unique(seg_table$bin)
  non_host_bins <- all_bins[!all_bins %in% host_bins]
  non_host_colors <- setNames(rep("#CCCCCC", length(non_host_bins)), non_host_bins)
  
  # combine colors
  bin_colors <- c(host_colors, non_host_colors)
  seg_table$bin_color <- bin_colors[seg_table$bin]
  
  return(seg_table)
})

# Register segment color schemes
register_segment_colors(list(bin = "bin_color"))

# Register seg_adj function (segment adjacency matrices, normalize to generic src/tgt)
register_seg_adj_f(function(assembly) {
  count_mat <- get_data("BINNING_SEG_ADJ_count", tag = assembly, null.on.missing = TRUE)
  total_mat <- get_data("BINNING_SEG_ADJ_total_read_count", tag = assembly, null.on.missing = TRUE)
  assoc_mat <- get_data("BINNING_SEG_ADJ_associated_read_count", tag = assembly, null.on.missing = TRUE)
  
  # normalize column names to generic src/tgt
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

# register csegments function (keep csegment column names)
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

# register cseg adjacency function (normalize to generic src/tgt)
register_cseg_adj_f(function(assembly) {
  count_mat <- get_data("BINNING_CSEG_ADJ_count", tag = assembly, null.on.missing = TRUE)
  total_mat <- get_data("BINNING_CSEG_ADJ_total_read_count", tag = assembly, null.on.missing = TRUE)
  assoc_mat <- get_data("BINNING_CSEG_ADJ_associated_read_count", tag = assembly, null.on.missing = TRUE)
  
  # normalize column names to generic src/tgt
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

view_file <- "configs/c60/c60_main_view.r"

view_register("early", view_file, timepoints = c("early"), show_self_align = FALSE)
view_register("pre", view_file, timepoints = c("pre"), show_self_align = FALSE)
view_register("post", view_file, timepoints = c("post"), show_self_align = FALSE)
view_register("late", view_file, timepoints = c("late"), show_self_align = FALSE)
view_register("pre / post", view_file, timepoints = c("pre", "post"), show_self_align = FALSE)
view_register("pre / post / late", view_file, timepoints = c("pre", "post", "late"), show_self_align = FALSE)
view_register("all", view_file, timepoints = c("early", "pre", "post", "late"), show_self_align = FALSE)

view_register("other subject", view_file, other_subject_ids = c("EBN"), timepoints = c("pre"), show_self_align = FALSE)

view_register("self-align", view_file, timepoints = NULL, show_self_align = TRUE)

########################################################
# genes
########################################################

# mge classification groups and colors
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

# helper function to generate taxonomy colors
get_tax_color <- function(tax_values) {
  unique_tax <- sort(unique(tax_values[!is.na(tax_values) & tax_values != "none"]))
  if (length(unique_tax) > 0) {
    tax_colors <- rainbow(length(unique_tax))
    names(tax_colors) <- unique_tax
    return(tax_colors)
  }
  return(NULL)
}

# helper function to generate mge colors based on gene descriptions (vectorized)
get_mge_color <- function(prot_desc_vector, mge_groups, mge_colors) {
  # initialize result vector with default color for non-mge genes
  result <- rep("#E8E8E8", length(prot_desc_vector))
  
  # process each group and its patterns
  for (group_name in names(mge_groups)) {
    patterns <- mge_groups[[group_name]]
    group_color <- mge_colors[[group_name]]
    
    # find genes matching any pattern in this group
    for (pattern in patterns) {
      matches <- grepl(pattern, prot_desc_vector, ignore.case = TRUE)
      # assign color to matches (later matches can override earlier ones)
      result[matches] <- group_color
    }
  }
  
  return(result)
}

# get_genes_f used by gene profile and the gene tab
get_genes_f <- function(assembly) {
  cache(paste0(assembly, "_genes"), {
    genes <- get_data("PRODIGAL_GENE_TABLE", tag = assembly)
    uniref <- get_data("UNIREF_GENE_TAX_TABLE", tag = assembly)
    ix <- match(genes$gene, uniref$gene)
    fields <- c("uniref", "identity", "coverage", "evalue", "bitscore", "prot_desc", "tax", "uniref_count")
    for (field in fields) {
      # set default value based on field type
      if (is.numeric(uniref[[field]])) {
        genes[[field]] <- ifelse(is.na(ix), 0, uniref[[field]][ix])
      } else {
        genes[[field]] <- ifelse(is.na(ix), "none", uniref[[field]][ix])
      }
    }
    
    # add tax color
    tax_color_map <- get_tax_color(genes$tax)
    genes$tax_color <- tax_color_map[genes$tax]
    
    # add mge color
    genes$mge_color <- get_mge_color(genes$prot_desc, mge_groups, mge_colors)
    
    # build label text for tooltips (same logic as default in profile)
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

# define shared get_map_tag function (used by both main_view and variants)
lib.table <- get_data("LONG_LIB_TABLE")
get_map_tag <- function(assembly, timepoint) {
  ix = lib.table$ASSEMBLY_ID == assembly & lib.table$SAMPLE_TYPE == timepoint
  if (sum(ix) == 0 || sum(ix) > 1) {
    return (NULL)
  }
  paste0(assembly, "_", lib.table$LIB_ID[ix])
}

# simple function to get alignment for any assembly/library_id combination
get_aln_f <- function(assembly, library_id) {
  tag <- get_map_tag(assembly, library_id)
  get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
}

# define get_gene_table_f function for variant queries
get_variants_gene_table_f <- function(assembly) {
  # get genes data using the existing get_genes_f function
  genes_data <- get_genes_f(assembly)

  # convert to alntools gene table format
  gene_table <- data.frame(
    gene = genes_data$gene,
    contig = genes_data$contig,
    start = genes_data$start,
    end = genes_data$end,
    strand = genes_data$strand,
    desc = genes_data$prot_desc,
    stringsAsFactors = FALSE
  )
  
  return(gene_table)
}

# load variants table with gene fields
get_variants_table_f <- function(assembly) {
  vars_table <- get_data("POLY_VARIANTS_TABLE", tag = assembly)
  vars_genic <- get_data("POLY_VARIANTS_GENIC", tag = assembly)
  
  if (is.null(vars_table)) {
    return(NULL)
  }
  
  if (is.null(vars_genic)) {
    return(vars_table)
  }
  
  # left join to add gene fields
  merged <- merge(vars_table, vars_genic, 
                  by = "variant_id", 
                  all.x = TRUE, sort = FALSE)
  merged$row_id = merged$variant_id  
  return(merged)
}

# variants tab - static mode (using pre-computed files)
register_tab(
  tab_id = "variants",
  tab_label = "Variants",
  tab_code = "tabs/variants/variants_tab.r",
  is.dynamic = FALSE,  # use pre-computed files
  library_ids = c("early", "pre", "post", "late"),
  get_variants_table_f = get_variants_table_f,
  get_variants_support_f = function(assembly) get_data("POLY_VARIANTS_SUPPORT", tag = assembly),
  get_variants_coverage_f = function(assembly) get_data("POLY_VARIANTS_COVERAGE", tag = assembly),
  use_genes = TRUE,
  supports_export = TRUE
)

# example dynamic mode configuration (commented out)
# register_tab(
#   tab_id = "variants",
#   tab_label = "Variants",
#   tab_code = "tabs/variants/variants_tab.r",
#   is.dynamic = TRUE,  # use alignment queries
#   min_reads = 2,
#   min_coverage = 10,
#   min_libraries = 0,
#   get_aln_f = get_aln_f,
#   library_ids = c("early", "pre", "post", "late"),
#   max_variants = 1000,
#   get_gene_table_f = get_variants_gene_table_f,
#   codon_table_path = "codon_tables/table11",
#   get_fasta_f = get_fasta_f,
#   use_genes = TRUE,
#   supports_export = TRUE
# )

########################################################
# register rearrangements tab
########################################################

# rearrangements tab - static mode (using pre-computed files)
register_tab(
  tab_id = "rearrangements",
  tab_label = "Rearrangements",
  tab_code = "tabs/rearrangements/rearrangements_tab.r",
  is.dynamic = FALSE,  # use pre-computed files
  library_ids = c("early", "pre", "post", "late"),
  get_rearrange_events_f = function(assembly) get_data("POLY_REARRANGE_EVENTS", tag = assembly),
  get_rearrange_support_f = function(assembly) get_data("POLY_REARRANGE_SUPPORT", tag = assembly),
  get_rearrange_coverage_f = function(assembly) get_data("POLY_REARRANGE_COVERAGE", tag = assembly),
  use_genes = TRUE,
  supports_export = TRUE
)

# example dynamic mode configuration (commented out)
# register_tab(
#   tab_id = "rearrangements",
#   tab_label = "Rearrangements",
#   tab_code = "tabs/rearrangements/rearrangements_tab.r",
#   is.dynamic = TRUE,  # use alignment queries
#   get_aln_f = get_aln_f,
#   library_ids = c("early", "pre", "post", "late"),
#   max_margin = 10,
#   min_element_length = 50,
#   min_anchor_length = 200,
#   max_anchor_mutations_percent = 1.0,
#   max_element_mutation_percent = 1.0,
#   supports_export = TRUE
# )

########################################################
# register associations tab
########################################################

# get library IDs from associations library table
get_associations_library_ids <- function(assembly) {
  lib_table <- get_data("BINNING_ASSOCIATIONS_LIBRARY_TABLE", tag = assembly, null.on.missing = TRUE)
  if (is.null(lib_table)) {
    return(c("early", "pre", "post", "late"))  # fallback to default
  }
  # filter by assembly and return LIB_ID column
  lib_table <- lib_table[lib_table$ASSEMBLY_ID == assembly, ]
  if (nrow(lib_table) == 0) {
    return(c("early", "pre", "post", "late"))  # fallback to default
  }
  return(lib_table$LIB_ID)
}

# associations tab
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

# poly tab - unified locals, rearrangements, and elements
register_tab(
  tab_id = "poly",
  tab_label = "Poly",
  tab_code = "tabs/poly/poly_tab.r",
  library_ids = c("early", "pre", "post", "late"),
  get_unify_table_f = function(assembly) get_data("EVO_UNIFY_TABLE_ASSEMBLY", tag = assembly, null.on.missing = TRUE),
  get_unify_support_f = function(assembly) get_data("EVO_UNIFY_SUPPORT_ASSEMBLY", tag = assembly, null.on.missing = TRUE),
  get_unify_coverage_f = function(assembly) get_data("EVO_UNIFY_COVERAGE_ASSEMBLY", tag = assembly, null.on.missing = TRUE),
  get_abundance_f = function(assembly) get_data("BINNING_ABUNDANCE_LR_MAT", tag = assembly, null.on.missing = TRUE)
)

########################################################
# bin segments profile function
########################################################

# function to get segments with bin colors for segments profile
get_bin_segments_f <- function(assembly) {
  seg_table <- get_seg_bins(assembly)  # now includes bin_color
  if (is.null(seg_table)) {
    return(NULL)
  }
  
  # format for interval_profile (needs: assembly, contig, start, end, desc, id)
  seg_table$assembly <- assembly
  seg_table$desc <- paste0("Segment: ", seg_table$segment, "\nBin: ", seg_table$bin)
  seg_table$id <- seg_table$segment
  
  return(seg_table)
}
