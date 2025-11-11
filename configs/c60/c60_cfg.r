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
# Register assemblies, genomes and contigs
########################################################

aids <- get_data("ASSEMBLY_TABLE")$ASSEMBLY_ID

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

# Register genomes function
register_genomes_f(function(assembly = NULL) {
    # df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
    df <- get_data("BINNING_HOST_TABLE", tag = assembly)
    
    if (is.null(df)) {
        return(NULL)
    }
    # df$gid <- df$contig
    rr = data.frame(gid = df$bin, df[,-match(c("bin","domain"), names(df))])
    rr = rr[order(-rr$length),]
})

# Register contig map function
register_contig_map_f(function(assembly = NULL) {
    # df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
    df <- get_data("BINNING_BIN_SEGMENT_TABLE", tag = assembly)
    if (is.null(df)) {
        return(NULL)
    }
    # data.frame(contig = df$contig, gid = df$contig)
    rr = data.frame(contig = df$contig, gid = df$bin)
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
get_genes_f <- function(cxt) {
  cache(paste0(cxt$assembly, "_genes"), {
    genes <- get_data("PRODIGAL_GENE_TABLE", tag = cxt$assembly)
    uniref <- get_data("UNIREF_GENE_TAX_TABLE", tag = cxt$assembly)
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
  cxt <- list(assembly = assembly)
  genes_data <- get_genes_f(cxt)

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
  vars_table <- get_data("MINIMAP_VARS_TABLE", tag = assembly)
  vars_genic <- get_data("MINIMAP_VARS_GENIC", tag = assembly)
  
  if (is.null(vars_table)) {
    return(NULL)
  }
  
  if (is.null(vars_genic)) {
    return(vars_table)
  }
  
  # left join to add gene fields
  merged <- merge(vars_table, vars_genic, 
                  by.x = "variant_id", by.y = "row_id", 
                  all.x = TRUE, sort = FALSE)
  
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
  get_variants_support_f = function(assembly) get_data("MINIMAP_VARS_SUPPORT", tag = assembly),
  get_variants_coverage_f = function(assembly) get_data("MINIMAP_VARS_COVERAGE", tag = assembly),
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
  get_rearrange_events_f = function(assembly) get_data("MINIMAP_REARRANGE_EVENTS", tag = assembly),
  get_rearrange_support_f = function(assembly) get_data("MINIMAP_REARRANGE_SUPPORT", tag = assembly),
  get_rearrange_coverage_f = function(assembly) get_data("MINIMAP_REARRANGE_COVERAGE", tag = assembly),
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
