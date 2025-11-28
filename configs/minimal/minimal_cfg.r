########################################################
# Minimal test configuration file
########################################################

# Set regions directory
set_regions_dir("configs/minimal/regions")

# example data directories
example_dir <- "examples"
tables_dir <- file.path(example_dir, "tables")

# helper function that caches tables
read_cached <- function(key, path, read_f = read.delim) {
  cache(key, {
    if (identical(read_f, read.delim)) {
      read_f(path, stringsAsFactors = FALSE)
    } else {
      read_f(path)
    }
  })
}

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
# get functions
########################################################

get_assembly_f <- function() {
  read_cached("assembly_table", file.path(tables_dir, "assembly_table.txt"))
}

get_contigs_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("contig_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("contigs_", assembly), path)
  data.frame(contig = df$contig, length = df$length, circular = df$circular)
}

get_segments_f <- function(assembly = NULL) {
  contigs <- get_contigs_f(assembly)
  if (is.null(contigs)) return(NULL)
  data.frame(
    segment = paste0("s", seq_len(nrow(contigs))),
    contig = contigs$contig,
    start = 1L,
    end = contigs$length,
    length = contigs$length,
    stringsAsFactors = FALSE
  )
}

get_genomes_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("genome_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("genomes_", assembly), path)
  data.frame(gid = df$gid, length = df$length)
}

get_segment_map_f <- function(assembly = NULL) {
  contigs <- get_contigs_f(assembly)
  if (is.null(contigs)) return(NULL)
  data.frame(
    segment = paste0("s", seq_len(nrow(contigs))),
    gid = contigs$contig,
    stringsAsFactors = FALSE
  )
}

get_aln_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path("examples", "aln", paste0(assembly, ".aln"))
  if (!file.exists(path)) return(NULL)
  read_cached(paste0("aln_", assembly), path, read_f = aln_load)
}

get_fasta_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("contig_fasta_", assembly, ".fasta"))
  if (!file.exists(path)) return(NULL)
  
  # check if seqinr is available
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("seqinr package not available, install with: install.packages('seqinr')")
  }
  
  read_cached(paste0("fasta_", assembly), path, read_f = function(file_path) {
    seqinr::read.fasta(file = file_path, seqtype = "DNA", as.string = TRUE)
  })
}

# helper function to generate taxonomy colors
get_tax_color <- function(tax_values) {
  unique_tax <- sort(unique(tax_values[!is.na(tax_values) & tax_values != ""]))
  if (length(unique_tax) > 0) {
    tax_colors <- rainbow(length(unique_tax))
    names(tax_colors) <- unique_tax
    return(tax_colors)
  }
  return(NULL)
}

# helper function to generate mge colors based on gene descriptions (vectorized)
get_mge_color <- function(prot_desc_vector, mge_groups, mge_colors) {
  # initialize result vector
  result <- rep(NA, length(prot_desc_vector))
  
  # process each group and its patterns
  for (group_name in names(mge_groups)) {
    patterns <- mge_groups[[group_name]]
    group_color <- mge_colors[[group_name]]
    
    # find genes matching any pattern in this group
    for (pattern in patterns) {
      matches <- grepl(pattern, prot_desc_vector, ignore.case = TRUE)
      # only assign color if not already assigned (first match wins)
      result[matches & is.na(result)] <- group_color
    }
  }
  
  return(result)
}

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

get_genes_f <- function(assembly) {
  gpath <- file.path(tables_dir, "gene_table.txt")
  if (!file.exists(gpath)) return(NULL)
  genes <- read_cached("gene_table", gpath)
  genes <- genes[genes$assembly == assembly, ]
  # select fields
  fields <- c("gene","contig","start","end","strand","uniref","identity","coverage","evalue","bitscore","prot_desc","tax","uniref_count")
  genes <- genes[, intersect(names(genes), fields)]
  
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
}

########################################################
# register assemblies, genomes and contigs
########################################################

atbl <- get_assembly_f()
aids <- atbl$assembly_id
set_assemblies(aids)

# register data functions
register_contigs_f(get_contigs_f)
register_segments_f(get_segments_f)
register_genomes_f(get_genomes_f)
register_segment_map_f(get_segment_map_f)
register_fasta_f(get_fasta_f)

########################################################
# register views
########################################################

view_file <- "configs/minimal/minimal_view.r"
view_register("example", view_file)

########################################################
# register gene tab
########################################################

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

