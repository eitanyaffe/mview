########################################################
# Minimal test configuration file
########################################################

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

get_genomes_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("genome_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("genomes_", assembly), path)
  data.frame(gid = df$gid, length = df$length)
}

get_contig_map_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("contig_map_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("contigmap_", assembly), path)
  data.frame(contig = df$contig, gid = df$gid)
}

get_aln_f <- function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path("examples", "aln", paste0(assembly, ".aln"))
  if (!file.exists(path)) return(NULL)
  read_cached(paste0("aln_", assembly), path, read_f = aln_load)
}

get_genes_f <- function(cxt) {
  gpath <- file.path(tables_dir, "gene_table.txt")
  if (!file.exists(gpath)) return(NULL)
  genes <- read_cached("gene_table", gpath)
  genes <- genes[genes$assembly == cxt$assembly, ]
  # select fields
  fields <- c("gene","contig","start","end","strand","uniref","identity","coverage","evalue","bitscore","prot_desc","tax","uniref_count")
  genes[, intersect(names(genes), fields)]
}

########################################################
# register assemblies, genomes and contigs
########################################################

atbl <- get_assembly_f()
aids <- atbl$assembly_id
set_assemblies(aids)

# register data functions
register_contigs_f(get_contigs_f)
register_genomes_f(get_genomes_f)
register_contig_map_f(get_contig_map_f)

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

