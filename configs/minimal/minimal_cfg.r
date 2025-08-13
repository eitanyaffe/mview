########################################################
# Example configuration for learning mview (no-lookup)
########################################################

# example data directories
example_dir <- "examples"
tables_dir <- file.path(example_dir, "tables")

# helper cached readers
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
# register assemblies, genomes and contigs
########################################################

atbl <- read_cached("assembly_table", file.path(tables_dir, "assembly_table.txt"))
aids <- atbl$assembly_id
set_assemblies(aids)

# register contigs function
register_contigs_f(function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("contig_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("contigs_", assembly), path)
  data.frame(contig = df$contig, length = df$length, circular = df$circular)
})

# register genomes function  
register_genomes_f(function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("genome_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("genomes_", assembly), path)
  data.frame(gid = df$gid, length = df$length)
})

# register contig map function
register_contig_map_f(function(assembly = NULL) {
  if (is.null(assembly)) return(NULL)
  path <- file.path(tables_dir, paste0("contig_map_table_", assembly, ".txt"))
  if (!file.exists(path)) return(NULL)
  df <- read_cached(paste0("contigmap_", assembly), path)
  data.frame(contig = df$contig, gid = df$gid)
})

########################################################
# register views
########################################################

view_file <- "configs/minimal/minimal_view.r"
view_register("example", view_file)

########################################################
# register gene tab
########################################################

# get_genes_f used by gene profile and the gene tab
get_genes_f <- function(cxt) {
  cache(paste0(cxt$assembly, "_genes"), {
    gpath <- file.path(tables_dir, "gene_table.txt")
    if (!file.exists(gpath)) return(NULL)
    genes <- read_cached("gene_table", gpath)
    genes <- genes[genes$assembly == cxt$assembly, ]
    # select fields
    keep <- c("gene","contig","start","end","strand","uniref","identity","coverage","evalue","bitscore","prot_desc","tax","uniref_count")
    genes[, intersect(names(genes), keep)]
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

