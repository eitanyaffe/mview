########################################################
# Register lookup files
########################################################

# project_name <- "pb4"
project_name <- "pb-b20"
# project_name <- "pb-pilot2"

# search for lookup files under MAKESHIFT_ROOT/export/long/pb4/
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
  data.frame(contig = df$contig, length = df$length, circular = df$circular)
})

# Register genomes function
register_genomes_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_GENOME_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  df$gid <- df$contig
  data.frame(gid = df$gid, length = df$length)
})

# Register contig map function
register_contig_map_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(contig = df$contig, gid = df$contig)
})

########################################################
# register views
########################################################

view_file <- "configs/c60/c60_main_view.r"
view_register("early", view_file, timepoints = c("early"))
view_register("pre", view_file, timepoints = c("pre"))
view_register("post", view_file, timepoints = c("post"))
view_register("late", view_file, timepoints = c("late"))
view_register("pre vs. post", view_file, timepoints = c("pre", "post"))
view_register("all", view_file, timepoints = c("early", "pre", "post", "late"))

########################################################
# register gene tab
########################################################

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