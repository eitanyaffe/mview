########################################################
# Set regions directory
########################################################

set_regions_dir("configs/c60/regions")

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
  data.frame(contig = df$contig, length = df$length, coverage = df$coverage, circular = df$circular)
})

# Register genomes function
register_genomes_f(function(assembly = NULL) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
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

# Register fasta function
register_fasta_f(function(assembly = NULL) {
  get_data("ASSEMBLY_CONTIG_FILE", tag = assembly, read_f = function(path) {
    # check if seqinr is available
    if (!requireNamespace("seqinr", quietly = TRUE)) {
      stop("seqinr package not available, install with: install.packages('seqinr')")
    }
    seqinr::read.fasta(file = path, seqtype = "DNA", as.string = TRUE)
  })
})

########################################################
# register views
########################################################

view_file <- "configs/c60/c60_main_view.r"
view_register("early", view_file, timepoints = c("early"))
view_register("pre", view_file, timepoints = c("pre"))
view_register("post", view_file, timepoints = c("post"))
view_register("late", view_file, timepoints = c("late"))
view_register("pre / post", view_file, timepoints = c("pre", "post"))
view_register("pre / post / late", view_file, timepoints = c("pre", "post", "late"))
view_register("all", view_file, timepoints = c("early", "pre", "post", "late"))

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