#!/usr/bin/env Rscript

# Script to create a simplified test configuration based on c60 config
# Creates example data with 2 assemblies, 10 shortest contigs each
# Includes genes only (no alignments), plus axis profile
# Creates 2 genomes per assembly (3 contigs + 7 contigs each)

library(dplyr)

cat("creating test configuration from c60...\n")

# source required modules
source("core/cache.r")
source("core/data.r")

# initialize cache system
cache_init("minimal_example")

# load only the data setup from c60 configuration
# extract the project name and setup data access
project_name <- "pb-b20"
dir <- paste(Sys.getenv("MAKESHIFT_ROOT"), "export", "long", project_name, "default", sep = "/")
fns <- list.files(dir, full.names = TRUE, pattern = "*.txt")
set_lookup(fns)

# create output directories
dir.create("examples/tables", recursive = TRUE, showWarnings = FALSE)

########################################################
# get original data and select assemblies
########################################################

# get assembly table and select first 2 assemblies
assembly_table <- get_data("ASSEMBLY_TABLE")
original_assemblies <- assembly_table$ASSEMBLY_ID[1:2]
cat(sprintf("selected assemblies: %s\n", paste(original_assemblies, collapse = ", ")))

########################################################
# process each assembly
########################################################

all_contig_data <- list()
all_gene_rows <- list()
all_genome_data <- list()
all_contig_map_data <- list()

gene_cols <- c("assembly", "gene", "contig", "start", "end", "strand")

for (assembly in original_assemblies) {
  cat(sprintf("processing assembly %s\n", assembly))
  
  # get contig data and select contigs above 10kb, then take 10 shortest of those
  contig_data <- get_data("ASSEMBLY_CONTIG_TABLE", tag = assembly)
  contig_data_filtered <- contig_data[contig_data$length > 10000, ]  # above 10kb
  contig_data_filtered <- contig_data_filtered[order(contig_data_filtered$length), ]  # shortest first
  selected_contigs <- head(contig_data_filtered, 10)
  
  # extract fasta sequences for selected contigs
  # get fasta file path (not content) by using a custom read function that returns the path
  contig_fasta_path <- get_data("ASSEMBLY_CONTIG_FILE", tag = assembly, read_f = function(path) path)
  if (is.null(contig_fasta_path) || !file.exists(contig_fasta_path)) {
    stop(sprintf("fasta file not found for assembly %s", assembly))
  }
  
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("seqinr package required for fasta processing. install with: install.packages('seqinr')")
  }
  
  # read original fasta
  original_sequences <- seqinr::read.fasta(file = contig_fasta_path, seqtype = "DNA", as.string = TRUE)
  
  # extract sequences for selected contigs
  selected_contigs_list <- selected_contigs$contig
  extracted_sequences <- list()
  for (contig in selected_contigs_list) {
    # find matching sequence
    seq_found <- FALSE
    for (seq_name in names(original_sequences)) {
      if (seq_name == contig || grepl(paste0("^", contig, "($|\\s)"), seq_name)) {
        extracted_sequences[[contig]] <- original_sequences[[seq_name]]
        seq_found <- TRUE
        break
      }
    }
    if (!seq_found) {
      stop(sprintf("sequence not found for contig %s", contig))
    }
  }
  
  # write extracted fasta
  output_fasta <- paste0("examples/tables/contig_fasta_", assembly, ".fasta")
  seqinr::write.fasta(extracted_sequences, names(extracted_sequences), output_fasta)
  cat(sprintf("  wrote %d sequences to %s\n", length(extracted_sequences), output_fasta))
  
  # keep only minimal, lowercase columns
  selected_contigs <- selected_contigs[, intersect(names(selected_contigs), c("contig", "length", "circular", "coverage"))]
  names(selected_contigs) <- tolower(names(selected_contigs))
  if (!"circular" %in% names(selected_contigs)) {
    selected_contigs$circular <- FALSE
  }
  selected_contigs <- selected_contigs[, c("contig", "length", "circular")]
  all_contig_data[[assembly]] <- selected_contigs
  
  # get gene data for selected contigs (single combined table later)
  gene_data <- get_data("PRODIGAL_GENE_TABLE", tag = assembly)
  original_contigs <- selected_contigs$contig
  gene_data_filtered <- gene_data[gene_data$contig %in% original_contigs, ]
  
  # keep minimal lowercase columns and add assembly
  gene_data_filtered <- gene_data_filtered[, intersect(names(gene_data_filtered), c("gene","contig","start","end","strand"))]
  names(gene_data_filtered) <- tolower(names(gene_data_filtered))
  
  # add annotation from uniref table (keep original case for now, lowercase column names later)
  uniref_data <- get_data("UNIREF_GENE_TAX_TABLE", tag = assembly)
  uniref_keep <- c("gene","uniref","identity","coverage","evalue","bitscore","prot_desc","tax","uniref_count")
  uniref_data <- uniref_data[, intersect(names(uniref_data), uniref_keep)]
  if (nrow(uniref_data) > 0) {
    ix <- match(gene_data_filtered$gene, uniref_data$gene)
    add_fields <- setdiff(names(uniref_data), "gene")
    for (fld in add_fields) {
      gene_data_filtered[[tolower(fld)]] <- ifelse(is.na(ix), NA, uniref_data[[fld]][ix])
    }
  }
  
  gene_data_filtered$assembly <- assembly
  gene_data_filtered <- gene_data_filtered[, unique(c(gene_cols, setdiff(names(gene_data_filtered), gene_cols)))]
  all_gene_rows[[assembly]] <- gene_data_filtered
  
  # create genomes (2 per assembly)
  genome1_contigs <- selected_contigs$contig[1:3]
  genome2_contigs <- selected_contigs$contig[1:7]
  genome1_length <- sum(selected_contigs$length[1:3])
  genome2_length <- sum(selected_contigs$length[1:7])
  genome_data <- data.frame(
    gid = c(paste0(assembly, "_genome1"), paste0(assembly, "_genome2")),
    length = c(genome1_length, genome2_length),
    stringsAsFactors = FALSE
  )
  all_genome_data[[assembly]] <- genome_data
  
  # contig map
  contig_map_data <- data.frame(
    contig = c(genome1_contigs, genome2_contigs),
    gid = c(rep(paste0(assembly, "_genome1"), 3), rep(paste0(assembly, "_genome2"), 7)),
    stringsAsFactors = FALSE
  )
  all_contig_map_data[[assembly]] <- contig_map_data
  
  cat(sprintf("  selected %d contigs, %d genes\n", 
              nrow(selected_contigs), nrow(gene_data_filtered)))
}

########################################################
# save tables (lowercase, minimal fields) and single gene table
########################################################

# assemblies (lowercase column, minimal)
final_assembly_table <- data.frame(
  assembly_id = original_assemblies,
  stringsAsFactors = FALSE
)

# single combined gene table
final_gene_table <- do.call(rbind, all_gene_rows)

# compute totals
total_genes <- nrow(final_gene_table)

# write files
write.table(final_assembly_table, "examples/tables/assembly_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

for (assembly in original_assemblies) {
  write.table(all_contig_data[[assembly]], paste0("examples/tables/contig_table_", assembly, ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(all_genome_data[[assembly]], paste0("examples/tables/genome_table_", assembly, ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(all_contig_map_data[[assembly]], paste0("examples/tables/contig_map_table_", assembly, ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
}

write.table(final_gene_table, "examples/tables/gene_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)