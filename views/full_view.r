########################################################
# helper functions
########################################################

get_contigs_f <- function(cxt) {
  df <- get_data("ASSEMBLY_CONTIG_TABLE", tag = cxt$assembly)
  if (is.null(df)) {
    return(NULL)
  }
  data.frame(contig = df$contig, coord = df$length / 2, value = df$length, desc = paste0(df$contig, " (", format_bp(df$length), ")"))
}

########################################################
# Align code
########################################################

library(Rcpp)

alntools_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/alntools/cpp")
cdir <- ".Rcpp_dir"
rebuild <- F
sourceCpp(paste0(alntools_dir, "/aln_R.cpp"), verbose = T, cacheDir = cdir, rebuild = rebuild)

########################################################
# Register profiles
########################################################

align_profile(
  id = "aln_L1",
  name = "Align L1",
  aln_f = make_get_data_f(
    id = "MINIMAP_LIB_ALN", tag = "BAM_L1",
    read_f = aln_load
  )
)

gene_profile(
  id = "genes",
  name = "Genes",
  gene_f = function(cxt) {
    cache(paste0(cxt$assembly, "_genes"), {
      genes <- get_data("PRODIGAL_GENE_TABLE", tag = cxt$assembly)
      uniref <- get_data("UNIREF_GENE_TAX_TABLE", tag = cxt$assembly)
      ix <- match(genes$gene, uniref$gene)
      fields <- c("uniref", "identity", "coverage", "evalue", "bitscore", "prot_desc", "tax", "uniref_count")
      for (field in fields) {
        # Set default value based on field type
        if (is.numeric(uniref[[field]])) {
          genes[[field]] <- ifelse(is.na(ix), 0, uniref[[field]][ix])
        } else {
          genes[[field]] <- ifelse(is.na(ix), "none", uniref[[field]][ix])
        }
      }
      genes
    })
  },
)

axis_profile()
