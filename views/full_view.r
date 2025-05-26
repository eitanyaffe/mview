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

get_map_tag <- function(assembly) {
  switch(assembly,
    "BAM" = "L1",
    "EBP" = "L2",
    "EBI" = "L9",
    "EBC" = "L11"
  )
}

align_profile(
  id = "align",
  name = "Align",
  aln_f = function(cxt) {
    tag <- paste0(cxt$assembly, "_", get_map_tag(cxt$assembly))
    get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
  },
  params = list(height_style = list(
    type = "select",
    choices = c("by_mutations", "by_coord"), default = "by_mutations"
  ))
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
