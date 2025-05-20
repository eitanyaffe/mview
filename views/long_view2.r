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
# Register profiles
########################################################

# Register a test points profile with minimal data
points_profile(
  id = "contig_length",
  name = "Contig Length",
  data = get_contigs_f,
  params = list(color = list(type = "string", default = "red"))
)

axis_profile()
