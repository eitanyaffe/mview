# tRNA gene profile
#
# Reads TRNA_TABLE (per assembly) via get_trna_f.
# Columns: trna_id, contig, start, end, strand, isotype, anticodon,
#          intron_start, intron_end, score
# Colors one track segment per tRNA; color keyed by isotype.

default_trna_colors <- local({
  standard <- c("Ala", "Arg", "Asn", "Asp", "Cys",
                "Gln", "Glu", "Gly", "His", "Ile", "Ile2",
                "Leu", "Lys", "Met", "fMet", "Phe", "Pro",
                "SeC", "Ser", "Thr", "Trp", "Tyr", "Val")
  cols <- rainbow(length(standard), s = 0.75, v = 0.85)
  names(cols) <- standard
  c(cols, Sup = "#AAAAAA", Undet = "#E8E8E8")
})

trna_profile <- function(id, name, get_trna_f, colors = default_trna_colors,
                          height = 30, is_fixed = TRUE,
                          auto_register = TRUE) {

  trna_f <- function(assembly) {
    trna <- cache(paste0(assembly, "_trna_", id), {
      get_trna_f()
    })

    if (is.null(trna) || nrow(trna) == 0)
      return(NULL)

    # gene_profile requires a 'gene' column for compatibility
    trna$gene <- trna$trna_id

    # color by isotype; unknown isotypes fall back to gray
    trna$trna_color <- ifelse(trna$isotype %in% names(colors),
                               colors[trna$isotype], "#E8E8E8")

    trna$label <- paste0(
      "tRNA: ",      trna$trna_id,   "\n",
      "Isotype: ",   trna$isotype,   "\n",
      "Anticodon: ", trna$anticodon, "\n",
      "Strand: ",    trna$strand,    "\n",
      "Score: ",     trna$score
    )

    trna
  }

  gene_profile(
    id            = id,
    name          = name,
    height        = height,
    is_fixed      = is_fixed,
    gene_f        = trna_f,
    color_field   = "trna",
    label_field   = "label",
    params        = NULL,
    auto_register = auto_register
  )
}
