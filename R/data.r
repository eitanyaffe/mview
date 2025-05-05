contigs <- data.frame(
  cid = paste0("c", 1:10),
  length = sample(100:500, 10),
  stringsAsFactors = FALSE
)

genomes <- data.frame(
  gid = paste0("g", 1:3),
  length = sample(1000:2000, 3),
  stringsAsFactors = FALSE
)

contig_map <- data.frame(
  cid = c("c1", "c2", "c3", "c4", "c5", "c6", "c1", "c7"),
  gid = c("g1", "g1", "g1", "g2", "g2", "g3", "g2", "g3"),
  stringsAsFactors = FALSE
)
