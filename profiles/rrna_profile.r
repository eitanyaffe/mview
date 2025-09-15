# default rRNA colors
default_rrna_colors <- c("16S" = "#7BB3D3", "23S" = "#C97BA8", "5S" = "#F8B563")

# helper function to read GFF files
read_gff_f <- function(path) {
  # read GFF file, skipping comment lines, no header
  df <- read.delim(path, sep = "\t", header = FALSE, comment.char = "#", 
                   stringsAsFactors = FALSE)
  # assign standard GFF column names
  colnames(df) <- c("contig", "source", "type", "start", "end", 
                    "score", "strand", "phase", "attributes")
  
  # parse Name from attributes (Name=16S_rRNA;product=...)
  df$name <- gsub(".*Name=([^;]+).*", "\\1", df$attributes)
  df$product <- gsub(".*product=([^;]*)", "\\1", df$attributes)
  
  # extract rRNA type for coloring
  df$rrna_type <- gsub("_rRNA", "", df$name)
  
  # create gene column for compatibility
  df$gene <- paste0(df$contig, "_", df$start, "_", df$end)
  
  return(df)
}

rrna_profile <- function(id, name, get_gff_f, colors = default_rrna_colors,
                         height = 30, is_fixed = TRUE,
                         auto_register = TRUE) {
  
  # create data function that handles caching and color/label assignment
  rrna_f <- function(cxt) {
    rrna <- cache(paste0(cxt$assembly, "_rrna_", id), {
      # get_gff_f should call get_data with read_f=read_gff_f
      get_gff_f(cxt)
    })
    
    if (is.null(rrna) || nrow(rrna) == 0) {
      return(NULL)
    }
    
    # assign colors and labels after cache (so they update when definitions change)
    rrna$rrna_color <- colors[rrna$rrna_type]
    
    # check for missing colors and assign default
    missing_colors <- is.na(rrna$rrna_color)
    if (any(missing_colors)) {
      rrna$rrna_color[missing_colors] <- "#E8E8E8"
    }
    
    # build label text for tooltips
    rrna$label <- paste0(
      "rRNA: ", rrna$name, "\n",
      "Type: ", rrna$rrna_type, "\n",
      "Strand: ", rrna$strand, "\n",
      "Product: ", rrna$product
    )
    
    return(rrna)
  }
  
  # use gene_profile with rRNA-specific settings
  gene_profile(
    id = id,
    name = name,
    height = height, is_fixed = is_fixed,
    gene_f = rrna_f,
    color_field = "rrna",
    label_field = "label",
    params = NULL,
    auto_register = auto_register
  )
}
