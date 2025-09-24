# tabs/variants/variants_utils.r
# standalone variant processing utilities

# query variants for a specific context (extracted from variants_tab.r)
query_variants_for_context <- function(assembly, contigs, zoom, tab_config) {
  # early return if no context
  if (is.null(assembly) || length(contigs) == 0) {
    return(NULL)
  }
  
  # extract required parameters from tab config
  min_reads <- tab_config$min_reads
  min_coverage <- tab_config$min_coverage
  min_libraries <- tab_config$min_libraries
  get_aln_f <- tab_config$get_aln_f
  library_ids <- tab_config$library_ids
  use_genes <- tab_config$use_genes
  get_gene_table_f <- tab_config$get_gene_table_f
  get_fasta_f <- tab_config$get_fasta_f
  codon_table_path <- tab_config$codon_table_path
  
  # get alignment stores from all configured library_ids
  stores <- list()
  
  for (library_id in library_ids) {
    tryCatch({
      aln <- get_aln_f(assembly, library_id)
      if (!is.null(aln) && inherits(aln, "externalptr")) {
        stores[[library_id]] <- aln
      } else {
        warning(sprintf("get_aln_f returned invalid alignment for library_id '%s'", library_id))
      }
    }, error = function(e) {
      warning(sprintf("error getting alignment for library_id '%s': %s", library_id, e$message))
    })
  }
  
  if (length(stores) == 0) {
    return(NULL)
  }
  
  # get contigs data and build context for intervals
  contigs_table <- get_contigs(assembly)
  if (is.null(contigs_table)) {
    return(NULL)
  }
  
  # build context for intervals
  cxt <- build_context(contigs, contigs_table, zoom, assembly)
  if (is.null(cxt) || is.null(cxt$intervals) || nrow(cxt$intervals) == 0) {
    return(NULL)
  }
  
  # get alignment filter parameters (use defaults if not available)
  get_filter_param <- function(param_id, default_value, type_converter = identity) {
    param_key <- paste0("align_filter_", param_id)
    if (!param_key %in% names(param_registry)) {
      return(default_value)
    }
    
    param_accessor <- get_param("align_filter", param_id)
    value <- param_accessor()
    return(type_converter(value))
  }
  
  clip_mode <- get_filter_param("clip_mode", "all")
  clip_margin <- get_filter_param("clip_margin", 10L, as.integer)
  min_mutations_percent <- get_filter_param("min_mutations_percent", 0.0, as.numeric)
  max_mutations_percent <- get_filter_param("max_mutations_percent", 10.0, as.numeric)
  min_alignment_length <- get_filter_param("min_alignment_length", 0L, as.integer)
  max_alignment_length <- get_filter_param("max_alignment_length", 0L, as.integer)
  min_indel_length <- get_filter_param("min_indel_length", 3L, as.integer)
  
  # create genes object for variant annotation
  genes_object <- NULL
  if (use_genes) {
    if (is.function(get_gene_table_f) && is.function(get_fasta_f)) {
      # create genes object with caching
      genes_object <- cache(paste0(assembly, "_genes_object"), {
        # get gene table data directly
        gene_table_data <- get_gene_table_f(assembly)
        if (!is.null(gene_table_data) && nrow(gene_table_data) > 0) {
          # get reference sequences directly
          reference_sequences_data <- get_fasta_f(assembly)
          
          genes_obj <- aln_genes(
            gene_table = gene_table_data,
            reference_sequences = reference_sequences_data,
            codon_table_path = codon_table_path
          )
          
          cat(sprintf("using gene annotation: %d genes loaded, %d reference sequences\n", 
                     nrow(gene_table_data), length(reference_sequences_data)))
          
          genes_obj
        } else {
          cat("gene table function returned no data, proceeding without gene annotation\n")
          NULL
        }
      })
    }
  }
  
  # call aln_query_variants with alignment filter parameters
  variants_results <- aln_query_variants(
    store_list = stores,
    intervals_df = cxt$intervals,
    min_variants_variant_support = min_reads,
    min_variants_library_support = min_libraries,
    min_variants_coverage_support = min_coverage,
    clip_mode_str = clip_mode,
    clip_margin = clip_margin,
    min_mutations_percent = min_mutations_percent,
    max_mutations_percent = max_mutations_percent,
    min_alignment_length = min_alignment_length,
    max_alignment_length = max_alignment_length,
    min_indel_length = min_indel_length,
    genes = genes_object
  )
  
  # enhance variants with gene annotation if genes were used
  if (use_genes) {
    variants_df <- variants_results$variants
    
    # only add gene annotation columns if there are variants
    if (nrow(variants_df) > 0) {
      # check if gene annotation data is available
      if (!is.null(variants_results$genic) && nrow(variants_results$genic) > 0) {
        genic_data <- variants_results$genic
        # match variant_id with row_id in genic data
        ix <- match(variants_df$variant_id, genic_data$row_id)
        
        # add gene columns using vectorized operations
        variants_df$gene_desc <- ifelse(!is.na(ix), genic_data$gene_desc[ix], "none")
        variants_df$mutation_desc <- ifelse(!is.na(ix), genic_data$mutation_desc[ix], "")
        
        cat(sprintf("enhanced %d variants with gene annotation data\n", sum(!is.na(ix))))
      } else {
        # no gene data available, add default columns
        variants_df$gene_desc <- "none"
        variants_df$mutation_desc <- ""
        cat("no gene annotation data available\n")
      }
      
      # update the variants results
      variants_results$variants <- variants_df
    } else {
      cat("no variants found, skipping gene annotation\n")
    }
  }
  
  return(variants_results)
}

# filter variants by span (frequency range) - extracted from variants_tab.r
filter_variants_by_span <- function(variant_data, min_span) {
  if (is.null(variant_data) || is.null(variant_data$support) || is.null(variant_data$coverage)) {
    return(variant_data)
  }
  
  # calculate frequency for each variant across libraries
  support_matrix <- variant_data$support
  coverage_matrix <- variant_data$coverage
  
  # avoid division by zero
  freq_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
  
  # calculate span (max - min frequency) for each variant
  variant_spans <- apply(freq_matrix, 1, function(row) {
    valid_freqs <- row[!is.na(row) & is.finite(row)]
    if (length(valid_freqs) == 0) return(0)
    max(valid_freqs) - min(valid_freqs)
  })
  
  # filter variants meeting span threshold
  keep_variants <- variant_spans >= min_span
  
  # filter all components
  filtered_data <- list(
    variants = variant_data$variants[keep_variants, ],
    support = variant_data$support[keep_variants, , drop = FALSE],
    coverage = variant_data$coverage[keep_variants, , drop = FALSE],
    library_ids = variant_data$library_ids
  )
  
  return(filtered_data)
}

# add color information to variants dataframe - extracted from variants_tab.r
add_variant_colors <- function(variants_df) {
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    return(variants_df)
  }
  
  # load alignment color functions if not available
  if (!exists("get_variant_type_colors") || !exists("get_mutation_type_colors")) {
    source("profiles/align/align_utils.r")
  }
  
  # get mutation color mode from alignment profiles
  mutation_color_mode <- profile_get_param("mutation_color_mode", "align", "detailed")
  
  # apply coloring based on mutation mode and add colors to variants dataframe
  if (mutation_color_mode == "detailed") {
    unique_descs <- unique(variants_df$desc)
    type_colors <- get_variant_type_colors(unique_descs)
    names(type_colors) <- unique_descs
    
    missing_descs <- unique_descs[!unique_descs %in% names(type_colors)]
    if (length(missing_descs) > 0) {
      black_colors <- rep("#000000", length(missing_descs))
      names(black_colors) <- missing_descs
      type_colors <- c(type_colors, black_colors)
    }
    
    # add color column to variants dataframe
    variants_df$color <- type_colors[variants_df$desc]
  } else {
    type_mapping <- c("sub" = "substitution", "ins" = "addition", "del" = "deletion")
    color_type <- type_mapping[variants_df$type]
    color_type[is.na(color_type)] <- "clip"
    
    color_types <- unique(color_type)
    type_colors <- get_mutation_type_colors(color_types)
    names(type_colors) <- color_types
    type_colors["clip"] <- "#000000"
    
    # add color column to variants dataframe
    variants_df$color <- type_colors[color_type]
  }
  
  return(variants_df)
}

