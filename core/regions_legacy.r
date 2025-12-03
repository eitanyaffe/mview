# core/regions_legacy.r
# Legacy region file loader - converts old format to new format
# This file can be removed once all region files are migrated

load_legacy_region_table <- function(file_path, regions_dir) {
  cat("[LEGACY] loading old format file:", basename(file_path), "\n")
  
  df <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Check if this is actually old format
  if (!"contigs" %in% colnames(df)) {
    return(NULL)  # Not old format
  }
  
  cat("[LEGACY] detected old format, converting...\n")
  
  # Required columns for old format
  required_cols <- c("id", "description", "assembly", "contigs", "zoom_start", "zoom_end")
  if (!all(required_cols %in% colnames(df))) {
    stop(paste("legacy file missing required columns:", paste(setdiff(required_cols, colnames(df)), collapse = ", ")))
  }
  
  # Create new format data frame
  new_df <- data.frame(
    id = df$id,
    level = if ("level" %in% colnames(df)) df$level else rep(1L, nrow(df)),
    description = df$description,
    assembly = df$assembly,
    segments = character(nrow(df)),
    xlim_start = df$zoom_start,
    xlim_end = df$zoom_end,
    stringsAsFactors = FALSE
  )
  
  # Convert each row: contigs + zoom -> segment IDs
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    
    # Parse contigs
    contigs_list <- if (row$contigs == "" || is.na(row$contigs)) {
      character(0)
    } else {
      trimws(strsplit(row$contigs, ",")[[1]])
    }
    
    if (length(contigs_list) == 0) {
      new_df$segments[i] <- ""
      next
    }
    
    # Get segments for this assembly
    all_segments <- tryCatch({
      get_segments(row$assembly)
    }, error = function(e) {
      cat("[LEGACY] warning: could not get segments for assembly", row$assembly, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(all_segments) || nrow(all_segments) == 0) {
      new_df$segments[i] <- ""
      next
    }
    
    # Filter segments by contigs
    region_segments <- all_segments[all_segments$contig %in% contigs_list, ]
    
    if (nrow(region_segments) == 0) {
      new_df$segments[i] <- ""
      next
    }
    
    # If zoom is specified, try to narrow down to segments within zoom
    if (!is.na(row$zoom_start) && !is.na(row$zoom_end)) {
      # Save context
      saved_assembly <- cxt_get_assembly()
      saved_segments <- tryCatch(get_state_segments(), error = function(e) NULL)
      saved_zoom <- tryCatch(cxt_get_xlim(), error = function(e) NULL)
      
      tryCatch({
        # Set context for this region
        cxt_set_assembly(row$assembly)
        cxt_set_view(region_segments)
        cxt_set_zoom(c(row$zoom_start, row$zoom_end))
        
        # Get view to see which segments are in zoom
        cdf <- cxt_get_entire_view()
        zoom_segments <- cdf[cdf$vstart < row$zoom_end & cdf$vend > row$zoom_start, ]
        
        if (nrow(zoom_segments) > 0 && "segment" %in% colnames(zoom_segments)) {
          # Use segments that overlap with zoom
          region_segments <- region_segments[region_segments$segment %in% zoom_segments$segment, ]
        } else if (nrow(zoom_segments) > 0 && "segment_ids" %in% colnames(zoom_segments)) {
          # Handle case where segment_ids is comma-separated
          zoom_segment_ids <- unique(unlist(strsplit(zoom_segments$segment_ids, ",")))
          region_segments <- region_segments[region_segments$segment %in% zoom_segment_ids, ]
        }
        
        # Restore context
        if (!is.null(saved_assembly) && !is.null(saved_segments)) {
          cxt_set_assembly(saved_assembly)
          cxt_set_view(saved_segments)
          if (!is.null(saved_zoom)) {
            cxt_set_zoom(saved_zoom)
          }
        }
      }, error = function(e) {
        cat("[LEGACY] warning: error computing zoom segments:", e$message, "\n")
        # Restore context on error
        if (!is.null(saved_assembly) && !is.null(saved_segments)) {
          cxt_set_assembly(saved_assembly)
          cxt_set_view(saved_segments)
          if (!is.null(saved_zoom)) {
            cxt_set_zoom(saved_zoom)
          }
        }
      })
    }
    
    # Convert segment IDs to comma-separated string
    if (nrow(region_segments) > 0) {
      new_df$segments[i] <- paste(region_segments$segment, collapse = ",")
    } else {
      new_df$segments[i] <- ""
    }
  }
  
  # Save converted file
  base_name <- sub("\\.txt$", "", basename(file_path))
  new_file_path <- file.path(regions_dir, paste0(base_name, ".txt"))
  
  # Backup old file
  old_backup_path <- file.path(regions_dir, paste0(base_name, "_old.txt"))
  if (file.exists(new_file_path)) {
    # If new format file already exists, use different backup name
    old_backup_path <- file.path(regions_dir, paste0(base_name, "_old_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  }
  
  file.copy(file_path, old_backup_path)
  cat("[LEGACY] backed up old file to:", basename(old_backup_path), "\n")
  
  # Write new format
  write.table(new_df, new_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[LEGACY] saved converted file to:", basename(new_file_path), "\n")
  
  return(new_df)
}

