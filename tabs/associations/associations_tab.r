# associations tab implementation

# source graph visualization code
source("tabs/associations/association_graph.r")

# validate and extract tab parameters
tab <- get_tab_by_id("associations")
if (is.null(tab)) {
  stop("associations tab not found during loading")
}

# extract required parameters
get_abundance_f <- tab$get_abundance_f
get_abundance_summary_f <- tab$get_abundance_summary_f
get_coverage_f <- tab$get_coverage_f
get_library_ids_f <- tab$get_library_ids_f
get_bin_adj_total_f <- tab$get_bin_adj_total_f
get_bin_adj_associated_f <- tab$get_bin_adj_associated_f
get_bin_adj_support_f <- tab$get_bin_adj_support_f
get_seg_adj_total_f <- tab$get_seg_adj_total_f
get_seg_adj_associated_f <- tab$get_seg_adj_associated_f
get_seg_adj_count_f <- tab$get_seg_adj_count_f
get_bin_segment_table_f <- tab$get_bin_segment_table_f
get_host_table_f <- tab$get_host_table_f

# validate functions
required_params <- c("get_abundance_f", "get_abundance_summary_f", "get_coverage_f", "get_library_ids_f", "get_bin_adj_total_f", 
                    "get_bin_adj_associated_f", "get_bin_adj_support_f",
                    "get_seg_adj_total_f", "get_seg_adj_associated_f", "get_seg_adj_count_f",
                    "get_bin_segment_table_f", "get_host_table_f")
missing_params <- required_params[!sapply(required_params, function(p) {
  param_name <- paste0("get_", gsub("_f$", "", p), "_f")
  if (p == "get_library_ids_f") param_name <- "get_library_ids_f"
  p %in% names(tab) || exists(param_name, envir = environment())
})]
if (length(missing_params) > 0) {
  stop(sprintf("associations tab missing required parameters: %s", paste(missing_params, collapse = ", ")))
}

# reactive value for selected bin
selected_bin <- reactiveVal(cache_get_if_exists("associations.selected_bin", ""))

# reactive value for show bins vs segments
show_bins <- reactiveVal(cache_get_if_exists("associations.show_bins", TRUE))

# reactive value for coverage field
coverage_field <- reactiveVal(cache_get_if_exists("associations.coverage_field", "total"))

# reactive values for graph color schemes
host_color_scheme <- reactiveVal(cache_get_if_exists("associations.host_color_scheme", "none"))
taxa_level <- reactiveVal(cache_get_if_exists("associations.taxa_level", "species"))
element_color_scheme <- reactiveVal(cache_get_if_exists("associations.element_color_scheme", "none"))
edge_color_scheme <- reactiveVal(cache_get_if_exists("associations.edge_color_scheme", "none"))
lib1 <- reactiveVal(cache_get_if_exists("associations.lib1", NULL))
lib2 <- reactiveVal(cache_get_if_exists("associations.lib2", NULL))
show_host_labels <- reactiveVal(cache_get_if_exists("associations.show_host_labels", FALSE))

# set the tab panel UI
set_tab_panel_f(function() {
  tabPanel(
    "Associations",
    fluidRow(
      column(12,
        wellPanel(
          h5("Controls"),
          fluidRow(
            column(3,
              textInput("associationsBinInput", "Bin:", 
                       value = cache_get_if_exists("associations.selected_bin", ""),
                       placeholder = "Enter bin ID (e.g., b1)")
            ),
            column(3,
              radioButtons("associationsShowBinsRadio", "Show:",
                         choices = list("Bins" = TRUE, "Segments" = FALSE),
                         selected = if (cache_get_if_exists("associations.show_bins", TRUE)) TRUE else FALSE,
                         inline = TRUE)
            ),
            column(3,
              selectInput("associationsCoverageField", "Field:",
                         choices = list("Total" = "total", "Associated" = "associated", "Support" = "support"),
                         selected = cache_get_if_exists("associations.coverage_field", "total"),
                         width = "100%")
            ),
            column(3,
              actionButton("associationsClearSelectionBtn", "Clear Selection", 
                         class = "btn btn-default",
                         style = "margin-top: 25px;")
            )
          )
        )
      )
    ),
    fluidRow(
      column(4,
        h4("Abundance Plot"),
        plotly::plotlyOutput("associationsAbundancePlot", height = "300px")
      ),
      column(4,
        h4("Coverage Plot"),
        plotly::plotlyOutput("associationsCoveragePlot", height = "300px")
      ),
      column(4,
        h4("Percents Plot"),
        plotly::plotlyOutput("associationsPercentsPlot", height = "300px")
      )
    ),
    fluidRow(
      column(12,
        tabsetPanel(
          id = "associationsTablesTabset",
          tabPanel("Host Table", DTOutput("associationsHostTable")),
          tabPanel("Segments", DTOutput("associationsSegmentsTable")),
          tabPanel("Graph", 
            fluidRow(
              column(3,
                wellPanel(
                  h5("Graph Controls"),
                  selectInput("associationsHostColorScheme", "Host Color:",
                    choices = list("None" = "none", "Taxa" = "taxa", "Abundance" = "abundance"),
                    selected = cache_get_if_exists("associations.host_color_scheme", "none")),
                  conditionalPanel(
                    condition = "input.associationsHostColorScheme == 'taxa'",
                    selectInput("associationsTaxaLevel", "Taxa Level:",
                      choices = list("Species" = "species", "Genus" = "genus", "Family" = "family", 
                                   "Order" = "order", "Class" = "class", "Phylum" = "phylum", "Domain" = "domain"),
                      selected = cache_get_if_exists("associations.taxa_level", "species"))
                  ),
                  selectInput("associationsElementColorScheme", "Element Color:",
                    choices = list("None" = "none", "Abundance" = "abundance"),
                    selected = cache_get_if_exists("associations.element_color_scheme", "none")),
                  selectInput("associationsEdgeColorScheme", "Edge Color:",
                    choices = list("None" = "none", "Percent" = "percent", "Percent Delta" = "percent_delta"),
                    selected = cache_get_if_exists("associations.edge_color_scheme", "none")),
                  conditionalPanel(
                    condition = "input.associationsEdgeColorScheme == 'percent_delta'",
                    selectInput("associationsLib1", "Library 1:",
                      choices = list(),
                      selected = NULL),
                    selectInput("associationsLib2", "Library 2:",
                      choices = list(),
                      selected = NULL)
                  ),
                  checkboxInput("associationsShowHostLabels", "Show Host Labels",
                    value = cache_get_if_exists("associations.show_host_labels", FALSE))
                )
              ),
              column(7,
                plotly::plotlyOutput("associationsGraphPlot", height = "1000px")
              ),
              column(2,
                plotly::plotlyOutput("associationsGraphLegend", height = "1000px")
              )
            )
          )
        )
      )
    )
  )
})

# ---- Data Loading Functions ----

# get library IDs for current assembly
get_library_ids <- function() {
  if (is.null(state$assembly)) {
    return(c("early", "pre", "post", "late"))  # fallback
  }
  lib_ids <- get_library_ids_f(state$assembly)
  if (is.null(lib_ids) || length(lib_ids) == 0) {
    return(c("early", "pre", "post", "late"))  # fallback
  }
  return(lib_ids)
}

# load abundance data for selected bin
load_abundance_data <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  abundance_mat <- get_abundance_f(state$assembly)
  if (is.null(abundance_mat)) {
    return(NULL)
  }
  
  # filter by bin
  if (!"bin" %in% names(abundance_mat)) {
    return(NULL)
  }
  
  bin_row <- abundance_mat[abundance_mat$bin == bin_id, ]
  if (nrow(bin_row) == 0) {
    return(NULL)
  }
  
  # get library columns (all columns except bin)
  lib_ids <- get_library_ids()
  lib_cols <- names(bin_row)[names(bin_row) != "bin"]
  
  # match library columns to library IDs
  values <- numeric(length(lib_ids))
  for (i in seq_along(lib_ids)) {
    lib_id <- lib_ids[i]
    if (lib_id %in% lib_cols) {
      values[i] <- as.numeric(bin_row[[lib_id]])
    } else {
      values[i] <- 0
    }
  }
  
  return(list(libraries = lib_ids, values = values))
}

# load segments for selected bin
load_bin_segments <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  seg_table <- get_bin_segment_table_f(state$assembly)
  if (is.null(seg_table)) {
    return(NULL)
  }
  
  # filter by bin
  if (!"bin" %in% names(seg_table)) {
    return(NULL)
  }
  
  bin_segments <- seg_table[seg_table$bin == bin_id, ]
  if (nrow(bin_segments) == 0) {
    return(NULL)
  }
  
  # exclude bin_index column if present
  if ("bin_index" %in% names(bin_segments)) {
    bin_segments <- bin_segments[, names(bin_segments) != "bin_index", drop = FALSE]
  }
  
  return(bin_segments)
}

# load associated segments (seg_tgt) with their bin information
load_associated_segments <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  # get segments for this bin
  bin_segments <- load_bin_segments(bin_id)
  if (is.null(bin_segments) || nrow(bin_segments) == 0) {
    return(NULL)
  }
  
  seg_ids <- bin_segments$segment
  
  # get any association matrix to find associated segments
  adj_mat <- get_seg_adj_total_f(state$assembly)
  if (is.null(adj_mat)) {
    return(NULL)
  }
  
  # filter by seg_src in bin segments to get associated segments (seg_tgt)
  if (!"seg_src" %in% names(adj_mat) || !"seg_tgt" %in% names(adj_mat)) {
    return(NULL)
  }
  
  seg_assoc <- adj_mat[adj_mat$seg_src %in% seg_ids, ]
  if (nrow(seg_assoc) == 0) {
    return(NULL)
  }
  
  # get unique associated segments (seg_tgt)
  associated_seg_ids <- unique(seg_assoc$seg_tgt)
  
  # get full segment table to look up bin for each associated segment
  seg_table <- get_bin_segment_table_f(state$assembly)
  if (is.null(seg_table)) {
    return(NULL)
  }
  
  # filter to associated segments and get their bin information
  associated_segments <- seg_table[seg_table$segment %in% associated_seg_ids, ]
  if (nrow(associated_segments) == 0) {
    return(NULL)
  }
  
  # filter out segments where bin == selected_bin (only show segments from other bins)
  if ("bin" %in% names(associated_segments)) {
    associated_segments <- associated_segments[associated_segments$bin != bin_id, ]
    if (nrow(associated_segments) == 0) {
      return(NULL)
    }
  }
  
  # exclude bin_index column if present
  if ("bin_index" %in% names(associated_segments)) {
    associated_segments <- associated_segments[, names(associated_segments) != "bin_index", drop = FALSE]
  }
  
  # reorder columns to put bin first, then segment
  if ("bin" %in% names(associated_segments) && "segment" %in% names(associated_segments)) {
    other_cols <- names(associated_segments)[!names(associated_segments) %in% c("bin", "segment")]
    associated_segments <- associated_segments[, c("bin", "segment", other_cols), drop = FALSE]
  }
  
  # sort by bin, then segment
  if ("bin" %in% names(associated_segments) && "segment" %in% names(associated_segments)) {
    associated_segments <- associated_segments[order(associated_segments$bin, associated_segments$segment), ]
  }
  
  return(associated_segments)
}

# load association data for bins
load_bin_associations <- function(bin_id, field) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  # select appropriate getter function based on field
  getter_f <- switch(field,
    "total" = get_bin_adj_total_f,
    "associated" = get_bin_adj_associated_f,
    "support" = get_bin_adj_support_f,
    NULL
  )
  
  if (is.null(getter_f)) {
    return(NULL)
  }
  
  adj_mat <- getter_f(state$assembly)
  if (is.null(adj_mat)) {
    return(NULL)
  }
  
  # filter by bin_src
  if (!"bin_src" %in% names(adj_mat)) {
    return(NULL)
  }
  
  bin_assoc <- adj_mat[adj_mat$bin_src == bin_id, ]
  if (nrow(bin_assoc) == 0) {
    return(NULL)
  }
  
  # get library columns
  lib_ids <- get_library_ids()
  lib_cols <- names(bin_assoc)[!names(bin_assoc) %in% c("bin_src", "bin_tgt")]
  
  # build result matrix: rows = associated bins, cols = libraries
  result <- list()
  for (i in seq_len(nrow(bin_assoc))) {
    bin_tgt <- bin_assoc$bin_tgt[i]
    values <- numeric(length(lib_ids))
    for (j in seq_along(lib_ids)) {
      lib_id <- lib_ids[j]
      if (lib_id %in% lib_cols) {
        values[j] <- as.numeric(bin_assoc[[lib_id]][i])
      } else {
        values[j] <- 0
      }
    }
    result[[bin_tgt]] <- values
  }
  
  return(list(bins = names(result), data = result, libraries = lib_ids))
}

# load association data for segments
load_seg_associations <- function(bin_id, field) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  # get segments for this bin
  bin_segments <- load_bin_segments(bin_id)
  if (is.null(bin_segments) || nrow(bin_segments) == 0) {
    return(NULL)
  }
  
  seg_ids <- bin_segments$segment
  
  # select appropriate getter function based on field
  getter_f <- switch(field,
    "total" = get_seg_adj_total_f,
    "associated" = get_seg_adj_associated_f,
    "support" = get_seg_adj_count_f,
    NULL
  )
  
  if (is.null(getter_f)) {
    return(NULL)
  }
  
  adj_mat <- getter_f(state$assembly)
  if (is.null(adj_mat)) {
    return(NULL)
  }
  
  # filter by seg_src in bin segments
  if (!"seg_src" %in% names(adj_mat)) {
    return(NULL)
  }
  
  seg_assoc <- adj_mat[adj_mat$seg_src %in% seg_ids, ]
  if (nrow(seg_assoc) == 0) {
    return(NULL)
  }
  
  # get segment table to filter out segments where bin == selected_bin
  seg_table <- get_bin_segment_table_f(state$assembly)
  if (!is.null(seg_table) && "bin" %in% names(seg_table) && "segment" %in% names(seg_table)) {
    # create a lookup: segment -> bin
    seg_to_bin <- setNames(seg_table$bin, seg_table$segment)
    
    # filter out associations where seg_tgt's bin == selected_bin
    seg_tgt_bins <- seg_to_bin[seg_assoc$seg_tgt]
    seg_assoc <- seg_assoc[!is.na(seg_tgt_bins) & seg_tgt_bins != bin_id, ]
    
    if (nrow(seg_assoc) == 0) {
      return(NULL)
    }
  }
  
  # get library columns (exclude key fields and metadata columns)
  lib_ids <- get_library_ids()
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  lib_cols <- names(seg_assoc)[!names(seg_assoc) %in% metadata_cols]
  
  # build result matrix: rows = associated segments, cols = libraries
  result <- list()
  for (i in seq_len(nrow(seg_assoc))) {
    seg_tgt <- seg_assoc$seg_tgt[i]
    seg_src <- seg_assoc$seg_src[i]
    key <- paste0(seg_src, "->", seg_tgt)
    values <- numeric(length(lib_ids))
    for (j in seq_along(lib_ids)) {
      lib_id <- lib_ids[j]
      if (lib_id %in% lib_cols) {
        values[j] <- as.numeric(seg_assoc[[lib_id]][i])
      } else {
        values[j] <- 0
      }
    }
    result[[key]] <- values
  }
  
  return(list(segments = names(result), data = result, libraries = lib_ids))
}

# load bin coverage data (for the bin itself in coverage plot)
load_bin_coverage <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  # use actual bin coverage from BINNING_COV_LR_MAT
  coverage_mat <- get_coverage_f(state$assembly)
  if (is.null(coverage_mat)) {
    return(NULL)
  }
  
  # filter by bin
  if (!"bin" %in% names(coverage_mat)) {
    return(NULL)
  }
  
  bin_row <- coverage_mat[coverage_mat$bin == bin_id, ]
  if (nrow(bin_row) == 0) {
    return(NULL)
  }
  
  # get library columns (all columns except bin)
  lib_ids <- get_library_ids()
  lib_cols <- names(bin_row)[names(bin_row) != "bin"]
  
  # match library columns to library IDs
  values <- numeric(length(lib_ids))
  for (i in seq_along(lib_ids)) {
    lib_id <- lib_ids[i]
    if (lib_id %in% lib_cols) {
      values[i] <- as.numeric(bin_row[[lib_id]])
    } else {
      values[i] <- 0
    }
  }
  
  return(list(libraries = lib_ids, values = values))
}

# load percent data for bins
load_bin_percents <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  # get count and associated_read_count
  count_mat <- get_bin_adj_support_f(state$assembly)
  associated_mat <- get_bin_adj_associated_f(state$assembly)
  
  if (is.null(count_mat) || is.null(associated_mat)) {
    return(NULL)
  }
  
  # filter by bin_src
  if (!"bin_src" %in% names(count_mat) || !"bin_src" %in% names(associated_mat)) {
    return(NULL)
  }
  
  count_rows <- count_mat[count_mat$bin_src == bin_id, ]
  associated_rows <- associated_mat[associated_mat$bin_src == bin_id, ]
  
  if (nrow(count_rows) == 0 || nrow(associated_rows) == 0) {
    return(NULL)
  }
  
  lib_ids <- get_library_ids()
  lib_cols <- names(count_rows)[!names(count_rows) %in% c("bin_src", "bin_tgt")]
  
  # build result: percent = support_read_count / associated_read_count
  result <- list()
  for (i in seq_len(nrow(count_rows))) {
    bin_tgt <- count_rows$bin_tgt[i]
    
    # find matching row in associated_mat
    assoc_row <- associated_rows[associated_rows$bin_tgt == bin_tgt, ]
    if (nrow(assoc_row) == 0) {
      next
    }
    
    values <- numeric(length(lib_ids))
    for (j in seq_along(lib_ids)) {
      lib_id <- lib_ids[j]
      if (lib_id %in% lib_cols) {
        count_val <- as.numeric(count_rows[[lib_id]][i])
        assoc_val <- as.numeric(assoc_row[[lib_id]][1])
        if (assoc_val > 0) {
          values[j] <- (count_val / assoc_val) * 100
        } else {
          values[j] <- 0
        }
      } else {
        values[j] <- 0
      }
    }
    result[[bin_tgt]] <- values
  }
  
  return(list(bins = names(result), data = result, libraries = lib_ids))
}

# load percent data for segments
load_seg_percents <- function(bin_id) {
  if (is.null(bin_id) || bin_id == "" || is.null(state$assembly)) {
    return(NULL)
  }
  
  # get segments for this bin
  bin_segments <- load_bin_segments(bin_id)
  if (is.null(bin_segments) || nrow(bin_segments) == 0) {
    return(NULL)
  }
  
  seg_ids <- bin_segments$segment
  
  # get count and associated_read_count
  count_mat <- get_seg_adj_count_f(state$assembly)
  associated_mat <- get_seg_adj_associated_f(state$assembly)
  
  if (is.null(count_mat) || is.null(associated_mat)) {
    return(NULL)
  }
  
  # filter by seg_src in bin segments
  if (!"seg_src" %in% names(count_mat) || !"seg_src" %in% names(associated_mat)) {
    return(NULL)
  }
  
  count_rows <- count_mat[count_mat$seg_src %in% seg_ids, ]
  associated_rows <- associated_mat[associated_mat$seg_src %in% seg_ids, ]
  
  if (nrow(count_rows) == 0 || nrow(associated_rows) == 0) {
    return(NULL)
  }
  
  # get segment table to filter out segments where bin == selected_bin
  seg_table <- get_bin_segment_table_f(state$assembly)
  if (!is.null(seg_table) && "bin" %in% names(seg_table) && "segment" %in% names(seg_table)) {
    # create a lookup: segment -> bin
    seg_to_bin <- setNames(seg_table$bin, seg_table$segment)
    
    # filter out associations where seg_tgt's bin == selected_bin
    seg_tgt_bins <- seg_to_bin[count_rows$seg_tgt]
    valid_rows <- !is.na(seg_tgt_bins) & seg_tgt_bins != bin_id
    count_rows <- count_rows[valid_rows, ]
    associated_rows <- associated_rows[associated_rows$seg_tgt %in% count_rows$seg_tgt, ]
    
    if (nrow(count_rows) == 0 || nrow(associated_rows) == 0) {
      return(NULL)
    }
  }
  
  lib_ids <- get_library_ids()
  metadata_cols <- c("seg_src", "seg_tgt", "side_src", "side_tgt", 
                     "contig_src", "start_src", "end_src", 
                     "contig_tgt", "start_tgt", "end_tgt")
  lib_cols <- names(count_rows)[!names(count_rows) %in% metadata_cols]
  
  # build result: percent = count / associated_read_count
  result <- list()
  for (i in seq_len(nrow(count_rows))) {
    seg_tgt <- count_rows$seg_tgt[i]
    seg_src <- count_rows$seg_src[i]
    key <- paste0(seg_src, "->", seg_tgt)
    
    # find matching row in associated_mat (match on seg_src, seg_tgt, and sides if present)
    if ("side_src" %in% names(count_rows) && "side_tgt" %in% names(count_rows) &&
        "side_src" %in% names(associated_rows) && "side_tgt" %in% names(associated_rows)) {
      side_src <- count_rows$side_src[i]
      side_tgt <- count_rows$side_tgt[i]
      assoc_row <- associated_rows[associated_rows$seg_src == seg_src & 
                                   associated_rows$seg_tgt == seg_tgt &
                                   associated_rows$side_src == side_src &
                                   associated_rows$side_tgt == side_tgt, ]
    } else {
      assoc_row <- associated_rows[associated_rows$seg_src == seg_src & associated_rows$seg_tgt == seg_tgt, ]
    }
    if (nrow(assoc_row) == 0) {
      next
    }
    
    values <- numeric(length(lib_ids))
    for (j in seq_along(lib_ids)) {
      lib_id <- lib_ids[j]
      if (lib_id %in% lib_cols) {
        count_val <- as.numeric(count_rows[[lib_id]][i])
        assoc_val <- as.numeric(assoc_row[[lib_id]][1])
        if (assoc_val > 0) {
          values[j] <- (count_val / assoc_val) * 100
        } else {
          values[j] <- 0
        }
      } else {
        values[j] <- 0
      }
    }
    result[[key]] <- values
  }
  
  return(list(segments = names(result), data = result, libraries = lib_ids))
}

# ---- Event Handlers ----

# observer for bin input change
observeEvent(input$associationsBinInput, {
  bin_id <- trimws(input$associationsBinInput)
  selected_bin(bin_id)
  cache_set("associations.selected_bin", bin_id)
})

# observer for show bins radio buttons
observeEvent(input$associationsShowBinsRadio, {
  show_bins(input$associationsShowBinsRadio)
  cache_set("associations.show_bins", input$associationsShowBinsRadio)
})

# observer for coverage field selector
observeEvent(input$associationsCoverageField, {
  coverage_field(input$associationsCoverageField)
  cache_set("associations.coverage_field", input$associationsCoverageField)
})

# observer for clear selection button
observeEvent(input$associationsClearSelectionBtn, {
  # clear selection in segments table
  proxy <- DT::dataTableProxy("associationsSegmentsTable")
  DT::selectRows(proxy, NULL)
})

# observers for graph color scheme inputs
observeEvent(input$associationsHostColorScheme, {
  host_color_scheme(input$associationsHostColorScheme)
  cache_set("associations.host_color_scheme", input$associationsHostColorScheme)
})

observeEvent(input$associationsTaxaLevel, {
  taxa_level(input$associationsTaxaLevel)
  cache_set("associations.taxa_level", input$associationsTaxaLevel)
})

observeEvent(input$associationsElementColorScheme, {
  element_color_scheme(input$associationsElementColorScheme)
  cache_set("associations.element_color_scheme", input$associationsElementColorScheme)
})

observeEvent(input$associationsEdgeColorScheme, {
  edge_color_scheme(input$associationsEdgeColorScheme)
  cache_set("associations.edge_color_scheme", input$associationsEdgeColorScheme)
})

observeEvent(input$associationsLib1, {
  lib1(input$associationsLib1)
  cache_set("associations.lib1", input$associationsLib1)
})

observeEvent(input$associationsLib2, {
  lib2(input$associationsLib2)
  cache_set("associations.lib2", input$associationsLib2)
})

observeEvent(input$associationsShowHostLabels, {
  show_host_labels(input$associationsShowHostLabels)
  cache_set("associations.show_host_labels", input$associationsShowHostLabels)
})

# observer to populate library dropdowns when assembly changes
observe({
  if (is.null(state$assembly)) {
    return()
  }
  
  lib_ids <- get_library_ids()
  if (is.null(lib_ids) || length(lib_ids) == 0) {
    return()
  }
  
  # create choices with assembly prefix (e.g., EBC_pre)
  lib_choices <- setNames(paste0(state$assembly, "_", lib_ids), lib_ids)
  
  # set default to pre and post if available
  default_lib1 <- NULL
  default_lib2 <- NULL
  if ("pre" %in% lib_ids) {
    default_lib1 <- paste0(state$assembly, "_pre")
  } else if (length(lib_ids) > 0) {
    default_lib1 <- lib_choices[1]
  }
  if ("post" %in% lib_ids) {
    default_lib2 <- paste0(state$assembly, "_post")
  } else if (length(lib_ids) > 1) {
    default_lib2 <- lib_choices[2]
  }
  
  # update dropdowns if session is available
  if (exists("session", envir = parent.frame())) {
    session <- get("session", envir = parent.frame())
    updateSelectInput(session, "associationsLib1", choices = lib_choices, selected = default_lib1)
    updateSelectInput(session, "associationsLib2", choices = lib_choices, selected = default_lib2)
  }
  
  # update reactive values if not already set
  if (is.null(lib1())) {
    lib1(default_lib1)
    cache_set("associations.lib1", default_lib1)
  }
  if (is.null(lib2())) {
    lib2(default_lib2)
    cache_set("associations.lib2", default_lib2)
  }
})

# observer for select host button
observeEvent(input$associationsSelectHostTrigger, {
  row_idx <- as.integer(input$associationsSelectHostTrigger)
  host_table <- get_host_table_f(state$assembly)
  
  if (is.null(host_table)) {
    showNotification("Invalid host selection", type = "error")
    return()
  }
  
  # apply same processing as in table renderer: merge abundance, remove domain, sort
  abundance_summary <- get_abundance_summary_f(state$assembly)
  if (!is.null(abundance_summary) && "bin" %in% names(abundance_summary) && "mean_abundance" %in% names(abundance_summary)) {
    host_table <- merge(host_table, abundance_summary[, c("bin", "mean_abundance")], 
                       by = "bin", all.x = TRUE, sort = FALSE)
  }
  
  # remove domain column if present
  display_table <- host_table
  if ("domain" %in% names(display_table)) {
    display_table <- display_table[, names(display_table) != "domain", drop = FALSE]
  }
  
  # sort by mean_abundance descending (high to low), handling NAs
  if ("mean_abundance" %in% names(display_table)) {
    display_table <- display_table[order(-display_table$mean_abundance, na.last = TRUE), ]
  }
  
  if (row_idx <= 0 || row_idx > nrow(display_table)) {
    showNotification("Invalid host selection", type = "error")
    return()
  }
  
  selected_host <- display_table[row_idx, ]
  if (!"bin" %in% names(selected_host)) {
    showNotification("Host table missing bin column", type = "error")
    return()
  }
  
  bin_id <- as.character(selected_host$bin)
  selected_bin(bin_id)
  cache_set("associations.selected_bin", bin_id)
  
  showNotification(sprintf("Selected bin: %s", bin_id), type = "message")
})

# observer for goto segment button
observeEvent(input$associationsGotoSegmentTrigger, {
  row_idx <- as.integer(input$associationsGotoSegmentTrigger)
  bin_id <- selected_bin()
  
  if (is.null(bin_id) || bin_id == "") {
    showNotification("No bin selected", type = "warning")
    return()
  }
  
  segments_table <- load_associated_segments(bin_id)
  if (is.null(segments_table) || row_idx <= 0 || row_idx > nrow(segments_table)) {
    showNotification("Invalid segment selection", type = "error")
    return()
  }
  
  selected_segment <- segments_table[row_idx, ]
  
  # check required columns
  if (!all(c("contig", "start", "end") %in% names(selected_segment))) {
    showNotification("Segment missing required columns (contig, start, end)", type = "error")
    return()
  }
  
  contig_name <- as.character(selected_segment$contig)
  seg_start <- as.numeric(selected_segment$start)
  seg_end <- as.numeric(selected_segment$end)
  
  # build context to convert to global coordinates
  contigs_table <- get_contigs(state$assembly)
  if (is.null(contigs_table)) {
    showNotification("Cannot get contig information", type = "error")
    return()
  }
  
  cxt <- build_context(c(contig_name), contigs_table, NULL, state$assembly)
  if (is.null(cxt)) {
    showNotification("Cannot build context for navigation", type = "error")
    return()
  }
  
  # convert to global coordinates
  gstart <- cxt$mapper$l2g(contig_name, seg_start)
  gend <- cxt$mapper$l2g(contig_name, seg_end)
  
  # calculate zoom with margin
  span <- gend - gstart
  margin <- span * 0.1
  zoom_start <- gstart - margin
  zoom_end <- gend + margin
  
  # push current region to undo before changing
  regions_module_output$push_undo_state()
  
  # set zoom to calculated range
  state$contigs <- c(contig_name)
  state$zoom <- c(zoom_start, zoom_end)
  
  showNotification(sprintf("Navigated to segment %s", selected_segment$segment[1]), type = "message")
})

# ---- Plot Renderers ----

# abundance plot
output$associationsAbundancePlot <- plotly::renderPlotly({
  bin_id <- selected_bin()
  abundance_data <- load_abundance_data(bin_id)
  
  if (is.null(abundance_data)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No abundance data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  plot_data <- data.frame(
    library = factor(abundance_data$libraries, levels = abundance_data$libraries),
    abundance = abundance_data$values
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = abundance, group = 1)) +
    ggplot2::geom_line(color = "#2c7bb6", size = 1.5) +
    ggplot2::geom_point(color = "#2c7bb6", size = 3) +
    ggplot2::labs(x = "Library", y = "Abundance (%)", title = sprintf("Bin %s Abundance", bin_id)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  plotly::ggplotly(p)
})

# coverage plot
output$associationsCoveragePlot <- plotly::renderPlotly({
  bin_id <- selected_bin()
  field <- coverage_field()
  show_bins_flag <- as.logical(show_bins())
  if (is.na(show_bins_flag)) show_bins_flag <- TRUE
  
  if (is.null(bin_id) || bin_id == "") {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Please select a bin", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # get bin coverage (thick dark grey line)
  bin_coverage <- load_bin_coverage(bin_id)
  
  # get associated data
  if (show_bins_flag) {
    assoc_data <- load_bin_associations(bin_id, field)
  } else {
    assoc_data <- load_seg_associations(bin_id, field)
  }
  
  # build plot data
  plot_data <- data.frame()
  
  # add bin line
  if (!is.null(bin_coverage)) {
    bin_df <- data.frame(
      library = factor(bin_coverage$libraries, levels = bin_coverage$libraries),
      value = bin_coverage$values,
      item = bin_id,
      item_type = "bin",
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, bin_df)
  }
  
  # add associated items
  if (!is.null(assoc_data)) {
    lib_ids <- assoc_data$libraries
    for (item_name in assoc_data$bins %||% assoc_data$segments) {
      values <- assoc_data$data[[item_name]]
      item_df <- data.frame(
        library = factor(lib_ids, levels = lib_ids),
        value = values,
        item = item_name,
        item_type = if (show_bins_flag) "associated_bin" else "associated_segment",
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, item_df)
    }
  }
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No coverage data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # get selected segment from segments table (if showing segments)
  selected_segment_id <- NULL
  if (!show_bins_flag) {
    selected_rows <- input$associationsSegmentsTable_rows_selected
    if (!is.null(selected_rows) && length(selected_rows) > 0) {
      segments_table <- load_associated_segments(bin_id)
      if (!is.null(segments_table) && selected_rows[1] <= nrow(segments_table)) {
        selected_segment_id <- as.character(segments_table$segment[selected_rows[1]])
      }
    }
  }
  
  # create colors: dark grey for bin, colors for associated items
  unique_items <- unique(plot_data$item)
  colors <- rainbow(length(unique_items))
  names(colors) <- unique_items
  colors[bin_id] <- "#666666"  # dark grey for bin
  
  # highlight selected segment if present (keep original color, just make it thicker)
  selected_item_key <- NULL
  if (!is.null(selected_segment_id) && !show_bins_flag) {
    # find the item key that matches the selected segment (check seg_tgt part)
    for (item_key in unique_items) {
      if (item_key != bin_id && grepl("->", item_key)) {
        seg_tgt <- strsplit(item_key, "->")[[1]][2]
        if (seg_tgt == selected_segment_id) {
          selected_item_key <- item_key
          # keep original color, don't change it
          break
        }
      }
    }
  }
  
  # add key field for plotly click events
  plot_data$key <- plot_data$item
  
  # separate selected segment from other items for highlighting
  if (!is.null(selected_item_key)) {
    # create filtered data frames
    bin_data <- plot_data[plot_data$item == bin_id, ]
    selected_data <- plot_data[plot_data$item == selected_item_key, ]
    other_data <- plot_data[plot_data$item != bin_id & plot_data$item != selected_item_key, ]
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = item, group = item, key = key)) +
      ggplot2::scale_color_manual(values = colors)
    
    # add bin line if present
    if (nrow(bin_data) > 0) {
      p <- p + ggplot2::geom_line(data = bin_data, size = 3, alpha = 0.8)
    }
    
    # add selected segment line if present
    if (nrow(selected_data) > 0) {
      p <- p + ggplot2::geom_line(data = selected_data, size = 3, alpha = 0.9, linetype = "solid") +
        ggplot2::geom_point(data = selected_data, size = 4)
    }
    
    # add other items lines if present
    if (nrow(other_data) > 0) {
      p <- p + ggplot2::geom_line(data = other_data, size = 1, alpha = 0.6) +
        ggplot2::geom_point(data = other_data, size = 2)
    }
    
    p <- p +
      ggplot2::labs(x = "Library", y = sprintf("Coverage (%s)", field), 
                  title = sprintf("Bin %s Coverage", bin_id), color = "Item") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "right"
      )
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = item, group = item, key = key)) +
      ggplot2::geom_line(data = plot_data[plot_data$item == bin_id, ], size = 3, alpha = 0.8) +
      ggplot2::geom_line(data = plot_data[plot_data$item != bin_id, ], size = 1, alpha = 0.6) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(x = "Library", y = sprintf("Coverage (%s)", field), 
                  title = sprintf("Bin %s Coverage", bin_id), color = "Item") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "right"
      )
  }
  
  plotly::ggplotly(p, source = "associations_coverage_plot")
})

# percents plot
output$associationsPercentsPlot <- plotly::renderPlotly({
  bin_id <- selected_bin()
  show_bins_flag <- as.logical(show_bins())
  if (is.na(show_bins_flag)) show_bins_flag <- TRUE
  
  if (is.null(bin_id) || bin_id == "") {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Please select a bin", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # get percent data (only associated items, not bin itself)
  if (show_bins_flag) {
    percent_data <- load_bin_percents(bin_id)
  } else {
    percent_data <- load_seg_percents(bin_id)
  }
  
  if (is.null(percent_data)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No percent data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # build plot data
  plot_data <- data.frame()
  lib_ids <- percent_data$libraries
  
  for (item_name in percent_data$bins %||% percent_data$segments) {
    values <- percent_data$data[[item_name]]
    item_df <- data.frame(
      library = factor(lib_ids, levels = lib_ids),
      percent = values,
      item = item_name,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, item_df)
  }
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No percent data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # get selected segment from segments table (if showing segments)
  selected_segment_id <- NULL
  if (!show_bins_flag) {
    selected_rows <- input$associationsSegmentsTable_rows_selected
    if (!is.null(selected_rows) && length(selected_rows) > 0) {
      segments_table <- load_associated_segments(bin_id)
      if (!is.null(segments_table) && selected_rows[1] <= nrow(segments_table)) {
        selected_segment_id <- as.character(segments_table$segment[selected_rows[1]])
      }
    }
  }
  
  # create colors for items
  unique_items <- unique(plot_data$item)
  colors <- rainbow(length(unique_items))
  names(colors) <- unique_items
  
  # highlight selected segment if present (keep original color, just make it thicker)
  selected_item_key <- NULL
  if (!is.null(selected_segment_id) && !show_bins_flag) {
    # find the item key that matches the selected segment (check seg_tgt part)
    for (item_key in unique_items) {
      if (grepl("->", item_key)) {
        seg_tgt <- strsplit(item_key, "->")[[1]][2]
        if (seg_tgt == selected_segment_id) {
          selected_item_key <- item_key
          # keep original color, don't change it
          break
        }
      }
    }
  }
  
  # add key field for plotly click events
  plot_data$key <- plot_data$item
  
  # separate selected segment from other items for highlighting
  if (!is.null(selected_item_key)) {
    # create filtered data frames
    selected_data <- plot_data[plot_data$item == selected_item_key, ]
    other_data <- plot_data[plot_data$item != selected_item_key, ]
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = percent, color = item, group = item, key = key)) +
      ggplot2::scale_color_manual(values = colors)
    
    # add selected segment line if present
    if (nrow(selected_data) > 0) {
      p <- p + ggplot2::geom_line(data = selected_data, size = 3, alpha = 0.9, linetype = "solid") +
        ggplot2::geom_point(data = selected_data, size = 4)
    }
    
    # add other items lines if present
    if (nrow(other_data) > 0) {
      p <- p + ggplot2::geom_line(data = other_data, size = 1, alpha = 0.7) +
        ggplot2::geom_point(data = other_data, size = 2)
    }
    
    p <- p +
      ggplot2::labs(x = "Library", y = "Percent (%)", 
                  title = sprintf("Bin %s Associated Items Percentages", bin_id), color = "Item") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "right"
      )
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = percent, color = item, group = item, key = key)) +
      ggplot2::geom_line(size = 1, alpha = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(x = "Library", y = "Percent (%)", 
                  title = sprintf("Bin %s Associated Items Percentages", bin_id), color = "Item") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "right"
      )
  }
  
  plotly::ggplotly(p, source = "associations_percents_plot")
})

# ---- Plot Click Handlers ----

# click observer for coverage plot
observeEvent(plotly::event_data("plotly_click", source = "associations_coverage_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "associations_coverage_plot")
  if (is.null(event_data) || is.null(event_data$key)) {
    return()
  }
  
  clicked_item <- as.character(event_data$key)
  bin_id <- selected_bin()
  show_bins_flag <- show_bins()
  
  if (show_bins_flag) {
    # clicked on a bin - select it in host table
    host_table <- get_host_table_f(state$assembly)
    if (is.null(host_table)) {
      return()
    }
    
    # apply same processing as in table renderer
    abundance_summary <- get_abundance_summary_f(state$assembly)
    if (!is.null(abundance_summary) && "bin" %in% names(abundance_summary) && "mean_abundance" %in% names(abundance_summary)) {
      host_table <- merge(host_table, abundance_summary[, c("bin", "mean_abundance")], 
                         by = "bin", all.x = TRUE, sort = FALSE)
    }
    
    display_table <- host_table
    if ("domain" %in% names(display_table)) {
      display_table <- display_table[, names(display_table) != "domain", drop = FALSE]
    }
    
    if ("mean_abundance" %in% names(display_table)) {
      display_table <- display_table[order(-display_table$mean_abundance, na.last = TRUE), ]
    }
    
    # find matching row
    matching_rows <- which(display_table$bin == clicked_item)
    if (length(matching_rows) > 0) {
      # switch to host table tab
      if (exists("session", envir = parent.frame())) {
        session <- get("session", envir = parent.frame())
        updateTabsetPanel(session, "associationsTablesTabset", selected = "Host Table")
      }
      proxy <- DT::dataTableProxy("associationsHostTable")
      DT::selectRows(proxy, matching_rows[1])
    }
  } else {
    # clicked on a segment - select it in segments table
    # key format is "seg_src->seg_tgt", extract seg_tgt
    if (grepl("->", clicked_item)) {
      seg_tgt <- strsplit(clicked_item, "->")[[1]][2]
    } else {
      seg_tgt <- clicked_item
    }
    
    segments_table <- load_associated_segments(bin_id)
    if (is.null(segments_table)) {
      return()
    }
    
    # find matching row by segment
    matching_rows <- which(segments_table$segment == seg_tgt)
    if (length(matching_rows) > 0) {
      # switch to segments table tab
      if (exists("session", envir = parent.frame())) {
        session <- get("session", envir = parent.frame())
        updateTabsetPanel(session, "associationsTablesTabset", selected = "Segments")
      }
      proxy <- DT::dataTableProxy("associationsSegmentsTable")
      DT::selectRows(proxy, matching_rows[1])
    }
  }
})

# click observer for percents plot
observeEvent(plotly::event_data("plotly_click", source = "associations_percents_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "associations_percents_plot")
  if (is.null(event_data) || is.null(event_data$key)) {
    return()
  }
  
  clicked_item <- as.character(event_data$key)
  bin_id <- selected_bin()
  show_bins_flag <- show_bins()
  
  if (show_bins_flag) {
    # clicked on a bin - select it in host table
    host_table <- get_host_table_f(state$assembly)
    if (is.null(host_table)) {
      return()
    }
    
    # apply same processing as in table renderer
    abundance_summary <- get_abundance_summary_f(state$assembly)
    if (!is.null(abundance_summary) && "bin" %in% names(abundance_summary) && "mean_abundance" %in% names(abundance_summary)) {
      host_table <- merge(host_table, abundance_summary[, c("bin", "mean_abundance")], 
                         by = "bin", all.x = TRUE, sort = FALSE)
    }
    
    display_table <- host_table
    if ("domain" %in% names(display_table)) {
      display_table <- display_table[, names(display_table) != "domain", drop = FALSE]
    }
    
    if ("mean_abundance" %in% names(display_table)) {
      display_table <- display_table[order(-display_table$mean_abundance, na.last = TRUE), ]
    }
    
    # find matching row
    matching_rows <- which(display_table$bin == clicked_item)
    if (length(matching_rows) > 0) {
      # switch to host table tab
      if (exists("session", envir = parent.frame())) {
        session <- get("session", envir = parent.frame())
        updateTabsetPanel(session, "associationsTablesTabset", selected = "Host Table")
      }
      proxy <- DT::dataTableProxy("associationsHostTable")
      DT::selectRows(proxy, matching_rows[1])
    }
  } else {
    # clicked on a segment - select it in segments table
    # key format is "seg_src->seg_tgt", extract seg_tgt
    if (grepl("->", clicked_item)) {
      seg_tgt <- strsplit(clicked_item, "->")[[1]][2]
    } else {
      seg_tgt <- clicked_item
    }
    
    segments_table <- load_associated_segments(bin_id)
    if (is.null(segments_table)) {
      return()
    }
    
    # find matching row by segment
    matching_rows <- which(segments_table$segment == seg_tgt)
    if (length(matching_rows) > 0) {
      # switch to segments table tab
      if (exists("session", envir = parent.frame())) {
        session <- get("session", envir = parent.frame())
        updateTabsetPanel(session, "associationsTablesTabset", selected = "Segments")
      }
      proxy <- DT::dataTableProxy("associationsSegmentsTable")
      DT::selectRows(proxy, matching_rows[1])
    }
  }
})

# ---- Graph Renderer ----

# graph plot
output$associationsGraphPlot <- plotly::renderPlotly({
  render_associations_graph(
    state$assembly, 
    selected_bin(), 
    get_host_table_f, 
    get_bin_adj_total_f,
    get_abundance_summary_f,
    get_bin_adj_support_f,
    get_bin_adj_associated_f,
    get_bin_segment_table_f,
    host_color_scheme(),
    taxa_level(),
    element_color_scheme(),
    edge_color_scheme(),
    lib1(),
    lib2(),
    show_host_labels()
  )
})

# graph legend plot (refresh when graph is plotted)
output$associationsGraphLegend <- plotly::renderPlotly({
  # depend on same reactive values as graph to ensure refresh
  selected_bin()  # dependency to refresh when graph updates
  state$assembly
  
  render_associations_legend(
    state$assembly,
    get_host_table_f,
    get_abundance_summary_f,
    get_bin_adj_total_f,
    host_color_scheme(),
    taxa_level(),
    edge_color_scheme()
  )
})

# click observer for graph plot
observeEvent(plotly::event_data("plotly_click", source = "associations_graph_plot"), {
  event_data <- plotly::event_data("plotly_click", source = "associations_graph_plot")
  if (is.null(event_data) || is.null(event_data$key)) {
    return()
  }
  
  clicked_bin <- as.character(event_data$key)
  if (!is.null(clicked_bin) && clicked_bin != "") {
    # update selected bin
    selected_bin(clicked_bin)
    cache_set("associations.selected_bin", clicked_bin)
    
    # update text input if session is available
    if (exists("session", envir = parent.frame())) {
      session <- get("session", envir = parent.frame())
      updateTextInput(session, "associationsBinInput", value = clicked_bin)
    }
  }
})

# ---- Table Renderers ----

# segments table
output$associationsSegmentsTable <- renderDT({
  bin_id <- selected_bin()
  segments_table <- load_associated_segments(bin_id)
  
  if (is.null(segments_table) || nrow(segments_table) == 0) {
    return(datatable(
      data.frame(Message = "No associated segments found. Please select a bin."),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # add goto button column
  display_table <- segments_table
  display_table$Actions <- sapply(seq_len(nrow(display_table)), function(i) {
    as.character(actionButton(
      paste0("goto_seg_btn_", i),
      "Goto",
      class = "btn btn-primary btn-sm",
      onclick = sprintf("Shiny.setInputValue('associationsGotoSegmentTrigger', %d, {priority: 'event'})", i)
    ))
  })
  
  datatable(
    display_table,
    rownames = FALSE,
    escape = FALSE,
    selection = "single",
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip",
      columnDefs = list(
        list(targets = ncol(display_table) - 1, orderable = FALSE, width = "80px")
      )
    ),
    filter = "top"
  )
})

# host table
output$associationsHostTable <- renderDT({
  if (is.null(state$assembly)) {
    return(datatable(
      data.frame(Message = "No assembly selected"),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  host_table <- get_host_table_f(state$assembly)
  
  if (is.null(host_table) || nrow(host_table) == 0) {
    return(datatable(
      data.frame(Message = "No host data available"),
      options = list(dom = "t"),
      rownames = FALSE,
      selection = "none"
    ))
  }
  
  # get abundance summary and merge mean_abundance
  abundance_summary <- get_abundance_summary_f(state$assembly)
  if (!is.null(abundance_summary) && "bin" %in% names(abundance_summary) && "mean_abundance" %in% names(abundance_summary)) {
    # merge mean_abundance by bin
    host_table <- merge(host_table, abundance_summary[, c("bin", "mean_abundance")], 
                       by = "bin", all.x = TRUE, sort = FALSE)
  } else {
    # add empty mean_abundance column if summary not available
    host_table$mean_abundance <- NA
  }
  
  # remove domain column if present
  display_table <- host_table
  if ("domain" %in% names(display_table)) {
    display_table <- display_table[, names(display_table) != "domain", drop = FALSE]
  }
  
  # sort by mean_abundance descending (high to low), handling NAs
  if ("mean_abundance" %in% names(display_table)) {
    display_table <- display_table[order(-display_table$mean_abundance, na.last = TRUE), ]
  }
  
  # add select button column
  display_table$Actions <- sapply(seq_len(nrow(display_table)), function(i) {
    as.character(actionButton(
      paste0("select_host_btn_", i),
      "Select",
      class = "btn btn-secondary btn-sm",
      onclick = sprintf("Shiny.setInputValue('associationsSelectHostTrigger', %d, {priority: 'event'})", i)
    ))
  })
  
  datatable(
    display_table,
    rownames = FALSE,
    escape = FALSE,
    selection = "single",
    options = list(
      pageLength = 15,
      lengthMenu = c(5, 10, 15, 25, 50),
      scrollX = TRUE,
      dom = "lftip",
      columnDefs = list(
        list(targets = ncol(display_table) - 1, orderable = FALSE, width = "80px")
      )
    ),
    filter = "top"
  )
})

