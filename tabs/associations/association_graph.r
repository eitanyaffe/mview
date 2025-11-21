# association graph visualization
# This file contains functions for rendering the host-element graph

# function to load graph data (hosts, elements, edges)
load_graph_data <- function(assembly, get_host_table_f, get_bin_adj_total_f) {
  if (is.null(assembly)) {
    return(NULL)
  }
  
  # get hosts
  host_table <- get_host_table_f(assembly)
  if (is.null(host_table) || nrow(host_table) == 0) {
    return(NULL)
  }
  host_ids <- unique(host_table$bin)
  
  # get adjacency matrix to find elements
  adj_mat <- get_bin_adj_total_f(assembly)
  if (is.null(adj_mat)) {
    return(NULL)
  }
  
  # filter to edges where bin_src is a host
  if (!"bin_src" %in% names(adj_mat) || !"bin_tgt" %in% names(adj_mat)) {
    return(NULL)
  }
  
  host_edges <- adj_mat[adj_mat$bin_src %in% host_ids, ]
  if (nrow(host_edges) == 0) {
    return(NULL)
  }
  
  # get unique element bins (bin_tgt that are not hosts)
  element_ids <- unique(host_edges$bin_tgt)
  element_ids <- element_ids[!element_ids %in% host_ids]
  
  # create nodes data frame
  nodes <- data.frame(
    id = c(host_ids, element_ids),
    type = c(rep("host", length(host_ids)), rep("element", length(element_ids))),
    stringsAsFactors = FALSE
  )
  
  # create edges data frame (only host->element edges)
  edges <- host_edges[host_edges$bin_tgt %in% element_ids, c("bin_src", "bin_tgt")]
  names(edges) <- c("from", "to")
  
  return(list(nodes = nodes, edges = edges, host_table = host_table))
}

# function to compute graph layout using igraph
compute_graph_layout <- function(graph_data) {
  if (is.null(graph_data) || nrow(graph_data$edges) == 0) {
    return(NULL)
  }
  
  # check for igraph package
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for graph visualization. Install with: install.packages('igraph')")
  }
  
  # create igraph object
  g <- igraph::graph_from_data_frame(graph_data$edges, directed = FALSE, vertices = graph_data$nodes)
  
  # compute force-directed layout
  layout <- igraph::layout_with_fr(g)
  
  # add layout coordinates to nodes
  nodes_with_layout <- graph_data$nodes
  nodes_with_layout$x <- layout[, 1]
  nodes_with_layout$y <- layout[, 2]
  
  # add coordinates to edges
  edges_with_coords <- graph_data$edges
  edges_with_coords$x0 <- nodes_with_layout$x[match(edges_with_coords$from, nodes_with_layout$id)]
  edges_with_coords$y0 <- nodes_with_layout$y[match(edges_with_coords$from, nodes_with_layout$id)]
  edges_with_coords$x1 <- nodes_with_layout$x[match(edges_with_coords$to, nodes_with_layout$id)]
  edges_with_coords$y1 <- nodes_with_layout$y[match(edges_with_coords$to, nodes_with_layout$id)]
  
  return(list(nodes = nodes_with_layout, edges = edges_with_coords, host_table = graph_data$host_table))
}

# function to create hover text for nodes
create_node_hover <- function(node_id, node_type, host_table, bin_lengths = NULL) {
  if (node_type == "host" && !is.null(host_table)) {
    host_info <- host_table[host_table$bin == node_id, ]
    if (nrow(host_info) > 0) {
      info_parts <- c(sprintf("Bin: %s", node_id))
      if ("species" %in% names(host_info) && !is.na(host_info$species[1]) && host_info$species[1] != "") {
        info_parts <- c(info_parts, sprintf("Species: %s", host_info$species[1]))
      }
      if ("length" %in% names(host_info)) {
        info_parts <- c(info_parts, sprintf("Length: %s", format(host_info$length[1], big.mark = ",")))
      }
      if ("gc_content" %in% names(host_info)) {
        info_parts <- c(info_parts, sprintf("GC: %.2f%%", host_info$gc_content[1] * 100))
      }
      return(paste(info_parts, collapse = "<br>"))
    }
  } else if (node_type == "element") {
    info_parts <- c(sprintf("Bin: %s", node_id), "Type: Element")
    if (!is.null(bin_lengths) && node_id %in% names(bin_lengths)) {
      info_parts <- c(info_parts, sprintf("Length: %s", format(bin_lengths[node_id], big.mark = ",")))
    }
    return(paste(info_parts, collapse = "<br>"))
  }
  return(sprintf("Bin: %s<br>Type: %s", node_id, node_type))
}

# function to process host_table: sort by taxa, assign colors, add color column
process_host_table_taxa <- function(host_table, taxa_level) {
  if (is.null(host_table) || !taxa_level %in% names(host_table)) {
    return(host_table)
  }
  
  # work on a copy
  result_table <- host_table
  
  # replace empty strings with "unknown"
  result_table[[taxa_level]][result_table[[taxa_level]] == ""] <- "unknown"
  
  # sort by the selected taxa level only
  result_table <- result_table[order(result_table[[taxa_level]]), ]
  
  # get unique taxa values in sorted order
  taxa_values_vec <- as.character(result_table[[taxa_level]])
  unique_taxa <- unique(taxa_values_vec)
  
  # assign colors to unique taxa
  if (length(unique_taxa) > 0) {
    taxa_colors <- rainbow(length(unique_taxa))
    names(taxa_colors) <- unique_taxa
    
    # set "unknown" to white
    if ("unknown" %in% names(taxa_colors)) {
      taxa_colors["unknown"] <- "#FFFFFF"
    }
    
    # add color column to result_table
    result_table$taxa_color <- taxa_colors[taxa_values_vec]
  } else {
    result_table$taxa_color <- "#CCCCCC"
  }
  
  return(result_table)
}

# helper function to compute host colors
compute_host_colors <- function(host_ids, host_table, host_color_scheme, taxa_level, abundance_summary, selected_bin_id) {
  colors <- rep("#CCCCCC", length(host_ids))  # default light gray
  names(colors) <- host_ids
  
  if (host_color_scheme == "none") {
    # light gray for all
    colors[] <- "#CCCCCC"
  } else if (host_color_scheme == "taxa") {
    # use processed host_table with taxa_color column
    if (!is.null(host_table) && "taxa_color" %in% names(host_table)) {
      colors[host_ids] <- host_table$taxa_color[match(host_ids, host_table$bin)]
      colors[is.na(colors)] <- "#CCCCCC"  # gray for missing
    }
  } else if (host_color_scheme == "abundance") {
    # color by mean abundance
    if (!is.null(abundance_summary) && "bin" %in% names(abundance_summary) && "mean_abundance" %in% names(abundance_summary)) {
      abundance_values <- abundance_summary$mean_abundance[match(host_ids, abundance_summary$bin)]
      # create color scale from low (light blue) to high (dark red)
      if (any(!is.na(abundance_values))) {
        max_abund <- max(abundance_values, na.rm = TRUE)
        min_abund <- min(abundance_values, na.rm = TRUE)
        if (max_abund > min_abund) {
          # normalize to 0-1
          normalized <- (abundance_values - min_abund) / (max_abund - min_abund)
          # create color gradient: light blue -> dark red
          colors[host_ids] <- rgb(
            red = normalized,
            green = 0.2,
            blue = 1 - normalized,
            maxColorValue = 1
          )
        }
        colors[is.na(colors)] <- "#CCCCCC"  # gray for missing abundance
      }
    }
  }
  
  # highlight selected bin in red
  if (!is.null(selected_bin_id) && selected_bin_id %in% host_ids) {
    colors[selected_bin_id] <- "#FF0000"
  }
  
  return(colors)
}

# helper function to compute element colors
compute_element_colors <- function(element_ids, element_color_scheme, abundance_summary, selected_bin_id) {
  colors <- rep("#E0E0E0", length(element_ids))  # default very light gray
  names(colors) <- element_ids
  
  if (element_color_scheme == "none") {
    # very light gray for all
    colors[] <- "#E0E0E0"
  } else if (element_color_scheme == "abundance") {
    # color by mean abundance
    if (!is.null(abundance_summary) && "bin" %in% names(abundance_summary) && "mean_abundance" %in% names(abundance_summary)) {
      abundance_values <- abundance_summary$mean_abundance[match(element_ids, abundance_summary$bin)]
      # create color scale from low (light green) to high (dark green)
      if (any(!is.na(abundance_values))) {
        max_abund <- max(abundance_values, na.rm = TRUE)
        min_abund <- min(abundance_values, na.rm = TRUE)
        if (max_abund > min_abund) {
          # normalize to 0-1
          normalized <- (abundance_values - min_abund) / (max_abund - min_abund)
          # handle NAs
          normalized[is.na(normalized)] <- 0
          # create color gradient: light green -> dark green
          for (i in seq_along(element_ids)) {
            if (!is.na(abundance_values[i])) {
              colors[element_ids[i]] <- rgb(
                red = 0.2,
                green = 0.5 + 0.5 * normalized[i],
                blue = 0.2,
                maxColorValue = 1
              )
            }
          }
        }
      }
    }
  }
  
  # highlight selected bin in red
  if (!is.null(selected_bin_id) && selected_bin_id %in% element_ids) {
    colors[selected_bin_id] <- "#FF0000"
  }
  
  return(colors)
}

# helper function to compute edge colors
compute_edge_colors <- function(edges, edge_color_scheme, assembly, get_bin_adj_support_f, get_bin_adj_associated_f, lib1, lib2, selected_bin_id) {
  colors <- rep("#808080", nrow(edges))  # default gray
  names(colors) <- seq_len(nrow(edges))
  
  if (edge_color_scheme == "none") {
    # gray for all
    colors[] <- "#808080"
  } else if (edge_color_scheme == "percent") {
    # color by mean percent (support / associated * 100)
    support_mat <- get_bin_adj_support_f(assembly)
    associated_mat <- get_bin_adj_associated_f(assembly)
    
    if (!is.null(support_mat) && !is.null(associated_mat)) {
      # get library columns
      lib_cols <- names(support_mat)[!names(support_mat) %in% c("bin_src", "bin_tgt")]
      
      for (i in seq_len(nrow(edges))) {
        bin_src <- edges$from[i]
        bin_tgt <- edges$to[i]
        
        # find matching rows
        support_row <- support_mat[support_mat$bin_src == bin_src & support_mat$bin_tgt == bin_tgt, ]
        associated_row <- associated_mat[associated_mat$bin_src == bin_src & associated_mat$bin_tgt == bin_tgt, ]
        
        if (nrow(support_row) > 0 && nrow(associated_row) > 0) {
          # compute mean percent across all libraries
          percents <- numeric(length(lib_cols))
          for (j in seq_along(lib_cols)) {
            lib_col <- lib_cols[j]
            if (lib_col %in% names(support_row) && lib_col %in% names(associated_row)) {
              support_val <- as.numeric(support_row[[lib_col]][1])
              associated_val <- as.numeric(associated_row[[lib_col]][1])
              if (associated_val > 0) {
                percents[j] <- (support_val / associated_val) * 100
              }
            }
          }
          mean_percent <- mean(percents, na.rm = TRUE)
          
          # color scale: low (light blue) to high (dark red)
          # assume percent range 0-100
          normalized <- pmin(pmax(mean_percent / 100, 0), 1)
          colors[i] <- rgb(
            red = normalized,
            green = 0.2,
            blue = 1 - normalized,
            maxColorValue = 1
          )
        }
      }
    }
  } else if (edge_color_scheme == "percent_delta") {
    # color by percent delta (lib2 - lib1)
    if (is.null(lib1) || is.null(lib2) || lib1 == "" || lib2 == "") {
      colors[] <- "#808080"
    } else {
      # extract library IDs from full names (e.g., "EBC_pre" -> "pre")
      lib1_id <- gsub(paste0("^", assembly, "_"), "", lib1)
      lib2_id <- gsub(paste0("^", assembly, "_"), "", lib2)
      
      support_mat <- get_bin_adj_support_f(assembly)
      associated_mat <- get_bin_adj_associated_f(assembly)
      
      if (!is.null(support_mat) && !is.null(associated_mat)) {
        for (i in seq_len(nrow(edges))) {
          bin_src <- edges$from[i]
          bin_tgt <- edges$to[i]
          
          # find matching rows
          support_row <- support_mat[support_mat$bin_src == bin_src & support_mat$bin_tgt == bin_tgt, ]
          associated_row <- associated_mat[associated_mat$bin_src == bin_src & associated_mat$bin_tgt == bin_tgt, ]
          
          if (nrow(support_row) > 0 && nrow(associated_row) > 0) {
            # compute percent for lib1 and lib2
            percent1 <- 0
            percent2 <- 0
            
            if (lib1_id %in% names(support_row) && lib1_id %in% names(associated_row)) {
              support_val1 <- as.numeric(support_row[[lib1_id]][1])
              associated_val1 <- as.numeric(associated_row[[lib1_id]][1])
              if (associated_val1 > 0) {
                percent1 <- (support_val1 / associated_val1) * 100
              }
            }
            
            if (lib2_id %in% names(support_row) && lib2_id %in% names(associated_row)) {
              support_val2 <- as.numeric(support_row[[lib2_id]][1])
              associated_val2 <- as.numeric(associated_row[[lib2_id]][1])
              if (associated_val2 > 0) {
                percent2 <- (support_val2 / associated_val2) * 100
              }
            }
            
            delta <- percent2 - percent1
            
            # color scale: negative (blue) to positive (red), zero (gray)
            if (delta < 0) {
              # blue for negative
              normalized <- pmin(abs(delta) / 100, 1)
              colors[i] <- rgb(
                red = 0.2,
                green = 0.2,
                blue = 0.5 + 0.5 * normalized,
                maxColorValue = 1
              )
            } else if (delta > 0) {
              # red for positive
              normalized <- pmin(delta / 100, 1)
              colors[i] <- rgb(
                red = 0.5 + 0.5 * normalized,
                green = 0.2,
                blue = 0.2,
                maxColorValue = 1
              )
            } else {
              # gray for zero
              colors[i] <- "#808080"
            }
          }
        }
      }
    }
  }
  
  # highlight edges connected to selected bin in red
  if (!is.null(selected_bin_id)) {
    for (i in seq_len(nrow(edges))) {
      if (edges$from[i] == selected_bin_id || edges$to[i] == selected_bin_id) {
        colors[i] <- "#FF0000"
      }
    }
  }
  
  return(colors)
}

# function to render the graph plot
render_associations_graph <- function(assembly, selected_bin_id, get_host_table_f, get_bin_adj_total_f,
                                      get_abundance_summary_f, get_bin_adj_support_f, get_bin_adj_associated_f,
                                      get_bin_segment_table_f,
                                      host_color_scheme, taxa_level, element_color_scheme, edge_color_scheme, lib1, lib2, show_host_labels) {
  if (is.null(assembly)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No assembly selected", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # check cache for graph data
  graph_cache_key <- paste0("associations.graph_data.", assembly)
  graph_data <- cache_get_if_exists(graph_cache_key, NULL)
  
  # load graph data if not cached
  if (is.null(graph_data)) {
    graph_data <- load_graph_data(assembly, get_host_table_f, get_bin_adj_total_f)
    if (is.null(graph_data) || nrow(graph_data$nodes) == 0) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No graph data available", size = 5) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
      return(plotly::ggplotly(p))
    }
    # cache the graph data
    cache_set(graph_cache_key, graph_data)
  }
  
  # check cache for layout data
  layout_cache_key <- paste0("associations.graph_layout.", assembly)
  layout_data <- cache_get_if_exists(layout_cache_key, NULL)
  
  # compute layout if not cached
  if (is.null(layout_data)) {
    layout_data <- compute_graph_layout(graph_data)
    if (is.null(layout_data)) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No edges found", size = 5) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
      return(plotly::ggplotly(p))
    }
    # cache the layout data
    cache_set(layout_cache_key, layout_data)
  }
  
  nodes <- layout_data$nodes
  edges <- layout_data$edges
  host_table <- layout_data$host_table
  
  # process host_table for taxa coloring if needed
  if (host_color_scheme == "taxa") {
    host_table <- process_host_table_taxa(host_table, taxa_level)
  }
  
  # compute bin lengths from bin segment table for elements
  bin_lengths <- NULL
  bin_segment_table <- get_bin_segment_table_f(assembly)
  if (!is.null(bin_segment_table) && "bin" %in% names(bin_segment_table) && "length" %in% names(bin_segment_table)) {
    # sum lengths by bin
    bin_lengths <- aggregate(length ~ bin, data = bin_segment_table, FUN = sum)
    bin_lengths <- setNames(bin_lengths$length, bin_lengths$bin)
  }
  
  # create hover text for nodes
  nodes$hover <- sapply(seq_len(nrow(nodes)), function(i) {
    create_node_hover(nodes$id[i], nodes$type[i], host_table, bin_lengths)
  })
  
  # get abundance summary for coloring
  abundance_summary <- get_abundance_summary_f(assembly)
  
  # compute node colors
  host_ids <- nodes$id[nodes$type == "host"]
  element_ids <- nodes$id[nodes$type == "element"]
  
  host_colors <- compute_host_colors(host_ids, host_table, host_color_scheme, taxa_level, abundance_summary, selected_bin_id)
  element_colors <- compute_element_colors(element_ids, element_color_scheme, abundance_summary, selected_bin_id)
  
  # assign colors to nodes
  nodes$color <- rep("", nrow(nodes))
  nodes$color[nodes$type == "host"] <- host_colors[host_ids]
  nodes$color[nodes$type == "element"] <- element_colors[element_ids]
  nodes$size <- ifelse(nodes$type == "host", 20, 16)  # increased sizes
  
  # compute edge colors
  edge_colors <- compute_edge_colors(edges, edge_color_scheme, assembly, get_bin_adj_support_f, get_bin_adj_associated_f, lib1, lib2, selected_bin_id)
  
  # create plotly plot
  p <- plotly::plot_ly(
    type = "scatter",
    mode = "markers+text",
    source = "associations_graph_plot"
  )
  
  # add edges as lines
  for (i in seq_len(nrow(edges))) {
    edge_width <- ifelse(edges$from[i] == selected_bin_id || edges$to[i] == selected_bin_id, 4, 2)  # thicker lines
    
    p <- p %>% plotly::add_trace(
      x = c(edges$x0[i], edges$x1[i]),
      y = c(edges$y0[i], edges$y1[i]),
      type = "scatter",
      mode = "lines",
      line = list(color = edge_colors[i], width = edge_width),
      showlegend = FALSE,
      hoverinfo = "skip"
    )
  }
  
  # add host nodes (circles, optionally with text)
  host_nodes <- nodes[nodes$type == "host", ]
  if (nrow(host_nodes) > 0) {
    trace_args <- list(
      x = host_nodes$x,
      y = host_nodes$y,
      type = "scatter",
      mode = if (show_host_labels) "markers+text" else "markers",
      marker = list(
        size = host_nodes$size,
        color = host_nodes$color,
        symbol = "circle",
        line = list(color = "black", width = 1)
      ),
      hovertext = host_nodes$hover,
      hoverinfo = "text",
      key = host_nodes$id,
      name = "Hosts",
      showlegend = FALSE
    )
    if (show_host_labels) {
      trace_args$text <- host_nodes$id
      trace_args$textposition <- "middle center"
      trace_args$textfont <- list(size = 10, color = "black")
    }
    p <- do.call(plotly::add_trace, c(list(p = p), trace_args))
  }
  
  # add element nodes (triangles, no labels)
  element_nodes <- nodes[nodes$type == "element", ]
  if (nrow(element_nodes) > 0) {
    p <- p %>% plotly::add_trace(
      x = element_nodes$x,
      y = element_nodes$y,
      type = "scatter",
      mode = "markers",
      marker = list(
        size = element_nodes$size,
        color = element_nodes$color,
        symbol = "triangle-up",
        line = list(color = "black", width = 1)
      ),
      hovertext = element_nodes$hover,
      hoverinfo = "text",
      key = element_nodes$id,
      name = "Elements",
      showlegend = FALSE
    )
  }
  
  # update layout
  p <- p %>% plotly::layout(
    title = "Host-Element Association Graph",
    xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    hovermode = "closest",
    showlegend = TRUE
  )
  
  return(p)
}

# function to render the legend plot
render_associations_legend <- function(assembly, get_host_table_f, get_abundance_summary_f, get_bin_adj_total_f,
                                        host_color_scheme, taxa_level, edge_color_scheme) {
  # get host_ids from graph data
  host_ids <- NULL
  graph_data <- load_graph_data(assembly, get_host_table_f, get_bin_adj_total_f)
  if (!is.null(graph_data) && nrow(graph_data$nodes) > 0) {
    host_ids <- graph_data$nodes$id[graph_data$nodes$type == "host"]
  }
  # create empty plot for legend
  p <- plotly::plot_ly(type = "scatter", mode = "markers")
  
  annotations <- list()
  shapes <- list()
  y_pos <- 0.98
  y_step <- 0.025
  
  # host color legend
  if (host_color_scheme == "taxa") {
    host_table <- get_host_table_f(assembly)
    if (!is.null(host_table)) {
      # process host_table to get sorted table with taxa_color column
      host_table <- process_host_table_taxa(host_table, taxa_level)
      
      if (!is.null(host_table) && "taxa_color" %in% names(host_table) && taxa_level %in% names(host_table)) {
        # get unique taxa values in sorted order (already sorted by process_host_table_taxa)
        taxa_values_vec <- as.character(host_table[[taxa_level]])
        unique_taxa <- unique(taxa_values_vec)
        
        if (length(unique_taxa) > 0) {
          # get colors from the processed table (first occurrence of each unique taxa)
          taxa_colors <- host_table$taxa_color[match(unique_taxa, taxa_values_vec)]
          names(taxa_colors) <- unique_taxa
          
          # count genomes per taxa
          taxa_counts <- as.numeric(table(taxa_values_vec))
          names(taxa_counts) <- names(table(taxa_values_vec))
          
          # add title
          annotations <- c(annotations, list(
            list(
              x = 0.05, y = y_pos,
              xref = "paper", yref = "paper",
              text = paste0("<b>Host Color (", taxa_level, ")</b>"),
              showarrow = FALSE,
              xanchor = "left",
              font = list(size = 11)
            )
          ))
          y_pos <- y_pos - y_step * 1.5
          
          # add taxa entries with colored squares
          for (i in seq_along(unique_taxa)) {
            taxa_val <- unique_taxa[i]
            if (y_pos < 0.05) break  # stop if running out of space
            
            # get count for this taxa
            count <- if (taxa_val %in% names(taxa_counts)) taxa_counts[taxa_val] else 0
            
            # add colored square shape
            shapes <- c(shapes, list(
              list(
                type = "rect",
                xref = "paper", yref = "paper",
                x0 = 0.05, y0 = y_pos - 0.008,
                x1 = 0.08, y1 = y_pos + 0.008,
                fillcolor = taxa_colors[taxa_val],
                line = list(color = "black", width = 0.5)
              )
            ))
            
            # add text with black font, including count in parentheses
            taxa_label <- paste0(substr(taxa_val, 1, 25), " (", count, ")")
            annotations <- c(annotations, list(
              list(
                x = 0.09, y = y_pos,
                xref = "paper", yref = "paper",
                text = taxa_label,
                showarrow = FALSE,
                xanchor = "left",
                yanchor = "middle",
                font = list(size = 9, color = "black")
              )
            ))
            y_pos <- y_pos - y_step
          }
          y_pos <- y_pos - y_step
        }
      }
    }
  } else if (host_color_scheme == "abundance") {
    abundance_summary <- get_abundance_summary_f(assembly)
    if (!is.null(abundance_summary) && "mean_abundance" %in% names(abundance_summary) && !is.null(host_ids)) {
      host_abundances <- abundance_summary$mean_abundance[abundance_summary$bin %in% host_ids]
      if (length(host_abundances) > 0 && any(!is.na(host_abundances))) {
        min_abund <- min(host_abundances, na.rm = TRUE)
        max_abund <- max(host_abundances, na.rm = TRUE)
        annotations <- c(annotations, list(
          list(
            x = 0.05, y = y_pos,
            xref = "paper", yref = "paper",
            text = "<b>Host Color (Abundance)</b>",
            showarrow = FALSE,
            xanchor = "left",
            font = list(size = 11)
          ),
          list(
            x = 0.05, y = y_pos - y_step,
            xref = "paper", yref = "paper",
            text = sprintf("Low: %.2f%%", min_abund * 100),
            showarrow = FALSE,
            xanchor = "left",
            font = list(size = 9, color = "rgb(0, 51, 204)")
          ),
          list(
            x = 0.05, y = y_pos - y_step * 2,
            xref = "paper", yref = "paper",
            text = sprintf("High: %.2f%%", max_abund * 100),
            showarrow = FALSE,
            xanchor = "left",
            font = list(size = 9, color = "rgb(204, 51, 0)")
          )
        ))
        y_pos <- y_pos - y_step * 3.5
      }
    }
  }
  
  # edge color legend
  if (edge_color_scheme == "percent") {
    annotations <- c(annotations, list(
      list(
        x = 0.05, y = y_pos,
        xref = "paper", yref = "paper",
        text = "<b>Edge Color (Percent)</b>",
        showarrow = FALSE,
        xanchor = "left",
        font = list(size = 11)
      ),
      list(
        x = 0.05, y = y_pos - y_step,
        xref = "paper", yref = "paper",
        text = "Low: 0%",
        showarrow = FALSE,
        xanchor = "left",
        font = list(size = 9, color = "rgb(0, 51, 204)")
      ),
      list(
        x = 0.05, y = y_pos - y_step * 2,
        xref = "paper", yref = "paper",
        text = "High: 100%",
        showarrow = FALSE,
        xanchor = "left",
        font = list(size = 9, color = "rgb(204, 51, 0)")
      )
    ))
  } else if (edge_color_scheme == "percent_delta") {
    annotations <- c(annotations, list(
      list(
        x = 0.05, y = y_pos,
        xref = "paper", yref = "paper",
        text = "<b>Edge Color (Delta)</b>",
        showarrow = FALSE,
        xanchor = "left",
        font = list(size = 11)
      ),
      list(
        x = 0.05, y = y_pos - y_step,
        xref = "paper", yref = "paper",
        text = "Negative (blue)",
        showarrow = FALSE,
        xanchor = "left",
        font = list(size = 9, color = "rgb(51, 51, 204)")
      ),
      list(
        x = 0.05, y = y_pos - y_step * 2,
        xref = "paper", yref = "paper",
        text = "Positive (red)",
        showarrow = FALSE,
        xanchor = "left",
        font = list(size = 9, color = "rgb(204, 51, 51)")
      )
    ))
  }
  
  # layout
  p <- p %>% plotly::layout(
    xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE, range = c(0, 1)),
    yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE, range = c(0, 1)),
    annotations = annotations,
    shapes = shapes,
    margin = list(l = 0, r = 0, t = 0, b = 0),
    plot_bgcolor = "rgba(0,0,0,0)",
    paper_bgcolor = "rgba(0,0,0,0)"
  )
  
  return(p)
}

