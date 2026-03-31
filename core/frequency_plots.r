# core/frequency_plots.r
# reusable frequency plotting functionality for rearrangements and variants

library(shiny)
library(plotly)
library(ggplot2)
library(scales)
library(DT)

# create UI for frequency plot controls
# returns a fluidRow with controls and plot output
create_frequency_plot_ui <- function(prefix, library_ids, plot_output_height = "400px") {
  
  # dynamic input IDs
  plot_type_id        <- paste0(prefix, "PlotType")
  plot_value_id       <- paste0(prefix, "PlotValue")
  x_lib_id            <- paste0(prefix, "XLib")
  y_lib_id            <- paste0(prefix, "YLib")
  jitter_id           <- paste0(prefix, "Jitter")
  matrix_color_by_id  <- paste0(prefix, "MatrixColorBy")
  plot_output_id      <- paste0(prefix, "FrequencyPlot")

  # cache keys
  cache_prefix <- paste0(prefix, "_frequency_plot_")
  default_plot_type       <- cache_get_if_exists(paste0(cache_prefix, "type"), "temporal")
  default_plot_value      <- cache_get_if_exists(paste0(cache_prefix, "value"), "frequency")
  default_x_lib           <- cache_get_if_exists(paste0(cache_prefix, "x_lib"), library_ids[1])
  default_y_lib           <- cache_get_if_exists(paste0(cache_prefix, "y_lib"),
                                                 if(length(library_ids) > 1) library_ids[2] else library_ids[1])
  default_jitter          <- cache_get_if_exists(paste0(cache_prefix, "jitter"), FALSE)
  default_matrix_color_by <- cache_get_if_exists(paste0(cache_prefix, "matrix_color_by"), "value")

  fluidRow(
    column(12,
      h4("Frequency Plot"),
      fluidRow(
        column(3,
          selectInput(plot_type_id, "Plot Type:",
                     choices = c("Temporal" = "temporal", "Scatter" = "scatter", "Matrix" = "matrix"),
                     selected = default_plot_type, width = "100%")
        ),
        column(3,
          selectInput(plot_value_id, "Plot Value:",
                     choices = c("Frequency" = "frequency",
                                "Support" = "support",
                                "Coverage" = "coverage"),
                     selected = default_plot_value, width = "100%")
        ),
        # x and y library selectors - conditionally shown for scatter
        column(2,
          conditionalPanel(
            condition = sprintf("input.%s == 'scatter'", plot_type_id),
            selectInput(x_lib_id, "X Library:",
                       choices = setNames(library_ids, library_ids),
                       selected = default_x_lib, width = "100%")
          )
        ),
        column(2,
          conditionalPanel(
            condition = sprintf("input.%s == 'scatter'", plot_type_id),
            selectInput(y_lib_id, "Y Library:",
                       choices = setNames(library_ids, library_ids),
                       selected = default_y_lib, width = "100%")
          )
        ),
        column(2,
          checkboxInput(jitter_id, "Jitter",
                       value = default_jitter, width = "100%")
        )
      ),
      conditionalPanel(
        condition = sprintf("input.%s == 'matrix'", plot_type_id),
        fluidRow(
          column(3,
            selectInput(matrix_color_by_id, "Color By:",
                       choices = c("Value" = "value", "Individuals" = "individuals"),
                       selected = default_matrix_color_by, width = "100%")
          )
        )
      ),
      tags$div(
        style = paste0(
          "resize: vertical; overflow: hidden; ",
          "min-height: 150px; height: ", plot_output_height, ";"
        ),
        plotly::plotlyOutput(plot_output_id, height = "100%")
      ),
      tags$script(HTML(sprintf("
        (function() {
          var plotId = '%s';
          var observerAttached = false;
          document.addEventListener('shiny:value', function(e) {
            if (e.target.id !== plotId) return;
            if (observerAttached) return;
            var checks = 0;
            function waitForPlotly() {
              var plotEl = document.getElementById(plotId);
              if (!plotEl || !plotEl._fullLayout) {
                if (++checks < 30) requestAnimationFrame(waitForPlotly);
                return;
              }
              var wrap = plotEl.parentElement;
              new ResizeObserver(function() {
                var el = document.getElementById(plotId);
                if (el && el._fullLayout && window.Plotly) Plotly.Plots.resize(el);
              }).observe(wrap);
              observerAttached = true;
            }
            requestAnimationFrame(waitForPlotly);
          });
        })();
      ", plot_output_id)))
    )
  )
}

# create frequency plot output renderer
# call this to create the output assignment
create_frequency_plot_output <- function(prefix, data_reactive, library_ids, 
                                        get_selected_items_func = NULL,
                                        empty_message = "No data available") {
  
    # dynamic input IDs and cache keys
  plot_type_id <- paste0(prefix, "PlotType")
  plot_value_id <- paste0(prefix, "PlotValue") 
  x_lib_id <- paste0(prefix, "XLib")
  y_lib_id <- paste0(prefix, "YLib")
  jitter_id <- paste0(prefix, "Jitter")
  plot_output_id <- paste0(prefix, "FrequencyPlot")

  # main plot renderer
  output[[plot_output_id]] <- plotly::renderPlotly({
    
    # get current data
    raw_data <- data_reactive()
    
    # get current settings with defaults
    plot_type <- input[[plot_type_id]] %||% "temporal"
    plot_value <- input[[plot_value_id]] %||% "frequency"
    x_lib <- input[[x_lib_id]] %||% library_ids[1]
    y_lib <- input[[y_lib_id]] %||% (if(length(library_ids) > 1) library_ids[2] else library_ids[1])
    jitter_enabled <- input[[jitter_id]] %||% FALSE
    
    # get selected items if function provided
    selected_items <- if (!is.null(get_selected_items_func)) get_selected_items_func() else NULL
    
    # render the plot
    render_frequency_plot_internal(raw_data, plot_type, plot_value, x_lib, y_lib, 
                                  jitter_enabled, selected_items, library_ids, empty_message)
  })
}

# create frequency plot observers for caching settings
# call this to set up the caching observers
create_frequency_plot_observers <- function(prefix) {
  
  # dynamic input IDs and cache keys
  plot_type_id <- paste0(prefix, "PlotType")
  plot_value_id <- paste0(prefix, "PlotValue") 
  x_lib_id <- paste0(prefix, "XLib")
  y_lib_id <- paste0(prefix, "YLib")
  jitter_id <- paste0(prefix, "Jitter")
  cache_prefix <- paste0(prefix, "_frequency_plot_")
  
  # observers for caching settings
  observeEvent(input[[plot_type_id]], {
    cache_set(paste0(cache_prefix, "type"), input[[plot_type_id]])
  })
  
  observeEvent(input[[plot_value_id]], {
    cache_set(paste0(cache_prefix, "value"), input[[plot_value_id]])
  })
  
  observeEvent(input[[x_lib_id]], {
    cache_set(paste0(cache_prefix, "x_lib"), input[[x_lib_id]])
  })
  
  observeEvent(input[[y_lib_id]], {
    cache_set(paste0(cache_prefix, "y_lib"), input[[y_lib_id]])
  })
  
  observeEvent(input[[jitter_id]], {
    cache_set(paste0(cache_prefix, "jitter"), input[[jitter_id]])
  })
}

# internal plot rendering function
# has_data: boolean flag indicating if we have data to plot
# items_df should have: id, label, color columns (can be NULL if has_data = FALSE)
# support_matrix: items × libraries matrix (can be NULL if has_data = FALSE)
# coverage_matrix: items × libraries matrix (can be NULL if has_data = FALSE)
# max_items: if not NULL, cap to top-N by total_support before plotting
render_frequency_plot_internal <- function(has_data, items_df, support_matrix, coverage_matrix, 
                                         plot_type, plot_value, x_lib, y_lib, 
                                         jitter_enabled, selected_items = NULL, 
                                         library_ids, empty_message = "No data available",
                                         max_items = NULL, sort_by = "id", col_sort_by = "id",
                                         sample_map = NULL, matrix_color_by = "value") {
  # early return for no data
  if (!has_data || is.null(items_df) || nrow(items_df) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = empty_message, size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # light data validation
  if (is.null(support_matrix) || is.null(coverage_matrix)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Invalid data structure", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # validate required columns in items_df
  required_cols <- c("id", "label", "color")
  missing_cols <- required_cols[!required_cols %in% names(items_df)]
  if (length(missing_cols) > 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = paste("Missing columns:", paste(missing_cols, collapse = ", ")), size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # reorder matrices to match configured library order
  lib_indices <- match(library_ids, colnames(support_matrix))
  support_matrix <- support_matrix[, lib_indices, drop = FALSE]
  coverage_matrix <- coverage_matrix[, lib_indices, drop = FALSE]
  
  # cap to top-N items by total_support
  if (!is.null(max_items) && nrow(items_df) > max_items && "total_support" %in% names(items_df)) {
    order_idx <- order(items_df$total_support, decreasing = TRUE)
    keep_idx <- order_idx[seq_len(min(max_items, length(order_idx)))]
    items_df <- items_df[keep_idx, ]
    support_matrix <- support_matrix[keep_idx, , drop = FALSE]
    coverage_matrix <- coverage_matrix[keep_idx, , drop = FALSE]
  }
  
  # calculate data matrix based on plot value
  if (plot_value == "support") {
    data_matrix <- support_matrix
    value_label <- "Support"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else if (plot_value == "coverage") {
    data_matrix <- coverage_matrix
    value_label <- "Coverage"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else { # frequency
    data_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
    value_label <- "Frequency"
    y_limits <- c(0, 1)
    y_format <- scales::percent_format()
  }
  
  # create plot based on plot type
  if (plot_type == "scatter") {
    return(create_scatter_plot(data_matrix, items_df, library_ids, x_lib, y_lib, 
                              value_label, y_limits, y_format, jitter_enabled, 
                              selected_items, plot_value))
  } else if (plot_type == "matrix") {
    return(create_matrix_plot(data_matrix, items_df, library_ids,
                              support_matrix, coverage_matrix,
                              value_label, selected_items, plot_value, sort_by, col_sort_by,
                              sample_map = sample_map, matrix_color_by = matrix_color_by))
  } else {
    return(create_temporal_plot(data_matrix, items_df, library_ids, 
                               value_label, y_limits, y_format, jitter_enabled, 
                               selected_items, plot_value))
  }
}

# create scatter plot (x vs y libraries)
create_scatter_plot <- function(data_matrix, items_df, library_ids, x_lib, y_lib, 
                               value_label, y_limits, y_format, jitter_enabled, 
                               selected_items, plot_value) {
  
  # get library indices
  x_lib_idx <- match(x_lib, library_ids)
  y_lib_idx <- match(y_lib, library_ids)
  
  if (is.na(x_lib_idx) || is.na(y_lib_idx)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Invalid library selection", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # create scatter plot data
  plot_data <- data.frame(
    item_id = items_df$id,
    item_type = items_df$type,
    label = items_df$label,
    color = items_df$color,
    x_value = data_matrix[, x_lib_idx],
    y_value = data_matrix[, y_lib_idx],
    stringsAsFactors = FALSE
  )
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$x_value) & !is.na(plot_data$x_value) &
                        is.finite(plot_data$y_value) & !is.na(plot_data$y_value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # apply bounds based on plot value type (always for frequency)
  if (plot_value == "frequency") {
    # frequency values must stay within [0, 1]
    plot_data$x_value <- pmax(0, pmin(1, plot_data$x_value))
    plot_data$y_value <- pmax(0, pmin(1, plot_data$y_value))
  } else {
    # count values (support/coverage) must not go negative
    plot_data$x_value <- pmax(0, plot_data$x_value)
    plot_data$y_value <- pmax(0, plot_data$y_value)
  }
  
  # check for selected items
  plot_data$is_selected <- FALSE
  if (!is.null(selected_items) && nrow(selected_items) > 0) {
    plot_data$is_selected <- plot_data$item_id %in% selected_items$id
  }
  
  # create hover text
  # format values for hover
  if (plot_value == "frequency") {
    x_text <- paste0(round(plot_data$x_value * 100, 1), "%")
    y_text <- paste0(round(plot_data$y_value * 100, 1), "%")
  } else {
    x_text <- as.character(plot_data$x_value)
    y_text <- as.character(plot_data$y_value)
  }
  
  plot_data$hover_text <- paste0(
    "ID: ", plot_data$item_id, "<br>",
    "Type: ", plot_data$item_type, "<br>",
    "Position: ", plot_data$label, "<br>",
    x_lib, " ", value_label, ": ", x_text, "<br>",
    y_lib, " ", value_label, ": ", y_text, "<br>"
  )
  
  # add jitter if enabled (5% of plot range) - after hover text
  if (jitter_enabled) {
    x_range <- max(plot_data$x_value) - min(plot_data$x_value)
    y_range <- max(plot_data$y_value) - min(plot_data$y_value)
    x_jitter <- x_range * 0.05
    y_jitter <- y_range * 0.05
    
    # apply jitter
    plot_data$x_value <- plot_data$x_value + runif(nrow(plot_data), -x_jitter, x_jitter)
    plot_data$y_value <- plot_data$y_value + runif(nrow(plot_data), -y_jitter, y_jitter)
    
    # reapply bounds after jitter
    if (plot_value == "frequency") {
      plot_data$x_value <- pmax(0, pmin(1, plot_data$x_value))
      plot_data$y_value <- pmax(0, pmin(1, plot_data$y_value))
    } else {
      plot_data$x_value <- pmax(0, plot_data$x_value)
      plot_data$y_value <- pmax(0, plot_data$y_value)
    }
  }
  
  # create plot
  non_selected_data <- plot_data[!plot_data$is_selected, ]
  selected_data <- plot_data[plot_data$is_selected, ]
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_value, y = y_value, color = color, text = hover_text, key = item_id))
  
  # draw non-selected items first
  if (nrow(non_selected_data) > 0) {
    p <- p + ggplot2::geom_point(data = non_selected_data, size = 2, alpha = 0.8)
  }
  
  # draw selected items on top
  if (nrow(selected_data) > 0) {
    p <- p + ggplot2::geom_point(data = selected_data, size = 4, alpha = 1, 
                                stroke = 1.5, shape = 21, fill = "white")
  }
  
  p <- p + ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = paste(x_lib, value_label),
      y = paste(y_lib, value_label)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  
  # add axis formatting
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_x_continuous(limits = y_limits, labels = y_format) +
             ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_x_continuous(labels = y_format) +
             ggplot2::scale_y_continuous(labels = y_format)
  }
  
  plotly::ggplotly(p, tooltip = "text", source = "scatter_plot") %>%
    plotly::layout(hovermode = "closest", showlegend = FALSE)
}

# create temporal plot (libraries on x-axis, connected lines)
create_temporal_plot <- function(data_matrix, items_df, library_ids, 
                                value_label, y_limits, y_format, jitter_enabled, 
                                selected_items, plot_value) {
  
  # convert to long format vectorized (column-major: items repeat per library)
  n_items <- nrow(data_matrix)
  n_libs  <- ncol(data_matrix)
  plot_data <- data.frame(
    item_id   = rep(items_df$id,    times = n_libs),
    item_type = rep(items_df$type,  times = n_libs),
    label     = rep(items_df$label, times = n_libs),
    color     = rep(items_df$color, times = n_libs),
    library   = rep(library_ids,    each  = n_items),
    value     = as.vector(data_matrix),
    stringsAsFactors = FALSE
  )
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$value) & !is.na(plot_data$value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # apply bounds based on plot value type (always for frequency)
  if (plot_value == "frequency") {
    # frequency values must stay within [0, 1]
    plot_data$value <- pmax(0, pmin(1, plot_data$value))
  } else {
    # count values (support/coverage) must not go negative
    plot_data$value <- pmax(0, plot_data$value)
  }
  
  # check for selected items
  plot_data$is_selected <- FALSE
  if (!is.null(selected_items) && nrow(selected_items) > 0) {
    plot_data$is_selected <- plot_data$item_id %in% selected_items$id
  }
  
  # set library order as factor
  plot_data$library <- factor(plot_data$library, levels = library_ids)
  
  # create hover text
  # format value for hover
  if (plot_value == "frequency") {
    value_text <- paste0(round(plot_data$value * 100, 1), "%")
  } else {
    value_text <- as.character(plot_data$value)
  }
  
  plot_data$hover_text <- paste0(
    "ID: ", plot_data$item_id, "<br>",
    "Type: ", plot_data$item_type, "<br>",
    "Position: ", plot_data$label, "<br>",
    "Library: ", plot_data$library, "<br>",
    value_label, ": ", value_text, "<br>"
  )
  
  # add jitter if enabled (5% of y-axis range) - after hover text
  if (jitter_enabled) {
    y_range <- max(plot_data$value) - min(plot_data$value)
    y_jitter <- y_range * 0.05
    
    # apply jitter
    plot_data$value <- plot_data$value + runif(nrow(plot_data), -y_jitter, y_jitter)
    
    # reapply bounds after jitter
    if (plot_value == "frequency") {
      plot_data$value <- pmax(0, pmin(1, plot_data$value))
    } else {
      plot_data$value <- pmax(0, plot_data$value)
    }
  }
  
  # create plot with conditional styling for selected items
  non_selected_data <- plot_data[!plot_data$is_selected, ]
  selected_data <- plot_data[plot_data$is_selected, ]
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = color, 
                                               group = item_id, text = hover_text, key = item_id))
  
  # draw non-selected items first with lower alpha
  if (nrow(non_selected_data) > 0) {
    p <- p + 
      ggplot2::geom_line(data = non_selected_data, alpha = 0.3, size = 0.5) +
      ggplot2::geom_point(data = non_selected_data, size = 2, alpha = 0.8)
  }
  
  # draw selected items on top with higher visibility
  if (nrow(selected_data) > 0) {
    p <- p + 
      ggplot2::geom_line(data = selected_data, alpha = 0.9, size = 2) +
      ggplot2::geom_point(data = selected_data, size = 4, alpha = 1, 
                         stroke = 1.5, shape = 21, fill = "white")
  }
  
  x_label_size <- if (length(library_ids) > 25) 6 else if (length(library_ids) > 15) 8 else NULL

  p <- p + ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Library",
      y = value_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = if (is.null(x_label_size)) ggplot2::element_text(angle = 45, hjust = 1)
                    else ggplot2::element_text(angle = 45, hjust = 1, size = x_label_size),
      legend.position = "none"
    )
  
  # add y-axis formatting
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = y_format)
  }
  
  plotly::ggplotly(p, tooltip = "text", source = "temporal_plot") %>%
    plotly::layout(hovermode = "closest", showlegend = FALSE)
}

# create matrix heatmap plot (interactive plotly)
# variants as rows (hclust-sorted), samples as columns, white-to-red color scale
# NA (support=0, coverage=0 in frequency mode) shown as light gray via plot background
create_matrix_plot <- function(data_matrix, items_df, library_ids,
                               support_matrix, coverage_matrix,
                               value_label, selected_items, plot_value, sort_by = "id",
                               col_sort_by = "id", sample_map = NULL,
                               matrix_color_by = "value") {
  
  n_variants <- nrow(data_matrix)
  n_libs <- length(library_ids)
  
  # for frequency, recompute with NA where coverage is 0 (true missing data)
  if (plot_value == "frequency") {
    value_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, NA_real_)
  } else {
    value_matrix <- data_matrix
  }
  
  # sort rows by user selection
  if (sort_by == "frequency" && "frequency" %in% names(items_df)) {
    row_order <- order(items_df$frequency, decreasing = TRUE, na.last = TRUE)
  } else {
    row_order <- order(items_df$id)
  }
  
  value_matrix_ordered    <- value_matrix[row_order, , drop = FALSE]
  support_ordered         <- support_matrix[row_order, , drop = FALSE]
  coverage_ordered        <- coverage_matrix[row_order, , drop = FALSE]
  items_ordered           <- items_df[row_order, ]

  # cluster columns (samples) by hclust when requested
  if (col_sort_by == "frequency" && n_libs > 1) {
    mat_for_col_clust <- value_matrix_ordered
    mat_for_col_clust[is.na(mat_for_col_clust)] <- 0
    hc_col <- hclust(dist(t(mat_for_col_clust)))
    col_order <- hc_col$order
    library_ids          <- library_ids[col_order]
    value_matrix_ordered <- value_matrix_ordered[, col_order, drop = FALSE]
    support_ordered      <- support_ordered[, col_order, drop = FALSE]
    coverage_ordered     <- coverage_ordered[, col_order, drop = FALSE]
  }

  show_labels <- n_variants <= 10

  # build hover text matrix (always shows frequency/support/coverage value)
  if (plot_value == "frequency") {
    val_str_vec <- paste0(round(as.vector(value_matrix_ordered) * 100, 1), "%")
  } else {
    val_str_vec <- as.character(round(as.vector(value_matrix_ordered)))
  }
  val_str_vec[is.na(as.vector(value_matrix_ordered))] <- "NA"
  hover_matrix <- matrix(
    paste0("Variant: ", rep(items_ordered$id, times = n_libs), "<br>",
           "Sample: ", rep(library_ids, each = n_variants), "<br>",
           value_label, ": ", val_str_vec),
    nrow = n_variants, ncol = n_libs
  )

  # choose z matrix and colorscale:
  # "individuals" mode: orange=1 individual, light blue=2+, gray=no support
  # falls back to value coloring when n_individuals column is absent
  if (matrix_color_by == "individuals" && "n_individuals" %in% names(items_ordered)) {
    # per-cell data.frame: one row per (variant, sample) cell
    cells <- data.frame(
      n_individuals = rep(items_ordered$n_individuals, times = n_libs),
      support       = as.vector(support_ordered),
      coverage      = as.vector(coverage_ordered),
      stringsAsFactors = FALSE
    )
    # z=0: detected, 1 individual (orange)
    # z=1: not detected but sampled (white)
    # z=2: detected, 2+ individuals (blue)
    # NA: no coverage (gray background)
    cells$z <- ifelse(cells$coverage == 0, NA_real_,
                 ifelse(cells$support > 0 & cells$n_individuals == 1, 0L,
                   ifelse(cells$support == 0, 1L, 2L)))

    z_for_plot <- matrix(cells$z, nrow = nrow(items_ordered), ncol = n_libs)
    # sharp transitions at 1/3 and 2/3 produce three discrete bands
    plot_colorscale <- list(
      list(0,    "#FFA500"), list(0.25, "#FFA500"),
      list(0.25, "white"),   list(0.75, "white"),
      list(0.75, "#add8e6"), list(1,    "#add8e6")
    )
    zmin_val      <- 0
    zmax_val      <- 2
    colorbar_opts <- list(title = "Individuals", tickvals = c(0, 1, 2),
                          ticktext = c("1 individual", "not detected", "2+"), len = 0.4)
  } else {
    z_for_plot      <- value_matrix_ordered
    plot_colorscale <- list(list(0, "#add8e6"), list(1, "red"))
    zmin_val        <- 0
    zmax_val        <- if (plot_value == "frequency") 1 else max(value_matrix_ordered, na.rm = TRUE)
    colorbar_opts   <- list(title = value_label)
  }

  # build rect shapes for selected rows (0-based index matches plotly categorical axis)
  shapes_list <- list()
  if (!is.null(selected_items) && nrow(selected_items) > 0) {
    selected_row_indices <- which(items_ordered$id %in% selected_items$id) - 1
    for (row_idx in selected_row_indices) {
      shapes_list[[length(shapes_list) + 1]] <- list(
        type = "rect",
        xref = "x", yref = "y",
        x0 = -0.5, x1 = n_libs - 0.5,
        y0 = row_idx - 0.5, y1 = row_idx + 0.5,
        line = list(color = "black", width = 0.6),
        fillcolor = "rgba(0,0,0,0)"
      )
    }
  }

  # vertical dividers between individuals when sorting by id
  if (col_sort_by == "id" && !is.null(sample_map) &&
      "sample" %in% names(sample_map) && "individual" %in% names(sample_map)) {
    lib_individuals <- sample_map$individual[match(library_ids, sample_map$sample)]
    for (j in seq_len(n_libs - 1)) {
      if (!is.na(lib_individuals[j]) && !is.na(lib_individuals[j + 1]) &&
          lib_individuals[j] != lib_individuals[j + 1]) {
        shapes_list[[length(shapes_list) + 1]] <- list(
          type = "line",
          xref = "x", yref = "paper",
          x0 = j - 0.5, x1 = j - 0.5,
          y0 = 0, y1 = 1,
          line = list(color = "black", width = 1.5)
        )
      }
    }
  }

  p <- plotly::plot_ly(
    z = z_for_plot,
    x = library_ids,
    y = items_ordered$id,
    type = "heatmap",
    colorscale = plot_colorscale,
    zmin = zmin_val,
    zmax = zmax_val,
    text = hover_matrix,
    hovertemplate = "%{text}<extra></extra>",
    colorbar = colorbar_opts,
    source = "matrix_plot"
  ) %>%
  plotly::layout(
    xaxis = list(title = "Sample", tickangle = -90, side = "bottom",
                 tickfont = list(size = if (n_libs > 25) 7 else if (n_libs > 15) 8 else 10)),
    yaxis = list(title = "", showticklabels = show_labels, autorange = "reversed"),
    # light gray background makes NA cells (no support) visually distinct
    paper_bgcolor = "#e8e8e8",
    plot_bgcolor = "#e8e8e8",
    shapes = shapes_list
  )
  
  return(p)
}

# create frequency plot for PDF export
# returns ggplot object (not plotly) with proper sizing for PDF
create_frequency_plot_for_export <- function(has_data, items_df, support_matrix, coverage_matrix,
                                           plot_type, plot_value, x_lib, y_lib, 
                                           jitter_enabled, library_ids, title = "Frequency Plot",
                                           max_items = NULL, col_sort_by = "id", sort_by = "id",
                                           matrix_color_by = "value", sample_map = NULL,
                                           selected_variant_id = NULL) {
  
  # early return for no data - create empty plot
  if (!has_data || is.null(items_df) || nrow(items_df) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
    
    # calculate dimensions based on plot type
    if (plot_type == "scatter") {
      width_inches <- 4  # square plot for scatter
      height_inches <- 4
    } else {
      width_inches <- max(3, length(library_ids) * 0.6)  # temporal: width based on libraries
      height_inches <- 4
    }
    
    return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
  }
  
  # data validation
  if (is.null(support_matrix) || is.null(coverage_matrix)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Invalid data structure", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
    
    width_inches <- if (plot_type == "scatter") 4 else max(3, length(library_ids) * 0.6)
    height_inches <- 4
    
    return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
  }
  
  # validate required columns in items_df
  required_cols <- c("id", "label", "color")
  missing_cols <- required_cols[!required_cols %in% names(items_df)]
  if (length(missing_cols) > 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = paste("Missing columns:", paste(missing_cols, collapse = ", ")), size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
    
    width_inches <- if (plot_type == "scatter") 4 else max(3, length(library_ids) * 0.6)
    height_inches <- 4
    
    return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
  }
  
  # reorder matrices to match configured library order
  lib_indices <- match(library_ids, colnames(support_matrix))
  support_matrix <- support_matrix[, lib_indices, drop = FALSE]
  coverage_matrix <- coverage_matrix[, lib_indices, drop = FALSE]
  
  # cap to top-N items by total_support
  if (!is.null(max_items) && nrow(items_df) > max_items && "total_support" %in% names(items_df)) {
    order_idx <- order(items_df$total_support, decreasing = TRUE)
    keep_idx <- order_idx[seq_len(min(max_items, length(order_idx)))]
    items_df <- items_df[keep_idx, ]
    support_matrix <- support_matrix[keep_idx, , drop = FALSE]
    coverage_matrix <- coverage_matrix[keep_idx, , drop = FALSE]
  }
  
  # calculate data matrix based on plot value
  if (plot_value == "support") {
    data_matrix <- support_matrix
    value_label <- "Support"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else if (plot_value == "coverage") {
    data_matrix <- coverage_matrix
    value_label <- "Coverage"
    y_limits <- NULL
    y_format <- scales::number_format()
  } else { # frequency
    data_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, 0)
    value_label <- "Frequency"
    y_limits <- c(0, 1)
    y_format <- scales::percent_format()
  }
  
  # create plot based on plot type
  if (plot_type == "scatter") {
    return(create_scatter_plot_for_export(data_matrix, items_df, library_ids, x_lib, y_lib, 
                                         value_label, y_limits, y_format, jitter_enabled, 
                                         plot_value, title))
  } else if (plot_type == "matrix") {
    return(create_matrix_plot_for_export(data_matrix, items_df, library_ids,
                                        support_matrix, coverage_matrix,
                                        value_label, plot_value, title,
                                        col_sort_by = col_sort_by, sort_by = sort_by,
                                        matrix_color_by = matrix_color_by,
                                        sample_map = sample_map,
                                        selected_variant_id = selected_variant_id))
  } else {
    return(create_temporal_plot_for_export(data_matrix, items_df, library_ids, 
                                          value_label, y_limits, y_format, jitter_enabled, 
                                          plot_value, title))
  }
}

# create scatter plot for PDF export (returns ggplot)
create_scatter_plot_for_export <- function(data_matrix, items_df, library_ids, x_lib, y_lib, 
                                          value_label, y_limits, y_format, jitter_enabled, 
                                          plot_value, title) {
  
  # get library indices
  x_lib_idx <- match(x_lib, library_ids)
  y_lib_idx <- match(y_lib, library_ids)
  
  if (is.na(x_lib_idx) || is.na(y_lib_idx)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Invalid library selection", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
    return(list(plot = p, width_inches = 4, height_inches = 4))
  }
  
  # create scatter plot data
  plot_data <- data.frame(
    item_id = items_df$id,
    item_type = items_df$type,
    label = items_df$label,
    color = items_df$color,
    x_value = data_matrix[, x_lib_idx],
    y_value = data_matrix[, y_lib_idx],
    stringsAsFactors = FALSE
  )
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$x_value) & !is.na(plot_data$x_value) &
                        is.finite(plot_data$y_value) & !is.na(plot_data$y_value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
    return(list(plot = p, width_inches = 4, height_inches = 4))
  }
  
  # add jitter if enabled (no bounds for PDF - keep clean)
  if (jitter_enabled) {
    x_range <- max(plot_data$x_value) - min(plot_data$x_value)
    y_range <- max(plot_data$y_value) - min(plot_data$y_value)
    x_jitter <- x_range * 0.02  # smaller jitter for PDF
    y_jitter <- y_range * 0.02
    
    plot_data$x_value <- plot_data$x_value + runif(nrow(plot_data), -x_jitter, x_jitter)
    plot_data$y_value <- plot_data$y_value + runif(nrow(plot_data), -y_jitter, y_jitter)
    
    # apply bounds based on plot value type
    if (plot_value == "frequency") {
      plot_data$x_value <- pmax(0, pmin(1, plot_data$x_value))
      plot_data$y_value <- pmax(0, pmin(1, plot_data$y_value))
    } else {
      plot_data$x_value <- pmax(0, plot_data$x_value)
      plot_data$y_value <- pmax(0, plot_data$y_value)
    }
  }
  
  # create plot (no selection highlighting for PDF)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_value, y = y_value, color = color)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = paste(x_lib, value_label),
      y = paste(y_lib, value_label),
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  # add axis formatting
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_x_continuous(limits = y_limits, labels = y_format) +
             ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_x_continuous(labels = y_format) +
             ggplot2::scale_y_continuous(labels = y_format)
  }
  
  return(list(plot = p, width_inches = 4, height_inches = 4))
}

# create temporal plot for PDF export (returns ggplot)
create_temporal_plot_for_export <- function(data_matrix, items_df, library_ids, 
                                           value_label, y_limits, y_format, jitter_enabled, 
                                           plot_value, title) {
  
  # convert to long format vectorized (column-major: items repeat per library)
  n_items <- nrow(data_matrix)
  n_libs  <- ncol(data_matrix)
  plot_data <- data.frame(
    item_id   = rep(items_df$id,    times = n_libs),
    item_type = rep(items_df$type,  times = n_libs),
    label     = rep(items_df$label, times = n_libs),
    color     = rep(items_df$color, times = n_libs),
    library   = rep(library_ids,    each  = n_items),
    value     = as.vector(data_matrix),
    stringsAsFactors = FALSE
  )
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$value) & !is.na(plot_data$value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
    
    width_inches <- max(3, length(library_ids) * 0.6)
    return(list(plot = p, width_inches = width_inches, height_inches = 4))
  }
  
  # add jitter if enabled (smaller for PDF)
  if (jitter_enabled) {
    y_range <- max(plot_data$value) - min(plot_data$value)
    y_jitter <- y_range * 0.02  # smaller jitter for PDF
    
    plot_data$value <- plot_data$value + runif(nrow(plot_data), -y_jitter, y_jitter)
    
    # apply bounds based on plot value type
    if (plot_value == "frequency") {
      plot_data$value <- pmax(0, pmin(1, plot_data$value))
    } else {
      plot_data$value <- pmax(0, plot_data$value)
    }
  }
  
  # set library order as factor
  plot_data$library <- factor(plot_data$library, levels = library_ids)
  
  x_label_size <- if (length(library_ids) > 25) 5 else if (length(library_ids) > 15) 7 else 8

  # create the plot (no selection highlighting for PDF)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = color, group = item_id)) +
    ggplot2::geom_line(alpha = 0.7, size = 0.8) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Library",
      y = value_label,
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = x_label_size),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  # add y-axis formatting
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = y_format)
  }
  
  # calculate dimensions: width based on library count
  width_inches <- max(3, length(library_ids) * 0.6)
  height_inches <- 4
  
  return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
}

# create matrix heatmap for PDF export (returns ggplot)
create_matrix_plot_for_export <- function(data_matrix, items_df, library_ids,
                                          support_matrix, coverage_matrix,
                                          value_label, plot_value, title,
                                          col_sort_by = "id", sort_by = "id",
                                          matrix_color_by = "value", sample_map = NULL,
                                          selected_variant_id = NULL) {
  
  n_variants <- nrow(data_matrix)
  
  # for frequency, use NA where coverage is 0
  if (plot_value == "frequency") {
    value_matrix <- ifelse(coverage_matrix > 0, support_matrix / coverage_matrix, NA_real_)
  } else {
    value_matrix <- data_matrix
  }
  
  if (sort_by == "frequency" && "frequency" %in% names(items_df)) {
    row_order <- order(items_df$frequency, decreasing = TRUE, na.last = TRUE)
  } else {
    row_order <- order(items_df$id)
  }
  
  value_matrix_ordered <- value_matrix[row_order, , drop = FALSE]
  support_ordered      <- support_matrix[row_order, , drop = FALSE]
  coverage_ordered     <- coverage_matrix[row_order, , drop = FALSE]
  items_ordered        <- items_df[row_order, ]

  if (col_sort_by == "frequency" && length(library_ids) > 1) {
    mat_for_col_clust <- value_matrix_ordered
    mat_for_col_clust[is.na(mat_for_col_clust)] <- 0
    hc_col <- hclust(dist(t(mat_for_col_clust)))
    col_order            <- hc_col$order
    library_ids          <- library_ids[col_order]
    value_matrix_ordered <- value_matrix_ordered[, col_order, drop = FALSE]
    support_ordered      <- support_ordered[, col_order, drop = FALSE]
    coverage_ordered     <- coverage_ordered[, col_order, drop = FALSE]
  }

  # build per-cell data.frame for fill assignment
  n_ind_col <- if ("n_individuals" %in% names(items_ordered)) items_ordered$n_individuals
               else rep(NA_integer_, n_variants)
  cells <- data.frame(
    row_i         = rep(seq_len(n_variants), times = length(library_ids)),
    col_j         = rep(seq_along(library_ids), each = n_variants),
    n_individuals = rep(n_ind_col, times = length(library_ids)),
    support       = as.vector(support_ordered),
    coverage      = as.vector(coverage_ordered),
    stringsAsFactors = FALSE
  )
  cells$variant_id <- factor(items_ordered$id[cells$row_i], levels = rev(items_ordered$id))
  cells$library    <- factor(library_ids[cells$col_j], levels = library_ids)

  use_individuals <- matrix_color_by == "individuals" && "n_individuals" %in% names(items_ordered)

  if (use_individuals) {
    cells$fill_cat <- ifelse(
      cells$coverage == 0, NA_character_,
      ifelse(cells$support > 0 & cells$n_individuals == 1, "1 individual",
        ifelse(cells$support == 0, "not detected", "2+ individuals")))
    cells$fill_cat <- factor(cells$fill_cat,
                             levels = c("1 individual", "not detected", "2+ individuals"))
  } else {
    cells$fill_val <- as.vector(value_matrix_ordered)
  }

  x_label_size <- 6

  if (use_individuals) {
    p <- ggplot2::ggplot(cells, ggplot2::aes(x = library, y = variant_id, fill = fill_cat)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(
        name = "Individuals",
        values = c("1 individual" = "#FFA500", "not detected" = "white", "2+ individuals" = "#add8e6"),
        na.value = "lightgray", drop = FALSE)
  } else {
    p <- ggplot2::ggplot(cells, ggplot2::aes(x = library, y = variant_id, fill = fill_val)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient(low = "#add8e6", high = "red", na.value = "lightgray",
                                   name = value_label)
  }

  # left margin accounts for variant ID labels
  max_id_nchar  <- max(nchar(as.character(items_ordered$id)), na.rm = TRUE)
  left_margin_in <- 0.1 + max_id_nchar * 0.045

  p <- p +
    ggplot2::labs(x = "Sample", y = "", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_label_size),
      axis.text.y = ggplot2::element_text(size = 6),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = left_margin_in * 72, unit = "pt")
    )

  # vertical dividers between individuals when sorting by id
  if (col_sort_by == "id" && !is.null(sample_map) &&
      "sample" %in% names(sample_map) && "individual" %in% names(sample_map)) {
    lib_individuals <- sample_map$individual[match(library_ids, sample_map$sample)]
    divider_positions <- c()
    for (j in seq_len(length(library_ids) - 1)) {
      if (!is.na(lib_individuals[j]) && !is.na(lib_individuals[j + 1]) &&
          lib_individuals[j] != lib_individuals[j + 1]) {
        divider_positions <- c(divider_positions, j + 0.5)
      }
    }
    if (length(divider_positions) > 0) {
      p <- p + ggplot2::geom_vline(xintercept = divider_positions,
                                   color = "black", linewidth = 0.6)
    }
  }

  # selection rect around the selected variant row
  if (!is.null(selected_variant_id) && selected_variant_id %in% items_ordered$id) {
    y_pos <- which(rev(items_ordered$id) == selected_variant_id)
    p <- p + ggplot2::annotate("rect",
      xmin = 0.5, xmax = length(library_ids) + 0.5,
      ymin = y_pos - 0.5, ymax = y_pos + 0.5,
      fill = NA, color = "black", linewidth = 0.3)
  }

  width_inches  <- max(3, length(library_ids) * 0.2) + left_margin_in
  height_inches <- max(3, n_variants * 0.15 + 1.5)
  
  return(list(plot = p, width_inches = width_inches, height_inches = height_inches))
}

