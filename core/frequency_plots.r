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
  plot_type_id <- paste0(prefix, "PlotType")
  plot_value_id <- paste0(prefix, "PlotValue") 
  x_lib_id <- paste0(prefix, "XLib")
  y_lib_id <- paste0(prefix, "YLib")
  jitter_id <- paste0(prefix, "Jitter")
  plot_output_id <- paste0(prefix, "FrequencyPlot")
  
  # cache keys
  cache_prefix <- paste0(prefix, "_frequency_plot_")
  default_plot_type <- cache_get_if_exists(paste0(cache_prefix, "type"), "temporal")
  default_plot_value <- cache_get_if_exists(paste0(cache_prefix, "value"), "frequency")
  default_x_lib <- cache_get_if_exists(paste0(cache_prefix, "x_lib"), library_ids[1])
  default_y_lib <- cache_get_if_exists(paste0(cache_prefix, "y_lib"), 
                                       if(length(library_ids) > 1) library_ids[2] else library_ids[1])
  default_jitter <- cache_get_if_exists(paste0(cache_prefix, "jitter"), FALSE)
  
  fluidRow(
    column(12,
      h4("Frequency Plot"),
      fluidRow(
        column(3,
          selectInput(plot_type_id, "Plot Type:", 
                     choices = c("Temporal" = "temporal", "Scatter" = "scatter"),
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
      plotly::plotlyOutput(plot_output_id, height = plot_output_height)
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
render_frequency_plot_internal <- function(has_data, items_df, support_matrix, coverage_matrix, 
                                         plot_type, plot_value, x_lib, y_lib, 
                                         jitter_enabled, selected_items = NULL, 
                                         library_ids, empty_message = "No data available") {
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
  
  # add jitter if enabled (5% of plot range)
  if (jitter_enabled) {
    x_range <- max(plot_data$x_value) - min(plot_data$x_value)
    y_range <- max(plot_data$y_value) - min(plot_data$y_value)
    x_jitter <- x_range * 0.05
    y_jitter <- y_range * 0.05
    
    # apply jitter
    plot_data$x_value <- plot_data$x_value + runif(nrow(plot_data), -x_jitter, x_jitter)
    plot_data$y_value <- plot_data$y_value + runif(nrow(plot_data), -y_jitter, y_jitter)
    
    # apply bounds based on plot value type
    if (plot_value == "frequency") {
      # frequency values must stay within [0, 1]
      plot_data$x_value <- pmax(0, pmin(1, plot_data$x_value))
      plot_data$y_value <- pmax(0, pmin(1, plot_data$y_value))
    } else {
      # count values (support/coverage) must not go negative
      plot_data$x_value <- pmax(0, plot_data$x_value)
      plot_data$y_value <- pmax(0, plot_data$y_value)
    }
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
  
  # create plot
  non_selected_data <- plot_data[!plot_data$is_selected, ]
  selected_data <- plot_data[plot_data$is_selected, ]
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_value, y = y_value, color = color, text = hover_text))
  
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
  
  plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(hovermode = "closest", showlegend = FALSE)
}

# create temporal plot (libraries on x-axis, connected lines)
create_temporal_plot <- function(data_matrix, items_df, library_ids, 
                                value_label, y_limits, y_format, jitter_enabled, 
                                selected_items, plot_value) {
  
  # convert to long format for plotting
  plot_data <- data.frame()
  for (i in seq_len(nrow(data_matrix))) {
    for (j in seq_len(ncol(data_matrix))) {
      plot_data <- rbind(plot_data, data.frame(
        item_id = items_df$id[i],
        item_type = items_df$type[i],
        label = items_df$label[i],
        color = items_df$color[i],
        library = library_ids[j],
        value = data_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # remove invalid values
  plot_data <- plot_data[is.finite(plot_data$value) & !is.na(plot_data$value), ]
  
  if (nrow(plot_data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data", size = 5) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }
  
  # add jitter if enabled (5% of y-axis range)
  if (jitter_enabled) {
    y_range <- max(plot_data$value) - min(plot_data$value)
    y_jitter <- y_range * 0.05
    
    # apply jitter
    plot_data$value <- plot_data$value + runif(nrow(plot_data), -y_jitter, y_jitter)
    
    # apply bounds based on plot value type
    if (plot_value == "frequency") {
      # frequency values must stay within [0, 1]
      plot_data$value <- pmax(0, pmin(1, plot_data$value))
    } else {
      # count values (support/coverage) must not go negative
      plot_data$value <- pmax(0, plot_data$value)
    }
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
  
  # create plot with conditional styling for selected items
  non_selected_data <- plot_data[!plot_data$is_selected, ]
  selected_data <- plot_data[plot_data$is_selected, ]
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = library, y = value, color = color, 
                                               group = item_id, text = hover_text))
  
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
  
  p <- p + ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Library",
      y = value_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # add y-axis formatting
  if (!is.null(y_limits)) {
    p <- p + ggplot2::scale_y_continuous(limits = y_limits, labels = y_format)
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = y_format)
  }
  
  plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(hovermode = "closest", showlegend = FALSE)
}

# create frequency plot for PDF export
# returns ggplot object (not plotly) with proper sizing for PDF
create_frequency_plot_for_export <- function(has_data, items_df, support_matrix, coverage_matrix,
                                           plot_type, plot_value, x_lib, y_lib, 
                                           jitter_enabled, library_ids, title = "Frequency Plot") {
  
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
  
  # convert to long format for plotting
  plot_data <- data.frame()
  for (i in seq_len(nrow(data_matrix))) {
    for (j in seq_len(ncol(data_matrix))) {
      plot_data <- rbind(plot_data, data.frame(
        item_id = items_df$id[i],
        item_type = items_df$type[i],
        label = items_df$label[i],
        color = items_df$color[i],
        library = library_ids[j],
        value = data_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
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
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
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

