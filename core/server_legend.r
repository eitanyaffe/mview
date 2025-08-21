# core/server_legend.r
# Server-side rendering for the legend tab

# Reactive observer to create all legend renderPlot outputs
observe({
  legends <- state$current_legends
  legend_scale <- if (!is.null(input$legend_scale)) input$legend_scale else 1.0
  
  if (!is.null(legends) && length(legends) > 0) {
    for (profile_id in names(legends)) {
      profile_legends <- legends[[profile_id]]
      
      if (!is.null(profile_legends) && length(profile_legends) > 0) {
        for (i in seq_along(profile_legends)) {
          legend_info <- profile_legends[[i]]
          if (!is.null(legend_info) && !is.null(legend_info$gg)) {
            output_id <- paste0("legend_", profile_id, "_", i)
            
            # calculate scaled dimensions
            height <- if (!is.null(legend_info$height)) legend_info$height * legend_scale else 200 * legend_scale
            width <- if (!is.null(legend_info$width)) legend_info$width * legend_scale else 300 * legend_scale
            
            # Create a local copy for proper scoping in renderPlot
            local({
              local_legend_info <- legend_info
              local_height <- height
              local_width <- width
              output[[output_id]] <- renderPlot({
                local_legend_info$gg
              }, height = local_height, width = local_width)
            })
          }
        }
      }
    }
  }
})

# render the legends content UI (no renderPlot calls here)
output$legendsContent <- renderUI({
  legends <- state$current_legends
  legend_scale <- if (!is.null(input$legend_scale)) input$legend_scale else 1.0
  
  # Check if legends are available
  
  if (is.null(legends) || length(legends) == 0) {
    return(div(
      p("No legends available for current profiles."),
      style = "padding: 20px; text-align: center; color: #666;"
    ))
  }
  
  # create UI elements for each profile's legends
  legend_elements <- list()
  
  for (profile_id in names(legends)) {
    profile_legends <- legends[[profile_id]]
    if (is.null(profile_legends) || length(profile_legends) == 0) {
      next
    }
    
    # add profile header
    legend_elements[[length(legend_elements) + 1]] <- h4(
      paste("Profile:", profile_id), 
      style = "margin-bottom: 15px; margin-top: 20px; color: #333; border-bottom: 1px solid #ddd; padding-bottom: 5px;"
    )
    
    # add each legend for this profile
    for (i in seq_along(profile_legends)) {
      legend_info <- profile_legends[[i]]
      if (is.null(legend_info) || is.null(legend_info$gg)) {
        next
      }
      
      # calculate scaled dimensions
      height <- if (!is.null(legend_info$height)) legend_info$height * legend_scale else 200 * legend_scale
      width <- if (!is.null(legend_info$width)) legend_info$width * legend_scale else 300 * legend_scale
      
      # create output id for this legend
      output_id <- paste0("legend_", profile_id, "_", i)

      # add legend title if provided
      legend_title <- if (!is.null(legend_info$title)) legend_info$title else paste("Legend", i)
      
      # add to UI elements (renderPlot is created in observer above)
      legend_elements[[length(legend_elements) + 1]] <- div(
        h6(legend_title, style = "margin-bottom: 10px; color: #555; font-weight: bold;"),
        plotOutput(output_id, height = paste0(height, "px"), width = paste0(width, "px")),
        style = "margin-bottom: 15px; margin-left: 20px; padding: 10px; border: 1px solid #ddd; border-radius: 4px; background-color: #f9f9f9;"
      )
    }
  }

  if (length(legend_elements) == 0) {
    return(div(
      p("No legends to display."),
      style = "padding: 20px; text-align: center; color: #666;"
    ))
  }
  
  # combine all legend elements
  do.call(tagList, legend_elements)
})
