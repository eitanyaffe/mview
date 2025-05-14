# ---- Dynamic Profile Rendering and Zoom ----

# Helper to build the context (cxt) for plotting
build_profile_context <- function() {
  valid_cids <- intersect(state$contigs, contigs$cid)
  cdf <- contigs[match(valid_cids, contigs$cid), ]
  if (nrow(cdf) == 0) {
    return(NULL)
  }
  cdf$start <- cumsum(c(0, head(cdf$length, -1)))
  cdf$end <- cdf$start + cdf$length
  mapper <- list(
    g2l = function(gcoords) {
      ii <- findInterval(gcoords, c(cdf$start, Inf), left.open = TRUE, all.inside = TRUE)
      if (any(ii <= 0 | ii > nrow(cdf))) stop("Some coordinate indices are out of range")
      data.frame(cid = cdf$cid[ii], coord = gcoords - cdf$start[ii], gcoord = gcoords)
    },
    l2g = function(df) {
      idx <- match(df$cid, cdf$cid)
      if (any(is.na(idx))) {
        stop(sprintf(
          "Some cids are not found in the contig table: %s",
          paste(df$cid[is.na(idx)], collapse = ", ")
        ))
      }
      cdf$start[idx] + df$coord
    },
    cdf = cdf
  )
  if (!is.null(state$zoom)) {
    mapper$xlim <- range(state$zoom)
  } else {
    mapper$xlim <- range(cdf$start, cdf$end)
  }
  list(mapper = mapper, state = state)
}

# Store info_df for each profile for efficient hover
profile_info_dfs <- reactiveValues()

# Track current view for reactivity
current_view <- reactiveVal(NULL)

# UI: dynamically generate plot outputs for all registered profiles
output$profilePlots <- renderUI({
  # Make this reactive to view changes
  current_view()
  profiles <- profiles_get_all()
  if (length(profiles) == 0) {
    return(NULL)
  }
  cat(sprintf("Profiles available in renderUI: %s\n", paste(names(profiles), collapse = ", ")))
  weights <- sapply(profiles, function(p) p$height)
  total <- sum(weights)
  brush_opts <- brushOpts(id = "zoomBrush", direction = "x")

  plot_list <- lapply(names(profiles), function(id) {
    prof <- profiles[[id]]
    plotOutput(
      outputId = paste0("plot_", id),
      height = paste0(300 * prof$height / total, "px"),
      hover = hoverOpts(id = paste0("hover_", id), delay = 50),
      brush = brush_opts,
      dblclick = "plot_dblclick"
    )
  })
  tagList(plot_list)
})

# Server: render plots for all registered profiles and cache info_df for hover
observe({
  cxt <- build_profile_context()
  if (is.null(cxt)) {
    return()
  }
  plots <- plot_profiles(cxt)
  for (id in names(plots)) {
    local({
      myid <- id
      myplot <- plots[[id]]$plot
      myinfo <- plots[[id]]$info_df
      output[[paste0("plot_", myid)]] <- renderPlot({
        myplot
      })
      profile_info_dfs[[myid]] <- myinfo
    })
  }
})

# Zoom logic
observeEvent(input$plot_dblclick, {
  b <- input$zoomBrush
  click_x <- input$plot_dblclick$x
  if (!is.null(b) && click_x >= b$xmin && click_x <= b$xmax) {
    state$zoom <- c(b$xmin, b$xmax)
    addLog(sprintf("Zoom set: %.1f - %.1f", b$xmin, b$xmax))
  } else {
    state$zoom <- NULL
    addLog("Zoom reset")
  }
})

# Hover info (efficient, uses cached info_df)
output$profileHoverInfo <- renderText({
  hover_data <- reactiveValuesToList(input)
  for (id in names(profile_info_dfs)) {
    hname <- paste0("hover_", id)
    if (!is.null(hover_data[[hname]])) {
      h <- hover_data[[hname]]
      gcoord <- as.numeric(h$x)
      info_df <- profile_info_dfs[[id]]
      if (!is.null(info_df) && nrow(info_df) > 0) {
        idx <- which.min(abs(info_df$gcoord - gcoord))
        if (length(idx) == 1 && is.finite(info_df$gcoord[idx])) {
          return(info_df$description[idx])
        }
      }
    }
  }
  NULL
})

# Update current_view when view changes
observeEvent(input$view_id,
  {
    current_view(input$view_id)
  },
  ignoreInit = TRUE
)
