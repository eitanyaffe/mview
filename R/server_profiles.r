# ---- Dynamic Profile Rendering and Zoom ----

output$profilePlots <- renderUI({
  req(length(profiles) > 0)

  weights <- sapply(profiles, function(p) p$size_f())
  total <- sum(weights)
  brush_opts <- brushOpts(id = "zoomBrush", direction = "x")

  plot_list <- lapply(seq_along(profiles), function(i) {
    prof <- profiles[[i]]
    plotOutput(
      outputId = paste0("plot_", prof$id),
      height = paste0(300 * weights[i] / total, "px"),
      hover = hoverOpts(id = paste0("hover_", prof$id), delay = 50),
      brush = brush_opts,
      dblclick = "plot_dblclick"
    )
  })

  tagList(plot_list)
})

# setup mapper and plot profiles
observe({
  valid_cids <- intersect(state$contigs, contigs$cid)
  cdf <- contigs[match(valid_cids, contigs$cid), ]
  if (nrow(cdf) == 0) {
    return()
  }

  cdf$start <- cumsum(c(0, head(cdf$length, -1)))
  cdf$end <- cdf$start + cdf$length

  mapper <- list(
    g2l = function(gcoord) {
      idx <- which(gcoord >= cdf$start & gcoord <= cdf$end)
      if (length(idx) == 0) {
        return(NULL)
      }
      list(cid = cdf$cid[idx], coord = gcoord - cdf$start[idx] + 1)
    },
    l2g = function(cid, coord) {
      idx <- which(cdf$cid == cid)
      if (length(idx) == 0) {
        return(NULL)
      }
      cdf$start[idx] + coord - 1
    }
  )

  df_slice <- cdf
  if (!is.null(state$zoom)) {
    df_slice <- df_slice[df_slice$end >= state$zoom[1] & df_slice$start <= state$zoom[2], ]
  }

  for (i in seq_along(profiles)) {
    local({
      prof <- profiles[[i]]
      myid <- prof$id
      output[[paste0("plot_", myid)]] <- renderPlot({
        g <- prof$plot_f(df_slice, mapper)
        g <- g + geom_vline(xintercept = df_slice$start, color = "gray")
        g <- g + labs(y = prof$spec_f()$ylabel)
        g <- g + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(2, 2, 2, 2)
        )
        if (prof$spec_f()$grid) {
          g <- g + theme(panel.grid.major.x = element_line(), panel.grid.major.y = element_line())
        }
        g
      })
    })
  }
})

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

output$profileHoverInfo <- renderText({
  hover_data <- reactiveValuesToList(input)
  for (id in names(profiles)) {
    hname <- paste0("hover_", id)
    if (!is.null(hover_data[[hname]])) {
      h <- hover_data[[hname]]
      gcoord <- as.integer(h$x)
      y <- round(h$y, 2)

      valid_cids <- intersect(state$contigs, contigs$cid)
      cdf <- contigs[match(valid_cids, contigs$cid), ]
      cdf$start <- cumsum(c(0, head(cdf$length, -1)))
      cdf$end <- cdf$start + cdf$length
      cid <- cdf$cid[which(gcoord >= cdf$start & gcoord <= cdf$end)]
      if (length(cid) > 0) {
        coord <- gcoord - cdf$start[cdf$cid == cid] + 1
        return(profiles[[id]]$info_f(cid, coord, y))
      }
    }
  }
  NULL
})
