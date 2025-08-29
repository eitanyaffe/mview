# Helper for creating a synteny-based profile
# Supports detail/summary modes based on style parameter

# get current binsize for synteny data
get_current_synteny_binsize <- function(xlim, binsize, target_bins = 200, 
    binsizes = c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)) 
{  
  if (binsize != "auto") {
    bs <- as.numeric(binsize)
    if (bs %in% binsizes) {
      return(bs)
    } else {
      warning(sprintf("binsize %d not available, using closest smaller binsize", bs))
      smaller_binsizes <- binsizes[binsizes <= bs]
      if (length(smaller_binsizes) > 0) {
        return(max(smaller_binsizes))
      } else {
        return(min(binsizes))
      }
    }
  }

  if (is.null(xlim) || length(xlim) != 2) {
    return(5000)  # default binsize
  }

  range_bp <- (xlim[2] + 1) - xlim[1]
  if (range_bp <= 0) {
    return(1000)
  }

  # calculate raw bin size to get close to target_bins
  raw_bin_size <- range_bp / target_bins
  
  # choose the closest smaller binsize to ensure we have enough bins
  smaller_binsizes <- binsizes[binsizes <= raw_bin_size]
  if (length(smaller_binsizes) > 0) {
    return(max(smaller_binsizes))
  } else {
    return(min(binsizes))
  }
}

# load synteny data using synteny_f function
load_synteny_data <- function(synteny_f, cxt, binsize, hide_self) {
  if (is.null(synteny_f)) {
    stop("synteny_f function not provided")
  }
  
  cache_key <- paste0("synteny_data_", cxt$assembly, "_", binsize, "_hideself_", hide_self)
  
  cache(cache_key, {
    if (is.function(synteny_f)) {
      synteny_f(cxt, binsize, hide_self)
    } else {
      synteny_f
    }
  })
}

# filter synteny data to current view range
filter_synteny_data <- function(data_list, cxt, xlim) {
  if (is.null(data_list) || is.null(data_list$sequenced_bp) || is.null(data_list$mutations)) {
    return(NULL)
  }
  
  sequenced_bp <- data_list$sequenced_bp
  mutations <- data_list$mutations
  
  # check that both tables have same structure
  if (!identical(sequenced_bp[, 1:3], mutations[, 1:3])) {
    warning("sequenced_bp and mutation tables have different genomic coordinates")
    return(NULL)
  }
  
  # filter to current contigs and range
  valid_rows <- sequenced_bp$contig %in% cxt$contigs
  
  if (!is.null(xlim) && length(xlim) == 2) {
    # filter by genomic range
    valid_rows <- valid_rows & 
                  sequenced_bp$start < xlim[2] & 
                  sequenced_bp$end > xlim[1]
  }
  
  if (!any(valid_rows)) {
    return(NULL)
  }
  
  # return filtered data
  list(
    sequenced_bp = sequenced_bp[valid_rows, ],
    mutations = mutations[valid_rows, ]
  )
}

# default synteny parameters
default_synteny_params <- list(
  style = list(
    group_id = "synteny",
    type = "select",
    choices = c("summary", "detail"),
    default = "summary"
  ),
  binsize = list(
    group_id = "synteny",
    type = "select",
    choices = c("auto", "1000", "2000", "5000", "10000", "20000", "50000", "100000", "200000", "500000"),
    default = "auto"
  ),
  target_bins = list(
    group_id = "synteny",
    type = "integer",
    default = 200
  ),
  min_xcov = list(
    group_id = "synteny",
    type = "double",
    default = 1.0
  ),
  color_style = list(
    group_id = "synteny",
    type = "select",
    choices = c("none", "mutations"),
    default = "none"
  ),
  height = list(
    group_id = "synteny",
    type = "integer",
    default = 400
  ),
  show_hover = list(
    group_id = "synteny",
    type = "boolean",
    default = FALSE
  ),
  hide_self = list(
    group_id = "synteny",
    type = "boolean",
    default = TRUE
  )
)

synteny_profile <- function(id, name,
                           synteny_f = NULL,
                           style = "summary",
                           binsize = "auto",
                           target_bins = 200,
                           min_xcov = 1.0,
                           color_style = "none",
                           height = 400,
                           show_hover = TRUE,
                           hide_self = TRUE,
                           binsizes = c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000),
                           params = default_synteny_params,
                           auto_register = TRUE) {
                           
  plot_f <- function(profile, cxt, gg) {
    # check that intervals dataframe has at least one row
    if (is.null(cxt$intervals) || nrow(cxt$intervals) == 0) {
      warning("intervals dataframe is empty")
      return(list(plot = gg, legends = list()))
    }

    # get current binsize
    current_binsize <- get_current_synteny_binsize(cxt$mapper$xlim, profile$binsize, profile$target_bins, profile$binsizes)
    cat(sprintf("using synteny binsize: %d\n", current_binsize))

    # load synteny data
    data_list <- load_synteny_data(profile$synteny_f, cxt, current_binsize, profile$hide_self)
    if (is.null(data_list)) {
      warning("failed to load synteny data")
      return(list(plot = gg, legends = list()))
    }

    # filter data to current view
    filtered_data <- filter_synteny_data(data_list, cxt, cxt$mapper$xlim)
    if (is.null(filtered_data)) {
      warning("no synteny data in current view")
      return(list(plot = gg, legends = list()))
    }

    # call mode-specific plot function
    if (profile$style == "detail") {
      return(synteny_profile_detail(profile, cxt, filtered_data, gg, current_binsize))
    } else { # style == "summary"
      return(synteny_profile_summary(profile, cxt, filtered_data, gg, current_binsize))
    }
  }

  # create profile
  profile_create(
    id = id, name = name, type = "synteny", height = height,
    params = params, plot_f = plot_f,
    auto_register = auto_register,
    synteny_f = synteny_f,
    style = style,
    binsize = binsize,
    target_bins = target_bins,
    min_xcov = min_xcov,
    color_style = color_style,
    show_hover = show_hover,
    hide_self = hide_self,
    binsizes = binsizes
  )
}
