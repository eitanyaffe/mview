get_server <- function() {
    server <- function(input, output, session) {
        source("core/server_state.r", local = TRUE)
        source("core/server_profiles.r", local = TRUE)
        source("core/server_tables.r", local = TRUE)
        source("core/server_buttons.r", local = TRUE)
        source("core/server_views.r", local = TRUE)
        source("core/server_parameters.r", local = TRUE)
        source("core/server_states.r", local = TRUE)
        source("core/server_navigation.r", local = TRUE)

        # load tabs
        load_tabs()
        
        # selector tab (UI function and server code - must be before server_tabs.r)
        source("core/selector_tab.r", local = TRUE)
        
        # segment colors (server-side, after selector_tab)
        source("core/server_segment_colors.r", local = TRUE)
        
        source("core/server_tabs.r", local = TRUE)
        source("core/server_legend.r", local = TRUE)
        
        # load gene regions module
        source("core/gene_regions.r", local = TRUE)

        regions_module_output <- regions_server(
            id = "regions_module",
            main_state_rv = state,
            session = session
        )
        
        gene_regions_module_output <- gene_regions_server(
            id = "gene_regions_module",
            main_state_rv = state,
            session = session
        )

        # load enhanced export dialogs with access to regions module
        source("core/export_dialogs.r", local = TRUE)

        keyboard_module_output <- keyboard_server(
            input = input,
            output = output,
            session = session,
            main_state_rv = state,
            regions_module_output = regions_module_output
        )
    }
}

rl_core <- function(project_id) {
    source("mview.r", local = FALSE)

    # utils
    source("core/utils.r", local = FALSE)
    source("core/plot_utils.r", local = FALSE)

    # cache
    source("core/cache.r", local = FALSE)

    # initialize cache
    cache_init(project_id = project_id)
    
    # reset variant cache on startup
    cache_set("variants.current", NULL)

    # data
    source("core/data.r", local = FALSE)

    # compile C++ context
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Rcpp package is required. Install with: install.packages('Rcpp')")
    }
    library(Rcpp)
    Rcpp::sourceCpp("cpp/R_context.cpp", cacheDir = ".Rcpp_dir", rebuild = FALSE)

    # Ensure .context_env exists globally before loading context.r
    if (!exists(".context_env", envir = .GlobalEnv)) {
      assign(".context_env", new.env(parent = emptyenv()), envir = .GlobalEnv)
      .context_env$contexts <- list()
      .context_env$current_assembly <- NULL
      .context_env$current_context <- NULL
      .context_env$current_segments <- NULL
    }

    # context
    source("core/context.r", local = FALSE)

    # parameters
    source("core/parameters.r", local = FALSE)

    # tabs system
    source("core/tabs.r", local = FALSE)

    # regions table
    source("core/regions.r", local = FALSE)

    # keyboard
    source("core/keyboard.r", local = FALSE)

    # states
    source("core/states.r", local = FALSE)

    # segment colors registry (before config loads)
    source("core/segment_colors.r", local = FALSE)

    # profiles
    source("core/profile_manager.r", local = FALSE)
    source("core/profile_viewer.r", local = FALSE)

    # source all profile files
    profile_files <- list.files("profiles", pattern = "\\.r$", recursive = TRUE, full.names = TRUE)
    for (profile_file in profile_files) {
        source(profile_file, local = FALSE)
    }
}

rl <- function(project_id = "minimal", cdir = "configs") {
    # load core modules
    rl_core(project_id = project_id)

    # load specific config
    cfg_file <- paste0(cdir, "/", project_id, "/", project_id, "_cfg.r")
    if (!file.exists(cfg_file)) {
        stop(paste0("Configuration file not found: ", cfg_file))
    }
    source(cfg_file, local = FALSE)

    # shiny app
    source("core/ui.r", local = TRUE)
    server <- get_server()

    shiny::shinyApp(ui = ui, server = server)
}
