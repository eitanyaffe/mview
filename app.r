get_server <- function() {
    server <- function(input, output, session) {
        source("core/server_state.r", local = TRUE)
        source("core/server_profiles.r", local = TRUE)
        source("core/server_tables.r", local = TRUE)
        source("core/server_buttons.r", local = TRUE)
        source("core/server_views.r", local = TRUE)
        source("core/server_parameters.r", local = TRUE)

        # load tabs
        load_tabs()
        source("core/server_tabs.r", local = TRUE)

        states_module_output <- states_server(
            id = "states_module",
            main_state_rv = state,
            session = session
        )

        keyboard_module_output <- keyboard_server(
            input = input,
            output = output,
            session = session,
            main_state_rv = state,
            states_module_output = states_module_output
        )
    }
}

rl_core <- function(project_id) {
    source("app.r", local = FALSE)

    # utils
    source("core/utils.r", local = FALSE)
    source("core/plot_utils.r", local = FALSE)

    # cache
    source("core/cache.r", local = FALSE)

    # initialize cache
    cache_init(project_id = project_id)

    # data
    source("core/data.r", local = FALSE)

    # context
    source("core/context.r", local = FALSE)

    # parameters
    source("core/parameters.r", local = FALSE)

    # tabs system
    source("core/tabs.r", local = FALSE)

    # state table
    source("core/states.R", local = FALSE)

    # keyboard
    source("core/keyboard.r", local = FALSE)

    # profiles
    source("core/profile_manager.r", local = FALSE)
    source("core/profile_viewer.r", local = FALSE)

    # source all profile files
    profile_files <- list.files("profiles", pattern = "\\.r$", recursive = TRUE, full.names = TRUE)
    for (profile_file in profile_files) {
        source(profile_file, local = FALSE)
    }

}

rl <- function(project_id = "c60") {
    # load core modules
    rl_core(project_id = "c60")

    # load specific config
    cfg_file <- paste0("configs/", project_id, "_cfg.r")
    source(cfg_file, local = FALSE)

    # shiny app
    source("core/ui.r", local = TRUE)
    server <- get_server()

    shiny::shinyApp(ui = ui, server = server)
}
