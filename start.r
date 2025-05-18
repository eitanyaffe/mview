# start.r

register_views <- function() {
    view_files <- list.files("views", pattern = "\\.r$", full.names = TRUE)
    for (view_file in view_files) {
        view_name <- tools::file_path_sans_ext(basename(view_file))
        view_register(view_name, view_file)
    }
}

get_server <- function() {
    server <- function(input, output, session) {
        # --- Source Core Server Logic Files ---
        source("core/server_state.r", local = TRUE)
        source("core/server_profiles.r", local = TRUE)
        source("core/server_tables.r", local = TRUE)
        source("core/server_buttons.r", local = TRUE)
        source("core/server_views.r", local = TRUE)
        source("core/server_parameters.r", local = TRUE)

        # Output for the last key pressed
        output$last_key_output <- renderText({
            req(keyboard_module_output$last_key_pressed())
            paste("Pressed:", keyboard_module_output$last_key_pressed())
        })

        states_module_output <- states_server(
            id = "states_module",
            main_state_rv = state,
            session = session
        )

        keyboard_module_output <- keyboard_server(
            input = input,
            session = session,
            states_module_output = states_module_output
        )
    }
}

rl <- function() {
    source("start.r", local = FALSE)

    # utils
    source("core/utils.r", local = FALSE)

    # data
    source("core/data.r", local = FALSE)

    # cache
    source("core/cache.r", local = FALSE)

    # context
    source("core/context.r", local = FALSE)

    # parameters
    source("core/parameters.r", local = FALSE)

    # states
    source("core/states.R", local = FALSE)

    # keyboard module
    source("core/keyboard.r", local = FALSE)

    # profile code
    source("core/profile_manager.r", local = FALSE)
    source("core/profile_viewer.r", local = FALSE)

    # source all profile files
    profile_files <- list.files("profiles", pattern = "\\.r$", recursive = TRUE, full.names = TRUE)
    for (profile_file in profile_files) {
        source(profile_file, local = FALSE)
    }

    # register all views
    register_views()

    # shiny app
    source("core/ui.r", local = TRUE)
    server <- get_server()

    shiny::shinyApp(ui = ui, server = server)
}
