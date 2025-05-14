# start.r

register_views <- function() {
    view_files <- list.files("views", pattern = "\\.r$", full.names = TRUE)
    for (view_file in view_files) {
        view_name <- tools::file_path_sans_ext(basename(view_file))
        view_register(view_name, view_file)
    }
}

rl <- function() {
    source("start.r", local = FALSE)
    source("core/data.r", local = FALSE)

    # profile manager
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
    source("core/server.r", local = TRUE)

    shiny::shinyApp(ui = ui, server = server)
}
