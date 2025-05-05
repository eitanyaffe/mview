# start.r
profiles <- list()

register_profile <- function(id, profile) {
    profiles[[id]] <<- profile
}

rl <- function() {
    source("start.r", local = FALSE)
    source("R/data.r", local = FALSE)
    source("R/profiles.r", local = FALSE)

    source("R/ui.r", local = TRUE)
    source("R/server.r", local = TRUE)
    shiny::shinyApp(ui = ui, server = server)
}
