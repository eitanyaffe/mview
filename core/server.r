server <- function(input, output, session) {
    source("core/server_state.r", local = TRUE)
    source("core/server_profiles.r", local = TRUE)
    source("core/server_tables.r", local = TRUE)
    source("core/server_events.r", local = TRUE)
    source("core/server_views.r", local = TRUE)
}
