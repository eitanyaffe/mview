server <- function(input, output, session) {
    source("R/server_state.r", local = TRUE)
    source("R/server_profiles.r", local = TRUE)
    source("R/server_tables.r", local = TRUE)
    source("R/server_events.r", local = TRUE)
}
