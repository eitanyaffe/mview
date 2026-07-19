########################################################
# RNA profile API
# Loads short_R.cpp (srt_load / srt_query) via Rcpp,
# mirroring the alntools align_profile_api.r pattern.
########################################################

library(Rcpp)

append_flag <- function(cur, add) {
    add <- trimws(add)
    if (nchar(add) == 0) return(cur)
    if (nchar(cur) == 0) return(add)
    paste(cur, add)
}

setup_short <- function(verbose = FALSE) {
    short_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/short/cpp")
    os_name   <- Sys.info()["sysname"]

    existing_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = "")
    existing_libs     <- Sys.getenv("PKG_LIBS",     unset = "")

    if (os_name == "Darwin") {
        if (file.exists("/opt/homebrew/opt/libomp/include/omp.h")) {
            existing_cppflags <- append_flag(existing_cppflags,
                "-I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp")
            existing_libs <- append_flag(existing_libs,
                "-L/opt/homebrew/opt/libomp/lib -lomp")
        }
    }

    Sys.setenv(PKG_CPPFLAGS = existing_cppflags)
    Sys.setenv(PKG_LIBS     = existing_libs)

    bridge_path <- paste0(short_dir, "/short_R.cpp")
    cat("loading short...\n")
    sourceCpp(bridge_path, verbose = verbose, cacheDir = ".Rcpp_dir", rebuild = FALSE)
    cat("short loaded\n")
}

init_short <- function(verbose = FALSE) {
    setup_short(verbose = verbose)
}
