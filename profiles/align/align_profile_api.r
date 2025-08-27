########################################################
# Shared Alignment Profile API
# Common R interface for alntools across mview configs
########################################################

library(Rcpp)

# Configure build flags helper function
append_flag <- function(cur, add) {
  add <- trimws(add)
  if (nchar(add) == 0) return(cur)
  if (nchar(cur) == 0) return(add)
  paste(cur, add)
}

# Setup alntools compilation
setup_alntools <- function(verbose = FALSE) {
  alntools_dir <- paste0(Sys.getenv("MAKESHIFT_ROOT"), "/tools/alntools/cpp")
  os_name <- Sys.info()["sysname"]
  
  # Reset environment variables
  existing_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = "")
  existing_libs <- Sys.getenv("PKG_LIBS", unset = "")
  existing_cxxflags <- Sys.getenv("PKG_CXXFLAGS", unset = "")
  
  if (os_name == "Darwin") {
    # OpenMP (Homebrew libomp)
    if (file.exists("/opt/homebrew/opt/libomp/include/omp.h")) {
      existing_cppflags <- append_flag(existing_cppflags, "-I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp")
      existing_libs <- append_flag(existing_libs, "-L/opt/homebrew/opt/libomp/lib -lomp")
    }
  }
  
  # Set environment variables
  Sys.setenv(PKG_CPPFLAGS = existing_cppflags)
  Sys.setenv(PKG_LIBS = existing_libs)
  if (nchar(existing_cxxflags) > 0) Sys.setenv(PKG_CXXFLAGS = existing_cxxflags)
  
  # Use the unified R bridge file
  bridge_file <- "aln_R.cpp"
  bridge_path <- paste0(alntools_dir, "/", bridge_file)
  
  cat("Loading alntools...\n")
  
  # Compile and load (avoid rebuilding unless necessary to prevent duplicates)
  sourceCpp(bridge_path, verbose = verbose, cacheDir = ".Rcpp_dir", rebuild = F)
  
  cat("✓ alntools loaded successfully\n")
}

# Initialize alntools (called by configs)
init_alntools <- function(verbose = FALSE) {
  setup_alntools(verbose = verbose)
}
