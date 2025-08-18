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

# Setup alntools compilation with optional GPU support
setup_alntools <- function(use_gpu = FALSE, verbose = FALSE) {
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
    
    # GPU support (Metal on macOS)
    if (use_gpu) {
      # Check for both Xcode tools AND actual Metal GPU support
      xcode_ok <- tryCatch(system("xcode-select -p >/dev/null 2>&1", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0, error = function(e) FALSE)
      metal_ok <- tryCatch(length(system("system_profiler SPDisplaysDataType | grep 'Metal'", intern = TRUE)) > 0, error = function(e) FALSE)
      
      metal_available <- xcode_ok && metal_ok
      if (metal_available) {
        existing_cppflags <- append_flag(existing_cppflags, "-DMETAL_SUPPORT")
        existing_libs <- append_flag(existing_libs, "-framework Metal -framework MetalPerformanceShaders -framework Foundation")
        # Force Objective-C++ so Foundation/Metal headers compile
        existing_cxxflags <- append_flag(existing_cxxflags, "-x objective-c++")
        cat("✓ GPU support enabled (Metal)\n")
      } else {
        if (!xcode_ok) {
          cat("⚠ GPU requested but Xcode Command Line Tools not available, using CPU-only\n")
        } else if (!metal_ok) {
          cat("⚠ GPU requested but Metal GPU not detected on this system, using CPU-only\n")
        } else {
          cat("⚠ GPU requested but not available, using CPU-only\n")
        }
        use_gpu <- FALSE
      }
    }
  }
  
  # Set environment variables
  Sys.setenv(PKG_CPPFLAGS = existing_cppflags)
  Sys.setenv(PKG_LIBS = existing_libs)
  if (nchar(existing_cxxflags) > 0) Sys.setenv(PKG_CXXFLAGS = existing_cxxflags)
  
  # Use the unified R bridge file (handles both GPU and CPU)
  bridge_file <- "aln_R.cpp"
  bridge_path <- paste0(alntools_dir, "/", bridge_file)
  
  cat("Loading alntools with", if (use_gpu) "GPU" else "CPU-only", "support...\n")
  
  # Compile and load (avoid rebuilding unless necessary to prevent duplicates)
  sourceCpp(bridge_path, verbose = verbose, cacheDir = ".Rcpp_dir", rebuild = F)
  
  cat("✓ alntools loaded successfully\n")
  return(use_gpu)
}

# Initialize alntools (called by configs)
# Default to CPU-only for compatibility
init_alntools <- function(enable_gpu = FALSE, verbose = FALSE) {
  setup_alntools(use_gpu = enable_gpu, verbose = verbose)
}
