# ==============================================================================
# Script: 00_install_package.R
# Purpose: Install dependencies and the source package
# ==============================================================================

required_deps <- c("igraph", "MASS", "Matrix")

new_pkgs <- required_deps[!(required_deps %in% installed.packages()[,"Package"])]

if(length(new_pkgs)) {
  message("Installing missing CRAN dependencies: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs)
}

pkg_filename <- "superviseNet_0.1.0.tar.gz"
target_path <- getwd() # Your working directory
pkg_path <- file.path(target_path, pkg_filename)

if(file.exists(pkg_path)) {
  message("----------------------------------------------------------")
  message("Found source package: ", pkg_path)
  message("Installing superviseNet from local source...")
  message("----------------------------------------------------------")

  install.packages(pkg_path, repos = NULL, type = "source")

  if(require("superviseNet")) {
    message("\n>>> SUCCESS: Package 'superviseNet' installed and loaded successfully. <<<")
  } else {
    stop("ERROR: Package installed but failed to load.")
  }

} else {
  stop(paste("ERROR: Could not find", pkg_filename, "in the parent directory."))
}

