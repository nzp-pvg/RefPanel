#!/usr/bin/env Rscript

# Resolve the repository root from the calling script rather than a machine-specific path.
release_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script_arg)) stop("Run this file with Rscript so its location can be resolved.")
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", ".."), mustWork = TRUE)
}

repo_root <- release_root()
data_root <- file.path(repo_root, "intermediate_data")
