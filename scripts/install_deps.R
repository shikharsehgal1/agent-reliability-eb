# Install the two non-base dependencies into a project-local library (.Rlib),
# leaving the user's global R library untouched. ggplot2 is also required for
# figures; it is installed into .Rlib only if not already available.
dir.create(".Rlib", showWarnings = FALSE)
.libPaths(c(normalizePath(".Rlib"), .libPaths()))
need <- c("lpSolve", "testthat", "ggplot2")
need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need, lib = ".Rlib", repos = "https://cloud.r-project.org")
cat("R", R.version$major, ".", R.version$minor, "\n",
    "lpSolve", as.character(packageVersion("lpSolve")), "\n",
    "testthat", as.character(packageVersion("testthat")), "\n",
    "ggplot2", as.character(packageVersion("ggplot2")), "\n")
