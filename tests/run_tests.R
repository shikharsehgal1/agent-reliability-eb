# Run the test suite from the project root: Rscript tests/run_tests.R
source("R/load.R")
library(testthat)
res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n%d tests, %d expectations, %d failed, %d errors, %d skipped\n",
            nrow(df), sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped)))
if (sum(df$failed) > 0 || any(df$error)) quit(status = 1)
