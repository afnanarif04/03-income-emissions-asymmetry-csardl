# =====================================================================
# 00_run_all.R — Master replication script
# Set the working directory to the repository root before running.
#
#   setwd("path/to/this/repository")
#   source("code/00_run_all.R")
#
# Runs every script in order and reports elapsed time for each. Output
# files are written to output/, one file per manuscript table (plus a
# small number of supporting "detail" files that back each table's
# summary numbers).
# =====================================================================
rm(list = ls())

code_dir <- "code"
scripts <- c(
  "01_packages_data.R",
  "02_estimator_functions.R",
  "03_descriptive_statistics_table1.R",
  "04_preliminary_tests_table2.R",
  "05_main_estimation_table3.R",
  "06_supplementary_dl_and_pmg.R",
  "07_subsample_sensitivity_table4.R",
  "08_temporal_stability_check.R",
  "09_csa_joint_test_table3.R",
  "10_bias_corrected_table5.R"
)

for (s in scripts) {
  cat(sprintf("--- Running %-40s", s))
  t0 <- proc.time()["elapsed"]
  tryCatch(
    source(file.path(code_dir, s), echo = FALSE),
    error = function(e) cat("\n  *** ERROR:", conditionMessage(e), "\n")
  )
  cat(sprintf(" [%.1f s]\n", proc.time()["elapsed"] - t0))
}

cat("\n=== All scripts complete. Check output/ for one file per manuscript table. ===\n")
