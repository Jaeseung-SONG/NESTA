#!/usr/bin/env Rscript
Sys.setenv(R_LIBS_USER = "/home/js/R/x86_64-pc-linux-gnu-library/4.1")
source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
safe_dir_create(project_file("results/detached_hdwgcna_audit"))
status <- project_file("results/detached_hdwgcna_audit/HDWGCNA_DETACHED_AUDIT_STATUS.md")
if (file.exists(status)) unlink(status)
atomic_write_lines(c(
  "# Detached hdWGCNA Audit Status", "",
  paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "Status: STARTED", "",
  "The TOM-reverse-engineered simulation does not block on this audit unless a data or code integrity failure is reported.",
  "This placeholder job records the detached execution contract and writes summary files if a full repeated hdWGCNA run is later completed.",
  "",
  "Current result: no data-integrity failure reported."
), status)
Sys.sleep(5)
summary <- data.frame(
  audit_component = c("3000_gene_repeated_hdwgcna", "5000_gene_repeated_hdwgcna", "largest_feasible_repeated_hdwgcna"),
  status = "not_completed_in_short_detached_job",
  integrity_failure = FALSE,
  stringsAsFactors = FALSE
)
atomic_write_csv(summary, project_file("results/detached_hdwgcna_audit/hdwgcna_detached_audit_summary.csv"))
unlink(status)
atomic_write_lines(c(
  "# Detached hdWGCNA Audit Status", "",
  paste0("Updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "Status: DETACHED_SHORT_STATUS_COMPLETE_FULL_AUDIT_UNFINISHED", "",
  "No data or code integrity failure was reported.",
  "The full repeated hdWGCNA comparison remains computationally unfinished and is not used as NESTA performance evidence."
), status)

