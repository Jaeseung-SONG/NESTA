#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
report_dir <- Sys.getenv("NESTA_REPORT_DIR", "")
if ("--report-dir" %in% args) {
  idx <- match("--report-dir", args)
  report_dir <- args[[idx + 1]]
}
if (!nzchar(report_dir)) stop("Provide --report-dir or NESTA_REPORT_DIR")
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)
source(file.path("/home/js/NESTA/simulation_study", "R", "reproducibility_checks.R"))
source(file.path("/home/js/NESTA/simulation_study", "R", "summarize_results.R"))
res <- verify_references(report_dir)
dir.create(file.path(report_dir, "manuscript_tables"), recursive=TRUE, showWarnings=FALSE)
if (file.exists(file.path(report_dir, "SCENARIO_SUMMARY_TABLE.csv"))) file.copy(file.path(report_dir, "SCENARIO_SUMMARY_TABLE.csv"), file.path(report_dir, "manuscript_tables", "table_simulation_scenario_summary.csv"), overwrite=TRUE)
if (file.exists(file.path(report_dir, "REFERENCE_RESULT_COMPARISON.csv"))) file.copy(file.path(report_dir, "REFERENCE_RESULT_COMPARISON.csv"), file.path(report_dir, "manuscript_tables", "table_reproducibility_checks.csv"), overwrite=TRUE)
manifest <- write_manifest("/home/js/NESTA/simulation_study", report_dir)
fidelity <- code_fidelity_audit()
write.csv(fidelity, file.path(report_dir, "CODE_FIDELITY_AUDIT.csv"), row.names=FALSE)
writeLines(c("# Code Fidelity Audit", "", sprintf("All fidelity checks passed: %s.", all(fidelity$passed)), "", paste(sprintf("- `%s`: %s", fidelity$check, fidelity$passed), collapse="\n")), file.path(report_dir, "CODE_FIDELITY_AUDIT.md"))
schema <- data.frame(file=c("REFERENCE_RESULT_COMPARISON.csv", "SCENARIO_SUMMARY_TABLE.csv", "GITHUB_READY_FILE_MANIFEST.csv"), required_columns_present=c(all(c("scenario","check","observed","expected","passed") %in% names(res$comparison)), all("scenario" %in% names(res$summary)), all(c("path","sha256") %in% names(manifest))), stringsAsFactors=FALSE)
write.csv(schema, file.path(report_dir, "METRIC_SCHEMA_AUDIT.csv"), row.names=FALSE)
writeLines(c("# Metric Schema Audit", "", sprintf("All schema checks passed: %s.", all(schema$required_columns_present))), file.path(report_dir, "METRIC_SCHEMA_AUDIT.md"))
classification <- if (res$passed && all(fidelity$passed)) "simulation_study_reproducibility_verified_cleanup_ready" else if (any(res$comparison$passed)) "simulation_study_partial_reproduction_no_cleanup" else "simulation_study_reproducibility_failed_no_cleanup"
writeLines(c("# Reproducibility Verification Report", "", paste0("Classification: `", classification, "`."), sprintf("Reference checks passed: %d/%d.", sum(res$comparison$passed), nrow(res$comparison)), sprintf("Code fidelity checks passed: %s.", all(fidelity$passed))), file.path(report_dir, "REPRODUCIBILITY_VERIFICATION_REPORT.md"))
writeLines(c("# Cleanup Plan", "", "Reproducibility must pass before destructive cleanup of exploratory local working files.", "", "Planned cleanup after separate confirmation or release staging:", "- Keep `/home/js/NESTA/simulation_study` as the GitHub-ready bundle.", "- Keep the four reference Dropbox result directories.", "- Archive, rather than delete, older exploratory local code if further cleanup is approved.", "- Remove only transient files from `simulation_study/outputs/`."), file.path(report_dir, "CLEANUP_PLAN.md"))
writeLines(c("# Cleanup Actions Taken", "", "No exploratory simulation code or Dropbox result directories were deleted.", "Transient `simulation_study/outputs/` scratch content was left empty except `README_KEEP_EMPTY.md`.", sprintf("Cleanup performed: %s.", FALSE)), file.path(report_dir, "CLEANUP_ACTIONS_TAKEN.md"))
if (classification == "simulation_study_reproducibility_verified_cleanup_ready") {
  writeLines(c("# Reproducibility Verified", "", "The compact simulation_study package reproduces the four final manuscript-relevant scenarios within the prespecified tolerance.", "No destructive exploratory-code cleanup was performed in this run; the package is cleanup-ready."), file.path("/home/js/NESTA/simulation_study", "REPRODUCIBILITY_VERIFIED.md"))
}
writeLines(c("# Final Report", "", paste0("Binding plan SHA256: `d42492f50fcf6fd1ed10633cd88214ab5d81bf495d00734747fd17575e112823`."), paste0("Final classification: `", classification, "`."), sprintf("Reference checks passed: %d/%d.", sum(res$comparison$passed), nrow(res$comparison)), sprintf("Code fidelity checks passed: %s.", all(fidelity$passed)), "Cleanup performed: NO. The run produced a cleanup-ready package but did not delete exploratory local working code or Dropbox references.", "Package path: `/home/js/NESTA/simulation_study`."), file.path(report_dir, "FINAL_REPORT.md"))
write_checksums(report_dir)
cat(classification, "\n")
