project_dir <- "/home/js/NESTA/simulation"
report_root <- "/home/js_subdir/Dropbox/NESTA_revision"
binding_plan <- "/home/js_subdir/Dropbox/NESTA_revision/STUDY_PLAN_0709_FALSE_POSITIVE_CONTROL_AUDIT.md"
binding_plan_sha256 <- "463a0eb9bfb2dc4ead568305c9d0df510301a14111bbffccbe8f7c405dfabb93"

project_file <- function(...) file.path(project_dir, ...)

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

file_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
  digest::digest(file = path, algo = "sha256")
}

verify_project_path <- function() {
  resolved <- normalizePath(project_dir, mustWork = TRUE)
  if (!identical(resolved, project_dir)) stop("Project path guard failed: ", resolved)
  info <- file.info(project_dir)
  if (isTRUE(info$isdir) && Sys.readlink(project_dir) != "") {
    stop("Project path is a symbolic link")
  }
  TRUE
}

verify_binding_plan <- function() {
  if (!file.exists(binding_plan) || file.access(binding_plan, 4) != 0) {
    stop("Binding plan missing or unreadable: ", binding_plan)
  }
  observed <- file_sha256(binding_plan)
  if (!identical(observed, binding_plan_sha256)) {
    stop("Binding plan checksum changed: ", observed)
  }
  observed
}

atomic_write_lines <- function(lines, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  safe_dir_create(dirname(path))
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".tmp")
  writeLines(lines, tmp, useBytes = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Atomic rename failed for ", path)
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  safe_dir_create(dirname(path))
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, quote = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Atomic rename failed for ", path)
  invisible(path)
}

atomic_save_rds <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  safe_dir_create(dirname(path))
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(x, tmp)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Atomic rename failed for ", path)
  invisible(path)
}

set_run_status <- function(status, pilot_started = "NO", confirmatory_started = "NO", reason = "") {
  path <- project_file("RUN_STATUS.md")
  if (file.exists(path)) unlink(path)
  atomic_write_lines(c(
    "# Run Status", "",
    paste0("Status: **", status, "**"), "",
    paste0("Pilot execution started: **", pilot_started, "**"), "",
    paste0("Confirmatory execution started: **", confirmatory_started, "**"), "",
    if (nzchar(reason)) paste0("Reason: `", reason, "`.") else "Reason: in progress.", "",
    paste0("Binding plan SHA256: `", binding_plan_sha256, "`")
  ), path)
}

non_overwrite_report_dir <- function(root = report_root) {
  stamp <- format(Sys.time(), "%d%m%y_%H%M")
  base <- file.path(root, paste0("NESTA_simulation_", stamp))
  out <- base
  i <- 1
  while (dir.exists(out) || file.exists(out)) {
    out <- sprintf("%s_%02d", base, i)
    i <- i + 1
  }
  out
}

copy_binding_plan_to_report <- function(report_dir) {
  safe_dir_create(report_dir)
  dst <- file.path(report_dir, basename(binding_plan))
  if (!file.exists(dst)) file.copy(binding_plan, dst, copy.mode = TRUE, copy.date = TRUE)
  if (!identical(file_sha256(dst), binding_plan_sha256)) {
    stop("Copied binding plan checksum mismatch")
  }
  invisible(dst)
}

read_report_dir <- function() {
  env <- Sys.getenv("NESTA_REPORT_DIR", "")
  if (nzchar(env)) return(env)
  for (path in c("/tmp/nesta_0709_fp_audit_report_dir", "/tmp/nesta_0709_strict_fraction_report_dir", "/tmp/nesta_0708_confirmatory_report_dir", "/tmp/nesta_0708_decision_repair_report_dir", "/tmp/nesta_0708_bidirectional_report_dir", "/tmp/nesta_0708_initial_signal_report_dir", "/tmp/nesta_0707_operator_audit_report_dir", "/tmp/nesta_0707_branch_isolation_report_dir", "/tmp/nesta_0707_dense_rescue_report_dir", "/tmp/nesta_0706_staged_repair_report_dir", "/tmp/nesta_0706_size_report_dir", "/tmp/nesta_0705_observed_adaptive_report_dir", "/tmp/nesta_0705_bounded_conductance_report_dir", "/tmp/nesta_0704_opposite_bridge_report_dir", "/tmp/nesta_0703_relay_a2_margin_report_dir", "/tmp/nesta_0703_qc_gate_report_dir", "/tmp/nesta_0703_sweet_spot_report_dir", "/tmp/nesta_restore_1607_report_dir", "/tmp/nesta_relay_restraint_report_dir", "/tmp/nesta_current_report_dir", "/tmp/nesta_relay_a2_compression_report_dir", "/tmp/nesta_relay_gate_report_dir", "/tmp/nesta_diffusion_retention_report_dir", "/tmp/nesta_degree_realigned_report_dir", "/tmp/nesta_network_specific_report_dir", "/tmp/nesta_recovery_first_report_dir", "/tmp/nesta_0628_report_dir", "/tmp/nesta_0626_report_dir", "/tmp/nesta_0626_degree_report_dir", "/tmp/nesta_0625_path_report_dir", "/tmp/nesta_0625_report_dir", "/tmp/nesta_0624_report_dir")) {
    if (file.exists(path)) return(trimws(readLines(path, warn = FALSE)[1]))
  }
  stop("Missing report directory record")
}

safe_sd <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) stop("Non-finite or zero sample SD")
  s
}

sample_or_all <- function(x, n) {
  if (length(x) <= n) x else sample(x, n)
}

ks_stat <- function(x, y) {
  suppressWarnings(as.numeric(stats::ks.test(x, y)$statistic))
}
