#!/usr/bin/env Rscript
Sys.setenv(R_LIBS_USER = "/home/js/R/x86_64-pc-linux-gnu-library/4.1")
source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))

report_dir <- read_report_dir()
safe_dir_create(report_dir)
safe_dir_create(file.path(report_dir, "summary_tables"))
safe_dir_create(file.path(report_dir, "figure_source_data"))
reports_dir <- project_file("results/reports")
safe_dir_create(reports_dir)
copy_binding_plan_to_report(report_dir)

copy_if_exists <- function(src, dst = file.path(report_dir, basename(src))) {
  if (file.exists(src)) {
    if (file.exists(dst)) unlink(dst)
    file.copy(src, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
}

copy_table <- function(src, name = basename(src)) {
  dst <- file.path(report_dir, name)
  copy_if_exists(src, dst)
  copy_if_exists(src, file.path(report_dir, "summary_tables", name))
}

pilot_dir <- project_file("results/pilot")
null_dir <- project_file("results/null")
topo_dir <- project_file("results/topology_qc")
for (nm in c("PRIMARY_FINAL_HEAT_METRICS.csv", "PRIMARY_FINAL_HEAT_CONTRASTS.csv",
             "DIRECTION_AWARE_METRICS.csv", "DIRECTION_AWARE_CONTRASTS.csv",
             "REPRIORITIZATION_METRICS.csv", "DELTA_INTERPRETATION_METRICS.csv",
             "DIFFUSION_RETENTION_AUDIT.csv",
             "DEGREE_BIAS_METRICS.csv", "DEGREE_BIAS_CONTRASTS.csv",
             "BENCHMARK_METRICS.csv", "BENCHMARK_CONTRASTS.csv",
             "BENCHMARK_DIRECTION_AWARE_METRICS.csv",
             "BENCHMARK_DEGREE_BIAS_METRICS.csv")) {
  copy_table(file.path(pilot_dir, nm), nm)
}
copy_table(file.path(null_dir, "NULL_BIAS_GUARDRAILS.csv"), "NULL_BIAS_GUARDRAILS.csv")
copy_table(project_file("results/calibration/CALIBRATION_CANDIDATE_AUDIT.csv"), "CALIBRATION_CANDIDATE_AUDIT.csv")
copy_if_exists(project_file("results/calibration/CALIBRATION_ADAPTIVE_TRACE.md"), file.path(report_dir, "CALIBRATION_ADAPTIVE_TRACE.md"))
copy_if_exists(project_file("results/calibration/CALIBRATION_ADAPTIVE_TRACE.md"), file.path(reports_dir, "CALIBRATION_ADAPTIVE_TRACE.md"))
for (nm in c("PATH_STRATIFICATION_AUDIT.csv", "DIRECTIONAL_SIGNAL_AUDIT.csv",
             "DEGREE_DISTRIBUTION_AUDIT.csv", "topology_qc_metrics.csv",
             "relevant_irrelevant_matching_qc.csv", "topology_qc_decision.csv",
             "TEMPLATE_OUTLIER_AUDIT.csv", "TARGET_INITIAL_SIGNAL_AUDIT.csv",
             "CONTROL_DISRUPTION_AUDIT.csv", "BRANCH_SPECIFICITY_AUDIT.csv",
             "BRANCH_CONDUCTANCE_AUDIT.csv", "RELAY_STRUCTURE_AUDIT.csv")) {
  copy_table(file.path(topo_dir, nm), nm)
}
copy_table(file.path(topo_dir, "TEMPLATE_OUTLIER_AUDIT.csv"), "TEMPLATE_EXCLUSION_AUDIT.csv")
copy_table(project_file("CLEANUP_MANIFEST.csv"), "CLEANUP_MANIFEST.csv")
copy_table(file.path(reports_dir, "summary_tables", "unit_test_results.csv"), "unit_test_results.csv")
copy_if_exists(project_file("RUN_STATUS.md"), file.path(report_dir, "RUN_STATUS.md"))

ensure_csv <- function(name) {
  p <- file.path(report_dir, name)
  if (!file.exists(p)) atomic_write_csv(data.frame(status = "not_available"), p)
  st <- file.path(report_dir, "summary_tables", name)
  if (!file.exists(st)) copy_if_exists(p, st)
}
for (nm in c("PRIMARY_FINAL_HEAT_METRICS.csv", "PRIMARY_FINAL_HEAT_CONTRASTS.csv",
             "DIRECTION_AWARE_METRICS.csv", "DIRECTION_AWARE_CONTRASTS.csv",
             "REPRIORITIZATION_METRICS.csv", "DELTA_INTERPRETATION_METRICS.csv",
             "DEGREE_BIAS_METRICS.csv", "DEGREE_BIAS_CONTRASTS.csv",
             "BENCHMARK_METRICS.csv", "BENCHMARK_CONTRASTS.csv",
             "BENCHMARK_DIRECTION_AWARE_METRICS.csv",
             "BENCHMARK_DEGREE_BIAS_METRICS.csv",
             "TARGET_INITIAL_SIGNAL_AUDIT.csv", "CONTROL_DISRUPTION_AUDIT.csv",
             "BRANCH_SPECIFICITY_AUDIT.csv", "BRANCH_CONDUCTANCE_AUDIT.csv",
             "RELAY_STRUCTURE_AUDIT.csv", "DIFFUSION_RETENTION_AUDIT.csv",
             "TEMPLATE_OUTLIER_AUDIT.csv", "TEMPLATE_EXCLUSION_AUDIT.csv",
             "NULL_BIAS_GUARDRAILS.csv", "CALIBRATION_CANDIDATE_AUDIT.csv")) ensure_csv(nm)

decision <- NULL
topology_decision <- if (file.exists(file.path(topo_dir, "topology_qc_decision.csv"))) read.csv(file.path(topo_dir, "topology_qc_decision.csv")) else data.frame()
path_qc <- if (file.exists(file.path(topo_dir, "PATH_STRATIFICATION_AUDIT.csv"))) read.csv(file.path(topo_dir, "PATH_STRATIFICATION_AUDIT.csv")) else data.frame()
direction_qc <- if (file.exists(file.path(topo_dir, "DIRECTIONAL_SIGNAL_AUDIT.csv"))) read.csv(file.path(topo_dir, "DIRECTIONAL_SIGNAL_AUDIT.csv")) else data.frame()
degree_qc <- if (file.exists(file.path(topo_dir, "DEGREE_DISTRIBUTION_AUDIT.csv"))) read.csv(file.path(topo_dir, "DEGREE_DISTRIBUTION_AUDIT.csv")) else data.frame()
target_signal_qc <- if (file.exists(file.path(topo_dir, "TARGET_INITIAL_SIGNAL_AUDIT.csv"))) read.csv(file.path(topo_dir, "TARGET_INITIAL_SIGNAL_AUDIT.csv")) else data.frame()
control_disruption_qc <- if (file.exists(file.path(topo_dir, "CONTROL_DISRUPTION_AUDIT.csv"))) read.csv(file.path(topo_dir, "CONTROL_DISRUPTION_AUDIT.csv")) else data.frame()
branch_qc <- if (file.exists(file.path(topo_dir, "BRANCH_SPECIFICITY_AUDIT.csv"))) read.csv(file.path(topo_dir, "BRANCH_SPECIFICITY_AUDIT.csv")) else data.frame()
conductance_qc <- if (file.exists(file.path(topo_dir, "BRANCH_CONDUCTANCE_AUDIT.csv"))) read.csv(file.path(topo_dir, "BRANCH_CONDUCTANCE_AUDIT.csv")) else data.frame()
relay_qc <- if (file.exists(file.path(topo_dir, "RELAY_STRUCTURE_AUDIT.csv"))) read.csv(file.path(topo_dir, "RELAY_STRUCTURE_AUDIT.csv")) else data.frame()
primary <- read.csv(file.path(report_dir, "PRIMARY_FINAL_HEAT_METRICS.csv"))
contrasts <- read.csv(file.path(report_dir, "PRIMARY_FINAL_HEAT_CONTRASTS.csv"))
direction <- read.csv(file.path(report_dir, "DIRECTION_AWARE_METRICS.csv"))
dir_contrasts <- read.csv(file.path(report_dir, "DIRECTION_AWARE_CONTRASTS.csv"))
reprior <- read.csv(file.path(report_dir, "REPRIORITIZATION_METRICS.csv"))
delta <- read.csv(file.path(report_dir, "DELTA_INTERPRETATION_METRICS.csv"))
diffusion_retention <- read.csv(file.path(report_dir, "DIFFUSION_RETENTION_AUDIT.csv"))
degree <- read.csv(file.path(report_dir, "DEGREE_BIAS_METRICS.csv"))
bench <- read.csv(file.path(report_dir, "BENCHMARK_METRICS.csv"))
bench_contrasts <- read.csv(file.path(report_dir, "BENCHMARK_CONTRASTS.csv"))
bench_direction <- read.csv(file.path(report_dir, "BENCHMARK_DIRECTION_AWARE_METRICS.csv"))
bench_degree <- read.csv(file.path(report_dir, "BENCHMARK_DEGREE_BIAS_METRICS.csv"))
null_guard <- read.csv(file.path(report_dir, "NULL_BIAS_GUARDRAILS.csv"))

has_cols <- function(x, cols) all(cols %in% names(x))
safe_aggregate <- function(formula, data, fun = mean, ...) {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  tryCatch(
    aggregate(formula, data, fun, ...),
    error = function(e) data.frame()
  )
}
fmt_bool <- function(x) if (!length(x) || is.na(x)) "UNKNOWN" else if (isTRUE(as.logical(x))) "PASS" else "FAIL"
fmt_num <- function(x) if (!length(x) || all(is.na(x))) "NA" else sprintf("%.4f", mean(x, na.rm = TRUE))
pick <- function(x, setting = "rare_target_detection", score = "NESTA_final_heat", target = "A2_intermediate_degree_capped") {
  if (!has_cols(x, c("difficulty_setting", "score_name", "network_label"))) return(data.frame())
  y <- x[x$difficulty_setting == setting & x$score_name == score & x$network_label == "relevant", , drop = FALSE]
  if ("target_set" %in% names(y)) y <- y[y$target_set == target, , drop = FALSE]
  y
}

run_status_lines <- if (file.exists(project_file("RUN_STATUS.md"))) readLines(project_file("RUN_STATUS.md"), warn = FALSE) else character()
decision_status <- if (any(grepl("Status: \\*\\*GO\\*\\*", run_status_lines))) "GO" else if (any(grepl("Status: \\*\\*IN_PROGRESS\\*\\*", run_status_lines))) "IN_PROGRESS" else "STOPPED"
decision_reason <- if (any(grepl("^Reason:", run_status_lines))) sub("^Reason: `?([^`.]+).*", "\\1", run_status_lines[grep("^Reason:", run_status_lines)[1]]) else "study_not_completed"
pilot_started <- if (any(grepl("Pilot execution started: \\*\\*YES\\*\\*", run_status_lines))) "YES" else "NO"
confirm_started <- if (any(grepl("Confirmatory execution started: \\*\\*YES\\*\\*", run_status_lines))) "YES" else "NO"
stop_go <- if (identical(decision_status, "GO")) "GO" else "STOP"
null_guard_status <- if (has_cols(null_guard, "guardrail_pass")) fmt_bool(all(null_guard$guardrail_pass)) else "UNKNOWN"

rare <- pick(primary)
rare_dir <- pick(direction)
rare_degree <- pick(degree)
rare_reprior <- if (has_cols(reprior, "difficulty_setting")) reprior[reprior$difficulty_setting == "rare_target_detection", , drop = FALSE] else data.frame()
rare_delta <- if (has_cols(delta, "difficulty_setting")) delta[delta$difficulty_setting == "rare_target_detection", , drop = FALSE] else data.frame()
rare_contrasts <- if (has_cols(contrasts, c("difficulty_setting", "target_set", "metric"))) {
  contrasts[contrasts$difficulty_setting == "rare_target_detection" &
              contrasts$target_set == "A2_intermediate_degree_capped" &
              contrasts$metric == "top100_recall", , drop = FALSE]
} else data.frame()

template_audit_source <- file.path(topo_dir, "TEMPLATE_OUTLIER_AUDIT.csv")
template_audit <- if (file.exists(template_audit_source)) {
  read.csv(template_audit_source)
} else if (nrow(path_qc) && "template_key" %in% names(path_qc)) {
  aggregate(cbind(direct_1hop_fraction, two_hop_fraction, three_plus_hop_fraction,
                  median_A1_A2_path, path_fallback_used) ~ difficulty_setting + template_key,
            path_qc, mean, na.rm = TRUE)
} else data.frame(template_key = character())
if (file.exists(file.path(report_dir, "TEMPLATE_OUTLIER_AUDIT.csv"))) unlink(file.path(report_dir, "TEMPLATE_OUTLIER_AUDIT.csv"))
atomic_write_csv(template_audit, file.path(report_dir, "TEMPLATE_OUTLIER_AUDIT.csv"))
copy_if_exists(file.path(report_dir, "TEMPLATE_OUTLIER_AUDIT.csv"),
               file.path(report_dir, "summary_tables", "TEMPLATE_OUTLIER_AUDIT.csv"))
copy_if_exists(file.path(report_dir, "TEMPLATE_OUTLIER_AUDIT.csv"),
               file.path(report_dir, "TEMPLATE_EXCLUSION_AUDIT.csv"))
copy_if_exists(file.path(report_dir, "TEMPLATE_EXCLUSION_AUDIT.csv"),
               file.path(report_dir, "summary_tables", "TEMPLATE_EXCLUSION_AUDIT.csv"))

write_report <- function(name, lines) {
  p <- file.path(report_dir, name)
  if (file.exists(p)) unlink(p)
  atomic_write_lines(lines, p)
  copy_if_exists(p, file.path(reports_dir, name))
}

write_csv_report <- function(x, path) {
  if (file.exists(path)) unlink(path)
  atomic_write_csv(x, path)
}

contrast_lines <- if (nrow(rare_contrasts)) {
  apply(rare_contrasts[, c("contrast", "base_mean", "comparator_mean", "mean_fold_improvement", "mean", "ci_low", "ci_high")], 1,
        function(x) {
          base_mean <- as.numeric(x[["base_mean"]])
          comparator_mean <- as.numeric(x[["comparator_mean"]])
          fold_txt <- if (!is.finite(comparator_mean) || comparator_mean < 0.05) {
            "fold not interpreted; sparse comparator mean near zero"
          } else {
            sprintf("mean ratio %.3f", base_mean / comparator_mean)
          }
          sprintf("- `%s`: final heat %.4f vs comparator %.4f, %s, paired difference %.4f [%.4f, %.4f]",
                  x[["contrast"]], base_mean, comparator_mean, fold_txt,
                  as.numeric(x[["mean"]]), as.numeric(x[["ci_low"]]), as.numeric(x[["ci_high"]]))
        })
} else "- Rare-setting contrasts unavailable."

write_report("CODE_FIDELITY_AUDIT.md", c(
  "# Code Fidelity Audit", "",
  paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),
  "Submitted implementation reference: `/home/js/NESTA/Analysis/Nesta.R`.",
  "Faithful M2 uses SCT row-mean expression inputs, uncentered sample-SD scaling, strict `TWAS.P < cutoff`, signed positive and absolute-negative diffusion, retained zero-weight edges, submitted self-loop behavior, and `diffuStats` `n.perm`.",
  "Synthetic TWAS P values are computed as `2 * pnorm(-abs(TWAS.Z))`.",
  "Primary ranking is `abs(final_NESTA_heat)`. `delta_NESTA` is auxiliary interpretation only.",
  paste0("Implementation audit status: ", if (file.exists(file.path(report_dir, "unit_test_results.csv"))) "tests_completed" else "tests_not_found", ".")
))

write_report("NETWORK_GENERATOR_AUDIT.md", c(
  "# Network Generator Audit", "",
  paste0("Topology QC: ", fmt_bool(topology_decision$topology_qc_pass), "."),
  paste0("Relevant hard-global pass fraction: ", fmt_num(topology_decision$relevant_hard_global_pass_fraction), "."),
  paste0("Rare relevant hard-global pass fraction: ", fmt_num(topology_decision$rare_relevant_hard_global_pass_fraction), "."),
  paste0("Rare relevant module-local pass fraction: ", fmt_num(topology_decision$rare_relevant_module_local_pass_fraction), "."),
  paste0("Path-stratification QC: ", fmt_bool(topology_decision$path_stratification_qc_pass), "."),
  paste0("Degree-distribution QC: ", fmt_bool(topology_decision$degree_distribution_qc_pass), "."),
  paste0("Directional QC: ", fmt_bool(topology_decision$directional_qc_pass), "."),
  paste0("Target initial-signal QC: ", fmt_bool(topology_decision$target_initial_signal_qc_pass), "."),
  paste0("Control-disruption QC: ", fmt_bool(topology_decision$control_disruption_qc_pass), "."),
  paste0("Branch-specificity QC: ", fmt_bool(topology_decision$branch_specificity_qc_pass), "."),
  paste0("Branch-conductance QC: ", fmt_bool(topology_decision$branch_conductance_qc_pass), "."),
  paste0("Relay-structure QC: ", fmt_bool(topology_decision$relay_structure_qc_pass), "."),
  "Binding settings generated for this round: rare target detection and extreme sparse target detection. Conventional module recovery is excluded from GO/STOP criteria."
))

write_report("PATH_STRATIFICATION_AUDIT.md", c(
  "# Path Stratification Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$path_stratification_qc_pass), "."),
  paste0("Mean direct 1-hop fraction: ", fmt_num(path_qc$direct_1hop_fraction), "."),
  paste0("Mean 2-hop fraction: ", fmt_num(path_qc$two_hop_fraction), "."),
  paste0("Mean 3-plus-hop fraction: ", fmt_num(path_qc$three_plus_hop_fraction), "."),
  paste0("Mean median A1-A2 path: ", fmt_num(path_qc$median_A1_A2_path), ".")
))

write_report("DIRECTIONAL_SIGNAL_AUDIT.md", c(
  "# Directional Signal Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$directional_qc_pass), "."),
  paste0("Mean A1 risk seeds: ", fmt_num(direction_qc$n_A1_risk), "."),
  paste0("Mean A1 protective seeds: ", fmt_num(direction_qc$n_A1_protective), "."),
  paste0("Mean opposite-sign decoy genes: ", fmt_num(direction_qc$n_D_opposite_sign_decoy), ".")
))

write_report("DEGREE_DISTRIBUTION_AUDIT.md", c(
  "# Degree Distribution Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$degree_distribution_qc_pass), "."),
  paste0("Mean A2 target count: ", fmt_num(degree_qc$n_A2_primary), "."),
  paste0("Mean fraction extreme-high-degree A2: ", fmt_num(degree_qc$fraction_A2_extreme_high_degree), "."),
  paste0("Mean median A2 degree percentile: ", fmt_num(degree_qc$median_A2_degree_percentile), "."),
  paste0("Degree hard-gate pass fraction: ", fmt_num(degree_qc$degree_hard_gate_pass), "."),
  paste0("Rare target-count pilot-eligible fraction: ",
         fmt_num(degree_qc$target_count_pilot_eligible[degree_qc$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Rare low-degree warning fraction: ",
         fmt_num(degree_qc$median_A2_degree_percentile_warning_low[degree_qc$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Rare low target-count warning fraction: ",
         fmt_num(degree_qc$target_count_warning_low[degree_qc$difficulty_setting == "rare_target_detection"]), "."),
  "Low median A2 degree percentile, low moderate-connector fraction, and A2 target counts between 20 and 29 are warning-only conditions in this 0705 observed-metric adaptive calibration round. Excessive extreme-high-degree A2 is the pre-pilot hard-stop degree risk."
))

write_report("BENCHMARK_IMPLEMENTATION_AUDIT.md", c(
  "# Benchmark Implementation Audit", "",
  "`RWR_abs_prior` and `PPR_abs_prior` are unsigned magnitude-only topology comparators.",
  "`RWR_signed_two_channel` and `PPR_signed_two_channel` propagate positive and negative TWAS-derived priors separately.",
  "All benchmarks use weighted networks and matched target labels."
))

rare_target_signal <- if (has_cols(target_signal_qc, "difficulty_setting")) {
  target_signal_qc[target_signal_qc$difficulty_setting == "rare_target_detection", , drop = FALSE]
} else data.frame()
write_report("TARGET_INITIAL_SIGNAL_AUDIT.md", c(
  "# Target Initial-Signal Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$target_initial_signal_qc_pass), "."),
  paste0("Rare fraction A2 in raw TWAS top 100: ", fmt_num(rare_target_signal$fraction_A2_in_raw_TWAS_top100), "."),
  paste0("Rare fraction A2 in raw TWAS top 200: ", fmt_num(rare_target_signal$fraction_A2_in_raw_TWAS_top200), "."),
  paste0("Rare fraction A2 in M2 initial top 100: ", fmt_num(rare_target_signal$fraction_A2_in_M2_initial_top100), "."),
  paste0("Rare fraction A2 in M2 initial top 200: ", fmt_num(rare_target_signal$fraction_A2_in_M2_initial_top200), "."),
  paste0("Rare median raw TWAS rank A2: ", fmt_num(rare_target_signal$median_raw_TWAS_rank_A2), "."),
  paste0("Rare median M2 initial rank A2: ", fmt_num(rare_target_signal$median_M2_initial_rank_A2), "."),
  paste0("Rare fallback target-selection fraction: ",
         fmt_num(rare_target_signal$target_initial_signal_fallback_fraction), ".")
))

rare_control <- if (has_cols(control_disruption_qc, "difficulty_setting")) {
  control_disruption_qc[control_disruption_qc$difficulty_setting == "rare_target_detection", , drop = FALSE]
} else data.frame()
write_report("CONTROL_DISRUPTION_AUDIT.md", c(
  "# Control Disruption Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$control_disruption_qc_pass), "."),
  paste0("Rare relevant median A1-A2 path: ", fmt_num(rare_control$relevant_median_A1_A2_path), "."),
  paste0("Rare I2 median A1-A2 path: ", fmt_num(rare_control$I2_median_A1_A2_path), "."),
  paste0("Rare I3 median A1-A2 path: ", fmt_num(rare_control$I3_median_A1_A2_path), "."),
  paste0("Rare I2 within-A TOM reduction: ", fmt_num(rare_control$within_A_TOM_reduction_fraction_I2), "."),
  paste0("Rare I3 within-A TOM reduction: ", fmt_num(rare_control$within_A_TOM_reduction_fraction_I3), "."),
  paste0("Rare signed branch preservation relevant/I2/I3: ",
         fmt_num(rare_control$signed_branch_preservation_relevant), " / ",
         fmt_num(rare_control$signed_branch_preservation_I2), " / ",
         fmt_num(rare_control$signed_branch_preservation_I3), ".")
))

rare_branch <- if (has_cols(branch_qc, "difficulty_setting")) {
  branch_qc[branch_qc$difficulty_setting == "rare_target_detection", , drop = FALSE]
} else data.frame()
write_report("BRANCH_SPECIFICITY_AUDIT.md", c(
  "# Branch Specificity Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$branch_specificity_qc_pass), "."),
  paste0("Rare within-risk branch TOM: ", fmt_num(rare_branch$within_risk_branch_TOM), "."),
  paste0("Rare within-protective branch TOM: ", fmt_num(rare_branch$within_protective_branch_TOM), "."),
  paste0("Rare between risk/protective TOM: ", fmt_num(rare_branch$between_risk_protective_TOM), "."),
  paste0("Rare A2-background TOM: ", fmt_num(rare_branch$A2_background_TOM), "."),
  paste0("Rare A2-high-degree-decoy TOM: ", fmt_num(rare_branch$A2_high_degree_decoy_TOM), "."),
  paste0("Rare within-to-between TOM ratio: ", fmt_num(rare_branch$within_to_between_TOM_ratio), "."),
  paste0("Rare within-to-background TOM ratio: ", fmt_num(rare_branch$within_to_background_TOM_ratio), "."),
  paste0("Rare fraction A2 with single dominant branch: ",
         fmt_num(rare_branch$fraction_A2_with_single_dominant_branch), ".")
))

rare_conductance <- if (has_cols(conductance_qc, "difficulty_setting")) {
  conductance_qc[conductance_qc$difficulty_setting == "rare_target_detection", , drop = FALSE]
} else data.frame()
write_report("BRANCH_CONDUCTANCE_AUDIT.md", c(
  "# Branch Conductance Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$branch_conductance_qc_pass), "."),
  paste0("Rare A-branch internal edge weight: ", fmt_num(rare_conductance$A_branch_internal_edge_weight), "."),
  paste0("Rare A-branch background cut weight: ", fmt_num(rare_conductance$A_branch_background_cut_weight), "."),
  paste0("Rare A-branch background cut fraction: ", fmt_num(rare_conductance$A_branch_background_cut_fraction), "."),
  paste0("Rare A-branch internal edge fraction: ", fmt_num(rare_conductance$A_branch_internal_edge_fraction), "."),
  paste0("Rare A1-relay edge weight: ", fmt_num(rare_conductance$A1_relay_edge_weight), "."),
  paste0("Rare relay-A2 edge weight: ", fmt_num(rare_conductance$relay_A2_edge_weight), "."),
  paste0("Rare A1-A2 effective path strength: ", fmt_num(rare_conductance$A1_A2_effective_path_strength), "."),
  paste0("Rare A2 local clustering: ", fmt_num(rare_conductance$A2_local_clustering), "."),
  paste0("Rare A2 branch volume fraction: ", fmt_num(rare_conductance$A2_branch_volume_fraction), "."),
  paste0("Rare seed-neighborhood background fraction: ", fmt_num(rare_conductance$seed_neighborhood_background_fraction), "."),
  paste0("Rare high-degree bridge count: ", fmt_num(rare_conductance$high_degree_bridge_count), "."),
  paste0("Rare opposite-sign bridge count: ", fmt_num(rare_conductance$opposite_sign_bridge_count), ".")
))

rare_relay <- if (has_cols(relay_qc, "difficulty_setting")) {
  relay_qc[relay_qc$difficulty_setting == "rare_target_detection", , drop = FALSE]
} else data.frame()
write_report("RELAY_STRUCTURE_AUDIT.md", c(
  "# Relay Structure Audit", "",
  paste0("QC status: ", fmt_bool(topology_decision$relay_structure_qc_pass), "."),
  paste0("Rare relay gene count: ", fmt_num(rare_relay$n_relay_genes), "."),
  paste0("Rare relay raw TWAS top-100 fraction: ", fmt_num(rare_relay$relay_raw_TWAS_top100_fraction), "."),
  paste0("Rare relay M2 initial top-100 fraction: ", fmt_num(rare_relay$relay_M2_initial_top100_fraction), "."),
  paste0("Rare relay M2 initial top-100 warning fraction: ", fmt_num(rare_relay$relay_M2_initial_top100_warning), "."),
  paste0("Rare low relay-count warning fraction: ", fmt_num(rare_relay$relay_count_warning_low), "."),
  paste0("Rare A1-relay mean TOM: ", fmt_num(rare_relay$A1_relay_mean_TOM), "."),
  paste0("Rare relay-A2 mean TOM: ", fmt_num(rare_relay$relay_A2_mean_TOM), "."),
  paste0("Rare relay-background mean TOM: ", fmt_num(rare_relay$relay_background_mean_TOM), "."),
  paste0("Rare relay-high-degree-decoy TOM: ", fmt_num(rare_relay$relay_high_degree_decoy_TOM), "."),
  paste0("Rare fraction A2 with relay path: ", fmt_num(rare_relay$fraction_A2_with_relay_path), "."),
  paste0("Rare median A1-relay-A2 path strength: ", fmt_num(rare_relay$median_A1_relay_A2_path_strength), ".")
))

rare_diffusion <- if (has_cols(diffusion_retention, "difficulty_setting")) {
  diffusion_retention[diffusion_retention$difficulty_setting == "rare_target_detection", , drop = FALSE]
} else data.frame()
write_report("DIFFUSION_RETENTION_AUDIT.md", c(
  "# Diffusion Retention Audit", "",
  if (identical(pilot_started, "NO")) {
    "Diffusion-retention metrics are unavailable because pre-pilot QC stopped the study before pilot execution."
  } else {
    "Diffusion-retention metrics are populated from pilot Final NESTA heat scores."
  },
  paste0("Rare fraction seed heat retained in A branch: ",
         fmt_num(rare_diffusion$fraction_seed_heat_retained_in_A_branch), "."),
  paste0("Rare fraction seed heat reaching A2: ",
         fmt_num(rare_diffusion$fraction_seed_heat_reaching_A2), "."),
  paste0("Rare fraction seed heat leaking to background: ",
         fmt_num(rare_diffusion$fraction_seed_heat_leaking_to_background), "."),
  paste0("Rare fraction seed heat leaking to high-degree decoys: ",
         fmt_num(rare_diffusion$fraction_seed_heat_leaking_to_high_degree_decoys), "."),
  paste0("Rare fraction seed heat leaking to opposite-sign decoys: ",
         fmt_num(rare_diffusion$fraction_seed_heat_leaking_to_opposite_sign_decoys), "."),
  paste0("Rare A2 final-heat rank compression top200-to-top100: ",
         fmt_num(rare_diffusion$A2_final_heat_rank_compression_from_top200_to_top100), "."),
  paste0("Rare relay gene top-100 rate: ", fmt_num(rare_diffusion$relay_gene_top100_rate), "."),
  paste0("Rare relay gene top-50 rate: ", fmt_num(rare_diffusion$relay_gene_top50_rate), "."),
  paste0("Rare relay-to-A2 top-100 ratio: ", fmt_num(rare_diffusion$relay_to_A2_top100_ratio), "."),
  paste0("Rare relay heat mass: ", fmt_num(rare_diffusion$relay_heat_mass), "."),
  paste0("Rare A2 heat mass: ", fmt_num(rare_diffusion$A2_heat_mass), "."),
  paste0("Rare relay-to-A2 heat ratio: ", fmt_num(rare_diffusion$relay_to_A2_heat_ratio), "."),
  paste0("Rare fraction A2 reached via relay: ", fmt_num(rare_diffusion$fraction_A2_reached_via_relay), ".")
))

write_report("NESTA_FINAL_HEAT_BENCHMARK_REPORT.md", c(
  "# NESTA Final Heat Benchmark Report", "",
  "Primary score: `abs(final_NESTA_heat)`.",
  "Rare target detection setting is the preferred manuscript-facing setting.",
  paste0("Rare top-100 recall: ", fmt_num(rare$top100_recall), "."),
  paste0("Rare top-150 recall: ", fmt_num(rare$top150_recall), "."),
  paste0("Rare top-200 recall: ", fmt_num(rare$top200_recall), "."),
  paste0("Rare top-100 fold enrichment over random: ", fmt_num(rare$top100_fold_enrichment_over_random), "."),
  paste0("Rare raw AUPRC: ", fmt_num(rare$raw_AUPRC), "."),
  "Rare top-100 contrasts:",
  contrast_lines,
  "",
  "Raw AUPRC is reported but intentionally not the only display metric because sparse target prevalence makes raw AUPRC numerically low."
))

write_report("NESTA_DELTA_INTERPRETATION_REPORT.md", c(
  "# NESTA Delta Interpretation Report", "",
  "`delta_NESTA = final_NESTA_heat - TWAS_Z` was retained as an auxiliary biological interpretation metric.",
  paste0("Rare median `delta_NESTA` among A2: ", fmt_num(rare_delta$median_delta_NESTA_A2), "."),
  paste0("Rare median absolute `delta_NESTA` among A2: ", fmt_num(rare_delta$median_abs_delta_NESTA_A2), "."),
  paste0("Rare direction-concordant `delta_NESTA` fraction: ", fmt_num(rare_delta$fraction_A2_delta_NESTA_direction_concordant), ".")
))

write_report("NESTA_REPRIORITIZATION_REPORT.md", c(
  "# NESTA Reprioritization Report", "",
  "A2_TWAS_blank definition: `abs(TWAS.Z) < 1.0`.",
  paste0("Rare mean promoted A2 count from outside raw TWAS top-200 to final heat top-100: ",
         fmt_num(rare_reprior$number_of_A2_genes_promoted_from_outside_raw_TWAS_top200_to_final_heat_top100), "."),
  paste0("Rare fraction of final heat top-100 A2 with TWAS.P > 0.10: ",
         fmt_num(rare_reprior$fraction_of_final_heat_top100_A2_with_TWAS.P_gt_0.10), "."),
  paste0("Rare median A2 rank improvement final heat vs raw TWAS: ",
         fmt_num(rare_reprior$median_rank_improvement_A2_final_heat_vs_raw_TWAS), ".")
))

write_report("NESTA_DEGREE_BIAS_REPORT.md", c(
  "# NESTA Degree Bias Report", "",
  paste0("Rare score-degree Spearman: ", fmt_num(rare_degree$score_degree_spearman), "."),
  paste0("Rare score-strength Spearman: ", fmt_num(rare_degree$score_strength_spearman), "."),
  paste0("Rare high-degree decoy top-100 rate: ", fmt_num(rare_degree$C_high_degree_decoy_top100_rate), "."),
  paste0("Rare top-100 degree percentile median: ", fmt_num(rare_degree$top100_degree_percentile_median), "."),
  paste0("Pre-pilot rare median A2 degree percentile: ",
         fmt_num(degree_qc$median_A2_degree_percentile[degree_qc$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Pre-pilot rare target-count warning fraction: ",
         fmt_num(degree_qc$target_count_warning_low[degree_qc$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Pre-pilot rare low-degree warning fraction: ",
         fmt_num(degree_qc$median_A2_degree_percentile_warning_low[degree_qc$difficulty_setting == "rare_target_detection"]), "."),
  "Pre-pilot low-degree warnings are diagnostic only. Pilot hard degree-bias guardrails are absolute score-degree Spearman <= 0.15 and high-degree decoy top-100 rate <= 0.12."
))

signed_rwr <- if (has_cols(direction, c("score_name", "network_label"))) {
  bench_direction[bench_direction$score_name == "RWR_signed_two_channel" & bench_direction$network_label == "relevant", , drop = FALSE]
} else data.frame()
signed_ppr <- if (has_cols(direction, c("score_name", "network_label"))) {
  bench_direction[bench_direction$score_name == "PPR_signed_two_channel" & bench_direction$network_label == "relevant", , drop = FALSE]
} else data.frame()
unsigned_rwr <- if (has_cols(bench, c("score_name", "network_label", "target_set"))) {
  bench[bench$score_name == "RWR_abs_prior" & bench$network_label == "relevant" &
          bench$target_set == "A2_intermediate_degree_capped", , drop = FALSE]
} else data.frame()
unsigned_ppr <- if (has_cols(bench, c("score_name", "network_label", "target_set"))) {
  bench[bench$score_name == "PPR_abs_prior" & bench$network_label == "relevant" &
          bench$target_set == "A2_intermediate_degree_capped", , drop = FALSE]
} else data.frame()
signed_rwr_abs <- if (has_cols(bench, c("score_name", "network_label", "target_set"))) {
  bench[bench$score_name == "RWR_signed_two_channel" & bench$network_label == "relevant" &
          bench$target_set == "A2_intermediate_degree_capped", , drop = FALSE]
} else data.frame()
signed_ppr_abs <- if (has_cols(bench, c("score_name", "network_label", "target_set"))) {
  bench[bench$score_name == "PPR_signed_two_channel" & bench$network_label == "relevant" &
          bench$target_set == "A2_intermediate_degree_capped", , drop = FALSE]
} else data.frame()
write_report("RWR_PPR_DIRECTIONAL_BENCHMARK_REPORT.md", c(
  "# RWR/PPR Directional Benchmark Report", "",
  paste0("NESTA rare sign-concordant top-100 recall: ", fmt_num(rare_dir$sign_concordant_top100_recall), "."),
  paste0("Signed RWR sign-concordant top-100 recall: ", fmt_num(signed_rwr$sign_concordant_top100_recall[signed_rwr$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Signed PPR sign-concordant top-100 recall: ", fmt_num(signed_ppr$sign_concordant_top100_recall[signed_ppr$difficulty_setting == "rare_target_detection"]), "."),
  paste0("NESTA rare opposite-sign decoy top-100 rate: ", fmt_num(rare_dir$opposite_sign_decoy_top100_rate), "."),
  paste0("Signed RWR opposite-sign decoy top-100 rate: ", fmt_num(signed_rwr$opposite_sign_decoy_top100_rate[signed_rwr$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Signed PPR opposite-sign decoy top-100 rate: ", fmt_num(signed_ppr$opposite_sign_decoy_top100_rate[signed_ppr$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Unsigned RWR absolute top-100 recall: ", fmt_num(unsigned_rwr$top100_recall[unsigned_rwr$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Unsigned PPR absolute top-100 recall: ", fmt_num(unsigned_ppr$top100_recall[unsigned_ppr$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Signed RWR absolute top-100 recall: ", fmt_num(signed_rwr_abs$top100_recall[signed_rwr_abs$difficulty_setting == "rare_target_detection"]), "."),
  paste0("Signed PPR absolute top-100 recall: ", fmt_num(signed_ppr_abs$top100_recall[signed_ppr_abs$difficulty_setting == "rare_target_detection"]), "."),
  "Unsigned PPR/RWR are magnitude-only topology comparators and are reported separately from direction-aware benchmarks."
))

gradient <- if (has_cols(primary, c("difficulty_setting", "score_name", "network_label", "target_set"))) {
  aggregate(cbind(top100_recall, top100_fold_enrichment_over_random, raw_AUPRC) ~ difficulty_setting,
            primary[primary$score_name == "NESTA_final_heat" & primary$network_label == "relevant" &
                      primary$target_set == "A2_intermediate_degree_capped", ],
            mean, na.rm = TRUE)
} else data.frame()
write_report("DIFFICULTY_GRADIENT_REPORT.md", c(
  "# Difficulty Gradient Report", "",
  if (nrow(gradient)) {
    apply(gradient, 1, function(x) sprintf("- `%s`: top-100 recall %.4f, fold enrichment %.4f, raw AUPRC %.4f",
                                           x[["difficulty_setting"]], as.numeric(x[["top100_recall"]]),
                                           as.numeric(x[["top100_fold_enrichment_over_random"]]),
                                           as.numeric(x[["raw_AUPRC"]])))
  } else {
    "Difficulty-gradient performance metrics are not available because QC stopped before pilot execution."
  }
))

contrast_diff <- function(nm) {
  if (!nrow(rare_contrasts) || !("contrast" %in% names(rare_contrasts))) return(NA_real_)
  x <- rare_contrasts$mean[rare_contrasts$contrast == nm]
  if (length(x)) x[1] else NA_real_
}
control_loss <- any(c(contrast_diff("M2_no_diffusion"),
                      contrast_diff("I2_module_disrupted"),
                      contrast_diff("I3_expression_matched_randomized"),
                      contrast_diff("median_weight_permutation")) < 0, na.rm = TRUE)
failure_mode <- if (has_cols(rare_diffusion, c("fraction_seed_heat_leaking_to_background")) &&
           mean(rare_diffusion$fraction_seed_heat_leaking_to_background, na.rm = TRUE) > 0.80) {
  "Seed heat still leaks excessively from the A branch into background."
} else if (has_cols(rare_diffusion, c("fraction_seed_heat_reaching_A2")) &&
           mean(rare_diffusion$fraction_seed_heat_reaching_A2, na.rm = TRUE) < 0.025) {
  "Seed heat does not reach A2 targets at the planned level."
} else if (has_cols(rare_diffusion, c("relay_to_A2_top100_ratio", "relay_to_A2_heat_ratio")) &&
           (mean(rare_diffusion$relay_to_A2_top100_ratio, na.rm = TRUE) > 0.50 ||
              mean(rare_diffusion$relay_to_A2_heat_ratio, na.rm = TRUE) > 0.95)) {
  "Relay genes still absorb too much heat or top-100 rank mass instead of passing it to A2."
} else if (has_cols(rare_diffusion, c("relay_gene_top100_rate")) &&
           mean(rare_diffusion$relay_gene_top100_rate, na.rm = TRUE) > 0.140) {
  "Relay diagnostic occupancy exceeded the 0703 diagnostic ceiling, indicating residual relay/intermediate-gene top-rank occupancy."
} else if (has_cols(rare_diffusion, c("A2_final_heat_rank_compression_from_top200_to_top100")) &&
           mean(rare_diffusion$A2_final_heat_rank_compression_from_top200_to_top100, na.rm = TRUE) < 0.145) {
  "A2 signal was not compressed from the top-200 into the top-100 at the 0703 QC-gate realignment recovery threshold."
} else if (has_cols(rare, c("top100_recall")) &&
           mean(rare$top100_recall, na.rm = TRUE) < 0.120) {
  "Final Heat did not meet the minimum A2 top-100 recovery threshold."
} else if (!control_loss) {
  "No primary control failure was detected in rare target detection."
} else if (has_cols(rare_target_signal, c("fraction_A2_in_M2_initial_top100")) &&
           mean(rare_target_signal$fraction_A2_in_M2_initial_top100, na.rm = TRUE) > 0.10) {
  "A2 targets remain too enriched by M2 initial weight."
} else if (has_cols(rare_branch, c("within_to_background_TOM_ratio")) &&
           mean(rare_branch$within_to_background_TOM_ratio, na.rm = TRUE) < 1.75) {
  "Branch specificity was insufficient relative to background TOM."
} else if (has_cols(rare_diffusion, c("A2_final_heat_rank_compression_from_top200_to_top100")) &&
           mean(rare_diffusion$A2_final_heat_rank_compression_from_top200_to_top100, na.rm = TRUE) < 0.15) {
  "A2 heat was retained in the top-200 but not compressed into the top-100."
} else if (has_cols(rare_control, c("control_disruption_qc_pass")) &&
           !all(rare_control$control_disruption_qc_pass, na.rm = TRUE)) {
  "I2/I3 controls did not sufficiently disrupt signed A1-A2 topology."
} else if (has_cols(rare_degree, c("score_degree_spearman")) &&
           abs(mean(rare_degree$score_degree_spearman, na.rm = TRUE)) > 0.20) {
  "Target placement or scoring appears degree/path driven rather than relevant-topology driven."
} else {
  "Diffusion may be over-smoothing signal away from A2 targets, or controls still retain enough topology to recover them."
}
control_sensitivity_lines <- c(
  "# Control Sensitivity Diagnosis", "",
  "Network-level ablation against I2/I3 is secondary but important in this network-specific recovery round.",
  paste0("Network-control sensitivity detected: ", if (control_loss) "YES" else "NO", "."),
  paste0("M2_no_diffusion top-100 paired difference: ", fmt_num(contrast_diff("M2_no_diffusion")), "."),
  paste0("I2 top-100 paired difference: ", fmt_num(contrast_diff("I2_module_disrupted")), "."),
  paste0("I3 top-100 paired difference: ", fmt_num(contrast_diff("I3_expression_matched_randomized")), "."),
  paste0("Median weight-permutation top-100 paired difference: ", fmt_num(contrast_diff("median_weight_permutation")), "."),
  paste0("Likely failure mode: ", failure_mode)
)
write_report("CONTROL_SENSITIVITY_DIAGNOSIS.md", control_sensitivity_lines)
write_report("CONTROL_FAILURE_DIAGNOSIS.md", control_sensitivity_lines)

if (identical(pilot_started, "NO")) {
  interpretation <- "Final NESTA heat recovery, degree-bias behavior, opposite-sign decoy suppression, and `delta_NESTA` interpretability could not be evaluated because pre-pilot QC stopped the study before pilot execution."
} else {
  interpretation <- paste0(
    "Final NESTA heat ",
    if (nrow(rare_contrasts) && any(rare_contrasts$contrast == "raw_TWAS_abs") &&
        rare_contrasts$mean[rare_contrasts$contrast == "raw_TWAS_abs"][1] > 0) "improved" else "did not clearly improve",
    " sparse TWAS-weak non-seed target recovery over raw TWAS in the rare target detection setting. ",
    "Degree bias was ",
    if (!nrow(rare_degree) || !"score_degree_spearman" %in% names(rare_degree) ||
        !is.finite(mean(rare_degree$score_degree_spearman, na.rm = TRUE))) {
      "not evaluated"
    } else if (abs(mean(rare_degree$score_degree_spearman, na.rm = TRUE)) <= 0.15) {
      "controlled"
    } else if (abs(mean(rare_degree$score_degree_spearman, na.rm = TRUE)) <= 0.20) {
      "partially controlled"
    } else {
      "not controlled"
    },
    " by the score-degree guardrail. `delta_NESTA` is useful as an interpretability metric for reprioritized genes, not as the primary endpoint."
  )
}

write_report("MANUSCRIPT_READY_SUMMARY.md", c(
  "# Manuscript-Ready Summary", "",
  "Because final NESTA heat is the direct output of network diffusion, it was used as the primary prioritization score. `delta_NESTA` was retained as an auxiliary interpretability metric to identify genes whose network-supported priority increased relative to their direct TWAS evidence.",
  "",
  interpretation,
  "",
  "Sparse target prevalence makes raw AUPRC numerically low; fold enrichment over random ranking and fold improvement over raw TWAS are the manuscript-facing display metrics."
))

write_report("STOP_GO_REPORT.md", c(
  "# STOP/GO Report", "",
  paste0("Decision: **", stop_go, "**."),
  paste0("Status: `", decision_status, "`."),
  paste0("Reason: `", decision_reason, "`."),
  paste0("Pilot started: **", pilot_started, "**."),
  paste0("Confirmatory execution started: **", confirm_started, "**."),
  paste0("Null-bias guardrails: ", null_guard_status, ".")
))

write_report("FINAL_REPORT.md", c(
  "# Final Report", "",
  paste0("Output directory: `", report_dir, "`."),
  paste0("Binding plan checksum: `", binding_plan_sha256, "`."),
  paste0("Topology QC: ", fmt_bool(topology_decision$topology_qc_pass), "."),
  paste0("Path-stratification QC: ", fmt_bool(topology_decision$path_stratification_qc_pass), "."),
  paste0("Degree-distribution QC: ", fmt_bool(topology_decision$degree_distribution_qc_pass), "."),
  paste0("Directional QC: ", fmt_bool(topology_decision$directional_qc_pass), "."),
  paste0("Target initial-signal QC: ", fmt_bool(topology_decision$target_initial_signal_qc_pass), "."),
  paste0("Control-disruption QC: ", fmt_bool(topology_decision$control_disruption_qc_pass), "."),
  paste0("Branch-specificity QC: ", fmt_bool(topology_decision$branch_specificity_qc_pass), "."),
  paste0("Branch-conductance QC: ", fmt_bool(topology_decision$branch_conductance_qc_pass), "."),
  paste0("Relay-structure QC: ", fmt_bool(topology_decision$relay_structure_qc_pass), "."),
  paste0("Rare relevant module-local pass fraction: ", fmt_num(topology_decision$rare_relevant_module_local_pass_fraction), "."),
  paste0("Null-bias guardrails: ", null_guard_status, "."),
  paste0("Pilot started: **", pilot_started, "**."),
  paste0("Confirmatory execution started: **", confirm_started, "**."),
  paste0("STOP/GO decision: **", stop_go, "**."),
  "",
  "Rare target detection, final NESTA heat:",
  paste0("Top-100 recall: ", fmt_num(rare$top100_recall), "."),
  paste0("Top-150 recall: ", fmt_num(rare$top150_recall), "."),
  paste0("Top-200 recall: ", fmt_num(rare$top200_recall), "."),
  paste0("Top100/top200 compression ratio: ",
         fmt_num(rare_diffusion$A2_final_heat_rank_compression_from_top200_to_top100), "."),
  paste0("Seed heat retained in A branch: ",
         fmt_num(rare_diffusion$fraction_seed_heat_retained_in_A_branch), "."),
  paste0("Seed heat reaching A2: ",
         fmt_num(rare_diffusion$fraction_seed_heat_reaching_A2), "."),
  paste0("Seed heat leaking to background: ",
         fmt_num(rare_diffusion$fraction_seed_heat_leaking_to_background), "."),
  paste0("Relay gene top-100 rate: ", fmt_num(rare_diffusion$relay_gene_top100_rate), "."),
  paste0("Relay-to-A2 top-100 ratio: ", fmt_num(rare_diffusion$relay_to_A2_top100_ratio), "."),
  paste0("Relay-to-A2 heat ratio: ", fmt_num(rare_diffusion$relay_to_A2_heat_ratio), "."),
  paste0("Fraction A2 reached via relay: ", fmt_num(rare_diffusion$fraction_A2_reached_via_relay), "."),
  paste0("Top-100 fold enrichment over random: ", fmt_num(rare$top100_fold_enrichment_over_random), "."),
  paste0("Raw AUPRC: ", fmt_num(rare$raw_AUPRC), "."),
  paste0("Sign-concordant top-100 recall: ", fmt_num(rare_dir$sign_concordant_top100_recall), "."),
  paste0("Opposite-sign decoy top-100 rate: ", fmt_num(rare_dir$opposite_sign_decoy_top100_rate), "."),
  paste0("Score-degree Spearman: ", fmt_num(rare_degree$score_degree_spearman), "."),
  "",
  "Main top-100 contrasts:",
  contrast_lines,
  "",
  paste0("Control failure diagnosis: ", failure_mode),
  "",
  interpretation
))

fig_dir <- file.path(report_dir, "figure_source_data")
panel_a <- if (nrow(path_qc)) aggregate(n_A2_primary ~ difficulty_setting, path_qc, mean, na.rm = TRUE) else data.frame()
write_csv_report(panel_a, file.path(fig_dir, "panel_A_simulation_design_counts.csv"))
panel_b <- if (has_cols(primary, c("difficulty_setting", "score_name", "network_label", "target_set", "prevalence"))) safe_aggregate(prevalence ~ difficulty_setting, primary[primary$score_name == "NESTA_final_heat" & primary$network_label == "relevant" & primary$target_set == "A2_intermediate_degree_capped", ], mean, na.rm = TRUE) else data.frame()
write_csv_report(panel_b, file.path(fig_dir, "panel_B_difficulty_prevalence.csv"))
panel_c <- if (has_cols(primary, c("difficulty_setting", "score_name", "network_label", "target_set", "top100_fold_enrichment_over_random"))) safe_aggregate(top100_fold_enrichment_over_random ~ difficulty_setting, primary[primary$score_name == "NESTA_final_heat" & primary$network_label == "relevant" & primary$target_set == "A2_intermediate_degree_capped", ], mean, na.rm = TRUE) else data.frame()
write_csv_report(panel_c, file.path(fig_dir, "panel_C_final_heat_fold_enrichment.csv"))
direction_panel <- if (has_cols(bench_direction, names(direction))) rbind(direction, bench_direction[, names(direction), drop = FALSE]) else direction
panel_d <- if (has_cols(direction_panel, c("difficulty_setting", "score_name", "network_label", "sign_concordant_top100_recall"))) safe_aggregate(sign_concordant_top100_recall ~ difficulty_setting + score_name, direction_panel[direction_panel$network_label == "relevant", ], mean, na.rm = TRUE) else data.frame()
write_csv_report(panel_d, file.path(fig_dir, "panel_D_sign_concordant_recovery.csv"))
degree_panel <- if (has_cols(bench_degree, names(degree))) rbind(degree, bench_degree[, names(degree), drop = FALSE]) else degree
panel_e <- if (has_cols(degree_panel, c("difficulty_setting", "score_name", "network_label", "C_high_degree_decoy_top100_rate", "score_degree_spearman"))) safe_aggregate(cbind(C_high_degree_decoy_top100_rate, score_degree_spearman) ~ difficulty_setting + score_name, degree_panel[degree_panel$network_label == "relevant", ], mean, na.rm = TRUE) else data.frame()
write_csv_report(panel_e, file.path(fig_dir, "panel_E_decoy_and_degree_bias.csv"))
panel_f <- if (has_cols(reprior, c("difficulty_setting", "number_of_A2_genes_promoted_from_outside_raw_TWAS_top200_to_final_heat_top100", "fraction_of_final_heat_top100_A2_with_TWAS.P_gt_0.10"))) safe_aggregate(cbind(number_of_A2_genes_promoted_from_outside_raw_TWAS_top200_to_final_heat_top100, fraction_of_final_heat_top100_A2_with_TWAS.P_gt_0.10) ~ difficulty_setting, reprior, mean, na.rm = TRUE) else data.frame()
write_csv_report(panel_f, file.path(fig_dir, "panel_F_reprioritized_blank_targets.csv"))

cat("Reports written to ", report_dir, "\n", sep = "")
