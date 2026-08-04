source(file.path("/home/js/NESTA/simulation_study", "R", "nesta_interface.R"))

read_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
extract_label <- function(final_report) {
  txt <- read_text(final_report)
  m <- regmatches(txt, regexpr("`[^`]+`", sub(".*classification[^`]*", "", txt, ignore.case = TRUE)))
  if (length(m)) gsub("`", "", m) else NA_character_
}
pass_tol <- function(observed, expected, tol = 0.02) is.finite(observed) && abs(observed - expected) <= tol
row1 <- function(x) if (nrow(x)) x[1, , drop = FALSE] else x

scenario_values <- function(ref) {
  base <- basename(ref)
  if (base == "NESTA_simulation_080726_1803") {
    m <- read.csv(file.path(ref, "PRIMARY_FINAL_HEAT_METRICS.csv"))
    s <- aggregate(cbind(top100_recall, top150_recall, top200_recall, risk_top100_recall, protective_top100_recall, opposite_sign_decoy_top100_rate, high_degree_decoy_top100_rate) ~ topology_arm + ranking_mode + rescue_arm, m, mean)
    f <- row1(s[s$topology_arm == "F" & s$ranking_mode == "NESTA_two_tail_balanced", ])
    htxt <- read_text(file.path(ref, "ROBUSTNESS_TOPOLOGY_H_SUMMARY.md"))
    return(list(label = "decision_rule_repair_passed_confirmatory_ready", metrics = c(F_top100=f$top100_recall, F_top150=f$top150_recall, F_top200=f$top200_recall, F_risk=f$risk_top100_recall, F_protective=f$protective_top100_recall, F_opp_decoy=f$opposite_sign_decoy_top100_rate, F_high_decoy=f$high_degree_decoy_top100_rate, H_pass=as.numeric(grepl("pass: TRUE", htxt, fixed=TRUE)))))
  }
  if (base == "NESTA_simulation_080726_1951") {
    p <- read.csv(file.path(ref, "PRIMARY_ENDPOINTS.csv"))
    f <- row1(p[p$topology_arm == "F" & p$ranking_mode == "NESTA_two_tail_balanced", ])
    h <- row1(p[p$topology_arm == "H" & p$ranking_mode == "NESTA_two_tail_balanced", ])
    return(list(label = "confirmatory_success_with_strong_unsigned_ppr_rwr", metrics = c(F_top100=f$top100_recall_mean, F_top150=f$top150_recall_mean, F_top200=f$top200_recall_mean, F_risk=f$risk_top100_recall_mean, F_protective=f$protective_top100_recall_mean, F_opp_decoy=f$opposite_sign_decoy_top100_rate_mean, F_high_decoy=f$high_degree_decoy_top100_rate_mean, H_top100=h$top100_recall_mean, H_top150=h$top150_recall_mean)))
  }
  if (base == "NESTA_simulation_090726_0045") {
    n <- read.csv(file.path(ref, "STRICT_TOP_FRACTION_METRICS.csv"))
    c <- read.csv(file.path(ref, "STRICT_COMPARATOR_METRICS.csv"))
    getn <- function(mode, frac, field) row1(n[n$topology_arm == "F" & n$ranking_mode == mode & abs(n$cutoff_fraction - frac) < 1e-12, ])[[field]]
    getc <- function(score, frac, field) row1(c[c$topology_arm == "F" & c$score_name == score & abs(c$cutoff_fraction - frac) < 1e-12, ])[[field]]
    return(list(label = "strict_fraction_audit_reveals_top100_leniency", metrics = c(PPR_top1=getc("PPR_abs_prior", .01, "total_A2_recall"), PPR_top5=getc("PPR_abs_prior", .05, "total_A2_recall"), RWR_top1=getc("RWR_abs_prior", .01, "total_A2_recall"), RWR_top5=getc("RWR_abs_prior", .05, "total_A2_recall"), NESTA_desc_risk_top50=getn("NESTA_signed_descending", .05, "risk_A2_recall"), NESTA_asc_prot_top50=getn("NESTA_signed_ascending", .05, "protective_A2_recall"))))
  }
  if (base == "NESTA_simulation_090726_0135") {
    t <- read.csv(file.path(ref, "TOPK_TOTAL_FP_METRICS.csv"))
    d <- read.csv(file.path(ref, "TOPK_DIRECTIONAL_FP_METRICS.csv"))
    gett <- function(method, field) row1(t[t$topology_arm == "F" & t$method == method & t$cutoff == 50, ])[[field]]
    getd <- function(method, tail, field) row1(d[d$topology_arm == "F" & d$method == method & d$tail == tail & d$cutoff == 50, ])[[field]]
    return(list(label = "false_positive_audit_completed_ppr_rwr_remain_strong_unsigned", metrics = c(NESTA_top50_recall=gett("NESTA_two_tail_balanced", "recall"), NESTA_top50_precision=gett("NESTA_two_tail_balanced", "precision"), NESTA_top50_FDR=gett("NESTA_two_tail_balanced", "FDR"), NESTA_top50_FPR=gett("NESTA_two_tail_balanced", "FPR"), risk_recall=getd("NESTA_signed_descending", "risk_tail", "recall"), risk_precision=getd("NESTA_signed_descending", "risk_tail", "precision"), risk_FDR=getd("NESTA_signed_descending", "risk_tail", "FDR"), prot_recall=getd("NESTA_signed_ascending", "protective_tail", "recall"), prot_precision=getd("NESTA_signed_ascending", "protective_tail", "precision"), prot_FDR=getd("NESTA_signed_ascending", "protective_tail", "FDR"), PPR_recall=gett("PPR_abs_prior", "recall"), PPR_precision=gett("PPR_abs_prior", "precision"), PPR_FDR=gett("PPR_abs_prior", "FDR"), RWR_recall=gett("RWR_abs_prior", "recall"), RWR_precision=gett("RWR_abs_prior", "precision"), RWR_FDR=gett("RWR_abs_prior", "FDR"))))
  }
  stop("Unknown reference directory: ", ref)
}

verify_references <- function(report_dir, tol = 0.02) {
  refs <- c(
    decision_rule_repair = "/home/js_subdir/Dropbox/NESTA_revision/NESTA_simulation_080726_1803",
    comparator_framed_confirmatory = "/home/js_subdir/Dropbox/NESTA_revision/NESTA_simulation_080726_1951",
    strict_top_fraction_audit = "/home/js_subdir/Dropbox/NESTA_revision/NESTA_simulation_090726_0045",
    false_positive_control_audit = "/home/js_subdir/Dropbox/NESTA_revision/NESTA_simulation_090726_0135"
  )
  expected <- list(
    decision_rule_repair = c(F_top100=.7913, F_top150=.9975, F_top200=1, F_risk=.7650, F_protective=.8175, F_opp_decoy=0, F_high_decoy=0, H_pass=1),
    comparator_framed_confirmatory = c(F_top100=.7795, F_top150=.9988, F_top200=1, F_risk=.7963, F_protective=.7628, F_opp_decoy=0, F_high_decoy=0),
    strict_top_fraction_audit = c(PPR_top1=.0835, PPR_top5=.6766, RWR_top1=.0221, RWR_top5=.3661, NESTA_desc_risk_top50=.7963, NESTA_asc_prot_top50=.7628),
    false_positive_control_audit = c(NESTA_top50_recall=.1255, NESTA_top50_precision=.1004, NESTA_top50_FDR=.8996, NESTA_top50_FPR=.0473, risk_recall=.7963, risk_precision=.3185, risk_FDR=.6815, prot_recall=.7628, prot_precision=.3051, prot_FDR=.6949, PPR_recall=.6766, PPR_precision=.5413, PPR_FDR=.4587, RWR_recall=.3661, RWR_precision=.2929, RWR_FDR=.7071)
  )
  expected_labels <- c(
    decision_rule_repair="decision_rule_repair_passed_confirmatory_ready",
    comparator_framed_confirmatory="confirmatory_success_with_strong_unsigned_ppr_rwr",
    strict_top_fraction_audit="strict_fraction_audit_reveals_top100_leniency",
    false_positive_control_audit="false_positive_audit_completed_ppr_rwr_remain_strong_unsigned"
  )
  rows <- list(); i <- 1
  summary <- list()
  for (nm in names(refs)) {
    vals <- scenario_values(refs[[nm]])
    label_pass <- identical(vals$label, expected_labels[[nm]])
    rows[[i]] <- data.frame(scenario=nm, check="classification", observed=vals$label, expected=expected_labels[[nm]], abs_diff=NA_real_, tolerance=NA_real_, passed=label_pass, stringsAsFactors=FALSE); i <- i+1
    for (metric in names(expected[[nm]])) {
      obs <- unname(vals$metrics[[metric]]); exp <- expected[[nm]][[metric]]
      rows[[i]] <- data.frame(scenario=nm, check=metric, observed=as.character(obs), expected=as.character(exp), abs_diff=abs(obs-exp), tolerance=tol, passed=pass_tol(obs, exp, tol), stringsAsFactors=FALSE); i <- i+1
    }
    summary[[nm]] <- data.frame(scenario=nm, classification=vals$label, metric=names(vals$metrics), value=as.numeric(vals$metrics), stringsAsFactors=FALSE)
  }
  comp <- do.call(rbind, rows)
  sumtab <- do.call(rbind, summary)
  write.csv(comp, file.path(report_dir, "REFERENCE_RESULT_COMPARISON.csv"), row.names=FALSE)
  write.csv(sumtab, file.path(report_dir, "SCENARIO_SUMMARY_TABLE.csv"), row.names=FALSE)
  list(comparison=comp, summary=sumtab, passed=all(comp$passed))
}
