#!/usr/bin/env Rscript
Sys.setenv(R_LIBS_USER = "/home/js/R/x86_64-pc-linux-gnu-library/4.1", NESTA_FP_AUDIT_SOURCE_ONLY = "1")
source(file.path("/home/js/NESTA/simulation", "R", "study_0709_false_positive_audit.R"))
safe_dir_create(project_file("results/reports/summary_tables"))
results <- list()
test <- function(name, expr) {
  msg <- ""; ok <- tryCatch({ force(expr); TRUE }, error = function(e) { msg <<- conditionMessage(e); FALSE })
  results[[length(results) + 1]] <<- data.frame(test = name, passed = ok, message = msg, stringsAsFactors = FALSE)
}
test("path guard", verify_project_path())
test("binding checksum", verify_binding_plan())
test("binding plan is false positive audit", {
  stopifnot(identical(binding_plan, "/home/js_subdir/Dropbox/NESTA_revision/STUDY_PLAN_0709_FALSE_POSITIVE_CONTROL_AUDIT.md"))
  stopifnot(identical(binding_plan_sha256, "463a0eb9bfb2dc4ead568305c9d0df510301a14111bbffccbe8f7c405dfabb93"))
})
test("faithful P conversion and n.perm enforcement", {
  stopifnot(all.equal(twas_p_from_z(c(-2, 0, 2)), 2 * pnorm(-abs(c(-2, 0, 2)))))
  stopifnot(inherits(try(reject_nperm(list(nperm = 1)), silent = TRUE), "try-error"))
  stopifnot(isTRUE(reject_nperm(list(n.perm = 1))))
})
test("exact faithful M2 initial arithmetic", {
  genes <- paste0("g", 1:5); adj <- matrix(.2, 5, 5, dimnames = list(genes, genes)); diag(adj) <- 1
  expr <- setNames(1:5, genes); z <- setNames(c(1, -1, 2, -2, 1), genes)
  tw <- data.frame(SYMBOL = genes, TWAS.Z = z, TWAS.P = twas_p_from_z(z))
  got <- no_diffusion_m2_scores(adj, expr, tw)$initial_weight
  stopifnot(all.equal(as.numeric(got), as.numeric((expr / sd(expr)) * (z / sd(z)))))
})
test("fixed top-k cutoffs are exact", {
  stopifnot(identical(as.integer(fp_cutoffs), c(10L, 50L, 100L)))
  stopifnot(identical(names(fp_cutoffs), c("top10", "top50", "top100")))
})
test("total FP metrics identities", {
  rep <- make_branch_isolation_rep("F", 1, 1)
  selected <- c(rep$A2[1:5], rep$background[1:5])
  row <- total_fp_metrics(rep, selected, 10, "x", "x", "test")
  stopifnot(row$TP == 5, row$FP == 5, row$FN == length(rep$A2) - 5)
  stopifnot(all.equal(row$precision, .5)); stopifnot(all.equal(row$FDR, .5))
  stopifnot(all.equal(row$specificity, 1 - row$FPR))
})
test("directional risk and protective FP definitions", {
  rep <- make_branch_isolation_rep("F", 1, 2)
  score <- setNames(rep(0, length(rep$genes)), rep$genes)
  score[rep$A2_risk[1:3]] <- 10:8
  score[rep$A2_protective[1:2]] <- -10:-9
  row <- directional_fp_metrics(rep, score, "NESTA_two_tail_direction_matched", 10)
  stopifnot(all(c("risk_tail", "protective_tail") %in% row$tail))
  stopifnot(all(c("opposite_direction_false_discoveries", "high_degree_false_discoveries") %in% names(row)))
})
test("decoy FP share is reported", {
  rep <- make_branch_isolation_rep("F", 1, 3)
  selected <- c(rep$A2[1:2], rep$D[1:2], rep$C[1:2], rep$background[1:4])
  row <- decoy_fp_metrics(rep, selected, 10, "x", "x", "test")
  stopifnot(row$total_FP == 8)
  stopifnot(row$opposite_direction_false_discoveries == 2)
  stopifnot(row$high_degree_false_discoveries == 2)
})
test("comparator sets separate unsigned and signed upper-bound", {
  stopifnot(all(c("PPR_abs_prior", "RWR_abs_prior") %in% fp_primary_comparators))
  stopifnot(!any(c("PPR_signed_two_channel", "RWR_signed_two_channel") %in% fp_primary_comparators))
  stopifnot(all(c("PPR_signed_two_channel", "RWR_signed_two_channel") %in% fp_upper_comparators))
})
test("frozen confirmatory seed schedule retained", {
  ss <- make_seed_schedule(); stopifnot(nrow(ss) == 400); stopifnot(all(table(ss$topology_arm) == 200))
})
test("P arm keeps sparse initial field without top100 leakage", {
  rep <- make_branch_isolation_rep("F", 1, 4); rep <- apply_bidirectional_arm(rep, confirmatory_arm, 5)
  ch <- faithful_m2_channels(rep$adj, rep$expr, rep$twas, n.perm = 2, seed = 6)
  no <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
  aud <- initial_field_audit_row(rep, ch, no)
  stopifnot(aud$raw_TWAS_top100_A2_fraction <= 0.05)
  stopifnot(aud$M2_top100_A2_fraction <= 0.10)
})
test("atomic writes and overwrite refusal", {
  p <- tempfile(tmpdir = project_file("results/reports/summary_tables")); atomic_write_lines("x", p)
  stopifnot(file.exists(p)); stopifnot(inherits(try(atomic_write_lines("y", p), silent = TRUE), "try-error")); unlink(p)
})
tab <- do.call(rbind, results)
out <- project_file("results/reports/summary_tables/unit_test_results.csv")
if (file.exists(out)) unlink(out)
atomic_write_csv(tab, out)
if (any(!tab$passed)) { print(tab); stop("Tests failed") }
cat(sprintf("Tests passed: %d; failed: %d\n", sum(tab$passed), sum(!tab$passed)))
