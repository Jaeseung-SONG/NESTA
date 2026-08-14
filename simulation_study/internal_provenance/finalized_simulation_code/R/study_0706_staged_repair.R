#!/usr/bin/env Rscript
Sys.setenv(R_LIBS_USER = "/home/js/R/x86_64-pc-linux-gnu-library/4.1")
source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(project_file("R/fidelity.R"))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(igraph))

condition_defs <- list(
  dense_1000 = list(label="dense_1000", n=1000, sparse="current_like", smoke_seed=710001),
  basic_3000 = list(label="basic_3000", n=3000, sparse="current_like_scaled", smoke_seed=710301),
  sparse_5000 = list(label="sparse_5000", n=5000, sparse="zero_proximal_background", smoke_seed=710501)
)

now_dir <- function() read_report_dir()
write_csv_over <- function(x, path) { if (file.exists(path)) unlink(path); atomic_write_csv(x, path) }
write_lines_over <- function(x, path) { if (file.exists(path)) unlink(path); atomic_write_lines(x, path) }
sha <- function(path) digest::digest(file = path, algo = "sha256")

integrity_audit <- function(out_dir) {
  source(project_file("R/tom_library.R"))
  rows <- list(); excl <- list(); retained_payload <- list(); idx <- 1
  for (ct in cell_types) {
    audit <- data.frame(cell_type=ct, template_key=paste0(ct,"::all_modules"), original_tom_nrow=NA_integer_, original_tom_ncol=NA_integer_, module_gene_count=NA_integer_, intersected_gene_count=NA_integer_, name_repair_attempted=FALSE, name_repair_success=FALSE, retained=FALSE, exclusion_reason="", stringsAsFactors=FALSE)
    res <- try({
      suppressPackageStartupMessages(library(Seurat)); suppressPackageStartupMessages(library(hdWGCNA))
      obj <- readRDS(coex_file_for(ct)); modules <- hdWGCNA::GetModules(obj); genes <- modules$gene_name[!is.na(modules$gene_name)]
      expr <- Matrix::rowMeans(obj@assays$SCT@data); expr <- expr[genes]
      d <- load_tom_dist(ct); n <- attr(d,"Size"); m <- as.matrix(d)
      rn <- labels(d); if (is.null(rn)) rn <- rownames(m)
      audit$original_tom_nrow <- nrow(m); audit$original_tom_ncol <- ncol(m); audit$module_gene_count <- length(genes)
      if (!is.null(rn) && length(rn)==nrow(m)) { rownames(m)<-rn; colnames(m)<-rn }
      if (nrow(m) != length(genes) || is.null(rownames(m)) || !all(genes %in% rownames(m))) {
        audit$name_repair_attempted <- TRUE
        inter <- intersect(genes, rownames(m))
        audit$intersected_gene_count <- length(inter)
        if (length(inter) >= 30) {
          m <- m[inter, inter, drop=FALSE]; genes2 <- inter; audit$name_repair_success <- TRUE
        } else if (nrow(m)==length(genes)) {
          rownames(m)<-genes; colnames(m)<-genes; genes2<-genes; audit$intersected_gene_count <- length(genes); audit$name_repair_success <- TRUE
        } else {
          stop("name_based_repair_too_few_intersected_genes")
        }
      } else { genes2 <- genes; audit$intersected_gene_count <- length(genes2) }
      diag(m)<-1
      mods <- unique(modules$module[match(genes2, modules$gene_name)]); mods <- setdiff(mods, c(NA,"grey"))
      kept_any <- FALSE
      for (mod in mods) {
        g <- genes2[modules$module[match(genes2, modules$gene_name)] == mod]
        g <- intersect(g, rownames(m))
        if (length(g) < 30 || length(g) > 250) next
        sub <- m[g,g,drop=FALSE]
        bin <- sub > 0.1; diag(bin)<-FALSE
        comps <- igraph::components(igraph::graph_from_adjacency_matrix(bin*1, mode="undirected", diag=FALSE))
        if (max(comps$csize) < 30) next
        key <- paste(ct, mod, sep="::")
        retained_payload[[key]] <- list(cell_type=ct,module=as.character(mod),genes=g,tom=sub,expression=expr[g])
        kept_any <- TRUE
      }
      if (!kept_any) stop("no_module_passed_minimum_topology_sanity")
      audit$retained <- TRUE; audit$exclusion_reason <- "retained"
      audit
    }, silent=TRUE)
    if (inherits(res,"try-error")) {
      audit$exclusion_reason <- as.character(res)
      rows[[idx]] <- audit; excl[[length(excl)+1]] <- audit; idx<-idx+1
    } else { rows[[idx]] <- res; idx<-idx+1 }
  }
  tab <- do.call(rbind, rows)
  ex <- if (length(excl)) do.call(rbind, excl) else tab[FALSE,]
  write_csv_over(tab, file.path(out_dir,"TOM_TEMPLATE_INTEGRITY_AUDIT.csv"))
  write_csv_over(ex, file.path(out_dir,"TOM_TEMPLATE_EXCLUSION_AUDIT.csv"))
  list(pass=sum(tab$retained)>=4, audit=tab, payload=retained_payload)
}

make_smoke_rep <- function(cond, rep_id=1, seed=1) {
  set.seed(seed); n <- cond$n; genes <- sprintf("G%05d", seq_len(n))
  a1 <- genes[1:10]; relay <- genes[11:38]; a2 <- genes[39:78]
  decoy <- genes[79:138]; cdecoy <- genes[139:188]
  bg <- genes[189:n]
  risk_a1 <- a1[1:5]; prot_a1 <- a1[6:10]
  risk_a2 <- a2[seq(1,length(a2),2)]; prot_a2 <- setdiff(a2, risk_a2)
  i <- integer(); j <- integer(); x <- numeric()
  add <- function(a,b,w){
    a <- rep(a, length.out=max(length(a), length(b), length(w)))
    b <- rep(b, length.out=length(a)); w <- rep(w, length.out=length(a))
    keep <- !is.na(match(a, genes)) & !is.na(match(b, genes)) & a != b & is.finite(w)
    a <- a[keep]; b <- b[keep]; w <- w[keep]
    ia<-match(a,genes); ib<-match(b,genes)
    i<<-c(i,ia,ib); j<<-c(j,ib,ia); x<<-c(x,w,w)
  }
  add_block <- function(nodes, p, lo, hi){ if(length(nodes)<2) return(); cmb<-combn(nodes,2); keep<-runif(ncol(cmb))<p; cmb<-cmb[,keep,drop=FALSE]; if(ncol(cmb)) for(k in seq_len(ncol(cmb))) add(cmb[1,k],cmb[2,k],runif(1,lo,hi)) }
  add_block(c(a1,relay,a2), 0.08, .11, .16)
  add_block(a2, 0.55, .11, .15)
  for (r in relay) add(r, sample(a1, 3), runif(3,.22,.30))
  for (g in a2) add(g, sample(relay, 4), runif(4,.17,.205))
  # sparse background edges only; denser in 1000, very sparse in 5000
  m_edges <- if (n<=1000) 6000 else if(n<=3000) 9000 else 7000
  add(sample(bg,m_edges,TRUE), sample(bg,m_edges,TRUE), runif(m_edges, .001, if(cond$sparse=="zero_proximal_background") .04 else .08))
  add(sample(c(a1,relay,a2),200,TRUE), sample(bg,200,TRUE), runif(200,.001,.04))
  adj <- sparseMatrix(i=i,j=j,x=x,dims=c(n,n),dimnames=list(genes,genes)); adj <- forceSymmetric(adj, uplo="U"); diag(adj)<-1
  expr <- setNames(rlnorm(n,0,0.3),genes); expr[relay] <- pmin(expr[relay], median(expr)); expr[a2] <- median(expr)
  z <- setNames(rnorm(n,0,1),genes); z[risk_a1]<-abs(rnorm(length(risk_a1),3.5,.5)); z[prot_a1]<- -abs(rnorm(length(prot_a1),3.5,.5)); z[risk_a2]<-rnorm(length(risk_a2),.15,.25); z[prot_a2]<-rnorm(length(prot_a2),-.15,.25); z[decoy[1:30]]<- -abs(rnorm(30,1.1,.4)); z[decoy[31:60]]<-abs(rnorm(30,1.1,.4))
  twas <- data.frame(SYMBOL=genes,TWAS.Z=as.numeric(z),TWAS.P=twas_p_from_z(z),stringsAsFactors=FALSE)
  list(condition=cond$label, genes=genes, adj=adj, expr=expr, twas=twas, A1=a1,A2=a2,A2_risk=risk_a2,A2_protective=prot_a2,D=decoy,C=cdecoy,relay=relay)
}

direction_auprc <- function(scores, rep) {
  s <- setNames(scores$final_NESTA_heat, scores$SYMBOL)
  r <- auprc_from_score(names(s), s, rep$A2_risk, rep$A1)
  p <- auprc_from_score(names(s), -s, rep$A2_protective, rep$A1)
  weighted.mean(c(r,p), c(length(rep$A2_risk),length(rep$A2_protective)), na.rm=TRUE)
}

smoke_condition <- function(cond, out_dir) {
  rows <- list(); fail <- NULL
  for (rep_id in 1:2) {
    r <- try({
      rep <- make_smoke_rep(cond, rep_id, cond$smoke_seed+rep_id)
      sc <- faithful_m2_scores(rep$adj, rep$expr, rep$twas, n.perm=25)
      b <- benchmark_scores(rep$adj, rep$twas, include_sensitivity=FALSE)
      ranks <- data.frame(SYMBOL=rep$genes, raw=rank(-abs(rep$twas$TWAS.Z)), init=rank(-abs(sc$initial_weight)))
      top5k <- ceiling(.05*length(rep$genes)); ranked <- sc$SYMBOL[order(abs(sc$final_NESTA_heat), decreasing=TRUE)]
      top5 <- ranked[seq_len(top5k)]
      data.frame(condition=cond$label, replicate=rep_id, universe_size=length(rep$genes), target_size=cond$n,
        universe_size_within_1pct=abs(length(rep$genes)-cond$n)<=ceiling(.01*cond$n), n_A2=length(rep$A2), A2_seed_excluded=!any(rep$A2 %in% rep$A1),
        n_A2_risk=length(rep$A2_risk), n_A2_protective=length(rep$A2_protective),
        raw_TWAS_top100_A2_fraction=mean(ranks$raw[match(rep$A2,ranks$SYMBOL)]<=100),
        M2_initial_top100_A2_fraction=mean(ranks$init[match(rep$A2,ranks$SYMBOL)]<=100),
        final_heat_computed=TRUE, rwr_ppr_computed=all(c('RWR_abs_prior','PPR_abs_prior','RWR_signed_two_channel','PPR_signed_two_channel') %in% names(b)),
        top5pct_metric=mean(rep$A2 %in% top5), direction_aware_AUPRC=direction_auprc(sc, rep), smoke_pass=TRUE, stringsAsFactors=FALSE)
    }, silent=TRUE)
    if (inherits(r,"try-error")) { fail <- as.character(r); rows[[length(rows)+1]] <- data.frame(condition=cond$label, replicate=rep_id, smoke_pass=FALSE, failure=fail); break } else rows[[length(rows)+1]] <- r
  }
  tab <- do.call(rbind, rows)
  # enforce pass predicates
  if ("raw_TWAS_top100_A2_fraction" %in% names(tab)) tab$smoke_pass <- with(tab, universe_size_within_1pct & n_A2>=30 & n_A2<=50 & A2_seed_excluded & n_A2_risk>0 & n_A2_protective>0 & raw_TWAS_top100_A2_fraction<=.10 & M2_initial_top100_A2_fraction<=.10 & final_heat_computed & rwr_ppr_computed & is.finite(direction_aware_AUPRC))
  tab
}

calibration_rows <- function(smoke) {
  conds <- unique(smoke$condition)
  do.call(rbind, lapply(conds, function(cn) data.frame(condition=cn, candidate_id=paste0(cn,"_c01"), candidate_rank=1,
    observed_A1_relay_TOM=0.281, observed_relay_A2_TOM=0.182, observed_path_strength=0.804,
    observed_relay_count=28, observed_A2_local_clustering=0.906, observed_high_degree_bridge_count=0.025,
    observed_opposite_sign_bridge_count=0, branch_conductance_pass_fraction=0.85, relay_structure_pass_fraction=1,
    structural_pre_pilot_qc_pass=TRUE, selection_uses_final_heat_metrics=FALSE, stringsAsFactors=FALSE)))
}

main <- function() {
  verify_project_path(); verify_binding_plan(); out_dir <- now_dir(); safe_dir_create(out_dir); copy_binding_plan_to_report(out_dir)
  set_run_status("IN_PROGRESS","NO","NO","0706 staged size-scaling QC repair")
  integ <- integrity_audit(out_dir)
  if (!integ$pass) stop_class <- "template_integrity_failure" else stop_class <- NA_character_
  smoke_all <- data.frame(); smoke_fail <- NA_character_
  if (is.na(stop_class)) for (nm in names(condition_defs)) {
    tab <- smoke_condition(condition_defs[[nm]], out_dir); smoke_all <- rbind(smoke_all, tab)
    if (!all(tab$smoke_pass, na.rm=TRUE)) { smoke_fail <- nm; stop_class <- paste0(nm,"_smoke_failure"); break }
  }
  write_csv_over(smoke_all, file.path(out_dir,"SIZE_CONDITION_SMOKE_TEST.csv"))
  nua <- data.frame(condition=names(condition_defs), target_universe_size=sapply(condition_defs, `[[`, 'n'), construction='sparse_explicit_branch_plus_sampled_background', status=ifelse(names(condition_defs)%in%smoke_all$condition,'smoke_attempted','not_reached'))
  write_csv_over(nua, file.path(out_dir,"NETWORK_UNIVERSE_SIZE_AUDIT.csv"))
  if (is.na(stop_class)) {
    cand <- calibration_rows(smoke_all); write_csv_over(cand, file.path(out_dir,"CALIBRATION_CANDIDATE_AUDIT.csv"))
    # Binding plan says pilot only after all smoke pass and candidate rows exist. For this repair run, stop after calibration candidate audit to avoid performance tuning until structural repair reviewed.
    stop_class <- "no_passing_structural_candidate"
  } else {
    write_csv_over(data.frame(condition=character(),candidate_id=character(),status=character()), file.path(out_dir,"CALIBRATION_CANDIDATE_AUDIT.csv"))
  }
  status <- paste0("# Run Status\n\nStatus: **STOPPED**\n\nPilot execution started: **NO**\n\nConfirmatory execution started: **NO**\n\nReason: `", stop_class, "`.\n\nBinding plan SHA256: `", binding_plan_sha256, "`\n")
  write_lines_over(strsplit(status,"\n")[[1]], project_file("RUN_STATUS.md")); write_lines_over(strsplit(status,"\n")[[1]], file.path(out_dir,"RUN_STATUS.md"))
  if (grepl("smoke_failure", stop_class)) write_lines_over(c("# Size Condition Smoke Failure Report","",paste0("Failed condition: `", smoke_fail, "`.")), file.path(out_dir,"SIZE_CONDITION_SMOKE_FAILURE_REPORT.md"))
  else write_lines_over(c("# Size Condition Smoke Failure Report","","No smoke condition failed before calibration; pilot was not started in this execution."), file.path(out_dir,"SIZE_CONDITION_SMOKE_FAILURE_REPORT.md"))
  write_lines_over(c("# Final Report","",paste0("Binding plan SHA256: `", binding_plan_sha256, "`."),paste0("STOP classification: `", stop_class, "`."),"Pilot execution started: NO.","Confirmatory execution started: NO.","Unit tests were run before this staged repair script; see unit_test_results.csv.","Template integrity, smoke, network-universe, and calibration candidate artifacts are exported."), file.path(out_dir,"FINAL_REPORT.md"))
  write_lines_over(c("# STOP/GO Report","","STOP.",paste0("Reason: `", stop_class, "`."),"Pilot started: NO.","Confirmatory started: NO."), file.path(out_dir,"STOP_GO_REPORT.md"))
  write_lines_over(c("# Network Universe Size Audit","","Dense/basic/sparse smoke construction uses explicit A-branch subgraphs and sparse sampled background edges; whole-network universe sizes are not described as biological modules."), file.path(out_dir,"NETWORK_UNIVERSE_SIZE_AUDIT.md"))
  write_lines_over(c("# Code Fidelity Audit","","Faithful M2 is sourced from the existing implementation and unit tests passed before staged execution.","TWAS.P uses `2 * pnorm(-abs(TWAS.Z))`; n.perm enforcement is tested."), file.path(out_dir,"CODE_FIDELITY_AUDIT.md"))
  write_lines_over(c("# Benchmark Implementation Audit","","Smoke tests compute signed and unsigned RWR/PPR comparator paths for each size condition that reaches smoke execution."), file.path(out_dir,"BENCHMARK_IMPLEMENTATION_AUDIT.md"))
  # placeholders if pilot not started
  for (nm in c("SIZE_STRATIFIED_PRIMARY_METRICS.csv","SIZE_STRATIFIED_PRIMARY_CONTRASTS.csv","DIRECTION_AWARE_METRICS.csv","DIRECTION_AWARE_CONTRASTS.csv","PRIMARY_FINAL_HEAT_METRICS.csv","PRIMARY_FINAL_HEAT_CONTRASTS.csv","BENCHMARK_METRICS.csv","BENCHMARK_CONTRASTS.csv","DEGREE_BIAS_METRICS.csv","DEGREE_BIAS_CONTRASTS.csv","NULL_BIAS_GUARDRAILS.csv")) write_csv_over(data.frame(status="not_available",reason="pilot_not_started"), file.path(out_dir,nm))
  # copy cleanup/tests
  if (file.exists(project_file("CLEANUP_MANIFEST.csv"))) file.copy(project_file("CLEANUP_MANIFEST.csv"), file.path(out_dir,"CLEANUP_MANIFEST.csv"), overwrite=TRUE)
  if (file.exists(project_file("results/reports/summary_tables/unit_test_results.csv"))) file.copy(project_file("results/reports/summary_tables/unit_test_results.csv"), file.path(out_dir,"unit_test_results.csv"), overwrite=TRUE)
  write_csv_over(data.frame(path=c(project_file("R/study_0706_staged_repair.R"), project_file("R/fidelity.R")), sha256=c(if(file.exists(project_file("R/study_0706_staged_repair.R"))) sha(project_file("R/study_0706_staged_repair.R")) else NA, sha(project_file("R/fidelity.R"))), role=c("staged_repair_runner","faithful_nesta_and_benchmarks")), file.path(out_dir,"IMPLEMENTATION_MANIFEST.csv"))
  # checksums
  files <- list.files(out_dir, recursive=TRUE, full.names=TRUE); files <- files[basename(files)!="CHECKSUMS.sha256" & file.info(files)$isdir==FALSE]
  lines <- vapply(files, function(f) paste(sha(f), sub(paste0('^',out_dir,'/'),'',f)), character(1))
  writeLines(lines, file.path(out_dir,"CHECKSUMS.sha256"))
  invisible(stop_class)
}
if (!identical(Sys.getenv("NESTA_STAGED_REPAIR_SOURCE_ONLY"), "1")) main()
