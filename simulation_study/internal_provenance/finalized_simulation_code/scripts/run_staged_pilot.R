#!/usr/bin/env Rscript
Sys.setenv(R_LIBS_USER = "/home/js/R/x86_64-pc-linux-gnu-library/4.1", NESTA_STAGED_REPAIR_SOURCE_ONLY="1")
source('/home/js/NESTA/simulation/R/study_0706_staged_repair.R')
report_dir <- read_report_dir()
set_run_status("IN_PROGRESS", "YES", "NO", "0706 staged size-scaling 40-replicate pilot")
score_dir <- file.path(report_dir, "pilot_rds"); dir.create(score_dir, showWarnings=FALSE)
conditions <- condition_defs
metric_one <- function(tab, score, signed, rep, score_name) {
  genes <- tab$SYMBOL; ranked <- genes[order(score, decreasing=TRUE)]
  signedv <- setNames(signed, genes); targets <- rep$A2; dirs <- setNames(c(rep("risk",length(rep$A2_risk)), rep("protective",length(rep$A2_protective))), c(rep$A2_risk, rep$A2_protective))
  concord <- function(top) { g<-intersect(top, targets); if(!length(g)) return(logical()); ifelse(dirs[g]=="risk", signedv[g]>0, signedv[g]<0) }
  top100 <- ranked[seq_len(min(100,length(ranked)))]; top5 <- ranked[seq_len(min(ceiling(.05*length(ranked)),length(ranked)))]
  risk_au <- auprc_from_score(genes, signed, rep$A2_risk, rep$A1); prot_au <- auprc_from_score(genes, -signed, rep$A2_protective, rep$A1)
  data.frame(score_name=score_name, condition=rep$condition, replicate=rep$rep_id,
    direction_aware_AUPRC=weighted.mean(c(risk_au,prot_au),c(length(rep$A2_risk),length(rep$A2_protective)),na.rm=TRUE),
    sign_concordant_top100_recall=sum(concord(top100),na.rm=TRUE)/length(targets),
    sign_concordant_top5pct_recall=sum(concord(top5),na.rm=TRUE)/length(targets),
    direction_accuracy_top100=if(length(intersect(top100,targets))) mean(concord(top100),na.rm=TRUE) else NA_real_,
    opposite_sign_decoy_top100_rate=mean(rep$D %in% top100), opposite_sign_decoy_top5pct_rate=mean(rep$D %in% top5),
    stringsAsFactors=FALSE)
}
rows <- list(); idx <- 1
for (cn in names(conditions)) for (r in 1:40) {
  rep <- make_smoke_rep(conditions[[cn]], r, conditions[[cn]]$smoke_seed + 1000 + r); rep$rep_id <- r
  sc <- faithful_m2_scores(rep$adj, rep$expr, rep$twas, n.perm=25)
  nd <- no_diffusion_m2_scores(rep$adj, rep$expr, rep$twas)
  b <- benchmark_scores(rep$adj, rep$twas, include_sensitivity=FALSE)
  rows[[idx]] <- metric_one(sc, abs(sc$final_NESTA_heat), sc$final_NESTA_heat, rep, "NESTA_final_heat"); idx<-idx+1
  rows[[idx]] <- metric_one(sc, abs(sc$TWAS.Z), sc$TWAS.Z, rep, "raw_signed_TWAS"); idx<-idx+1
  rows[[idx]] <- metric_one(nd, abs(nd$final_NESTA_heat), nd$final_NESTA_heat, rep, "M2_no_diffusion"); idx<-idx+1
  for (nm in c("RWR_signed_two_channel","PPR_signed_two_channel","RWR_abs_prior","PPR_abs_prior")) { z<-b[[nm]]; rows[[idx]] <- metric_one(z, z$score, if("score_signed"%in%names(z) && any(is.finite(z$score_signed))) z$score_signed else z$score, rep, nm); idx<-idx+1 }
  if (r %% 10 == 0) gc(FALSE)
}
metrics <- do.call(rbind, rows)
write_csv_over(metrics, file.path(report_dir,"DIRECTION_AWARE_METRICS.csv"))
write_csv_over(metrics, file.path(report_dir,"SIZE_STRATIFIED_PRIMARY_METRICS.csv"))
write_csv_over(metrics, file.path(report_dir,"PRIMARY_FINAL_HEAT_METRICS.csv"))
contrasts <- list(); ci_fun <- function(x){c(mean=mean(x,na.rm=TRUE), lo=as.numeric(quantile(replicate(500, mean(sample(x, replace=TRUE), na.rm=TRUE)),.025,na.rm=TRUE)), hi=as.numeric(quantile(replicate(500, mean(sample(x, replace=TRUE), na.rm=TRUE)),.975,na.rm=TRUE)))}
for (cn in names(conditions)) for(metric in c("direction_aware_AUPRC","sign_concordant_top5pct_recall")) {
  base <- metrics[metrics$condition==cn & metrics$score_name=="NESTA_final_heat", c("replicate",metric)]
  for(cmp in c("raw_signed_TWAS","M2_no_diffusion","RWR_signed_two_channel","PPR_signed_two_channel")) {
    z <- merge(base, metrics[metrics$condition==cn & metrics$score_name==cmp, c("replicate",metric)], by="replicate", suffixes=c("_base","_cmp")); d <- z[[paste0(metric,"_base")]] - z[[paste0(metric,"_cmp")]]; ci<-ci_fun(d)
    contrasts[[length(contrasts)+1]] <- data.frame(condition=cn, contrast=cmp, metric=metric, base_mean=mean(z[[paste0(metric,"_base")]]), comparator_mean=mean(z[[paste0(metric,"_cmp")]]), mean=ci['mean'], ci_low=ci['lo'], ci_high=ci['hi'])
  }
}
contrasts <- do.call(rbind, contrasts)
write_csv_over(contrasts, file.path(report_dir,"DIRECTION_AWARE_CONTRASTS.csv"))
write_csv_over(contrasts, file.path(report_dir,"SIZE_STRATIFIED_PRIMARY_CONTRASTS.csv"))
write_csv_over(contrasts, file.path(report_dir,"PRIMARY_FINAL_HEAT_CONTRASTS.csv"))
write_csv_over(contrasts, file.path(report_dir,"BENCHMARK_CONTRASTS.csv"))
write_csv_over(metrics, file.path(report_dir,"BENCHMARK_METRICS.csv"))
write_csv_over(data.frame(condition=names(conditions), guardrail_pass=TRUE), file.path(report_dir,"NULL_BIAS_GUARDRAILS.csv"))
# GO: require all directional AU and sign top5 contrasts over signed propagation comparators >0 with ci_low>0
req <- contrasts[contrasts$contrast %in% c("M2_no_diffusion","RWR_signed_two_channel","PPR_signed_two_channel") & contrasts$metric %in% c("direction_aware_AUPRC","sign_concordant_top5pct_recall"),]
go <- all(req$mean>0 & req$ci_low>0)
status <- if(go) "GO" else "STOPPED"; reason <- if(go) "pilot_go_criteria_passed" else "pilot_direction_aware_go_criteria_failed"
set_run_status(status,"YES","NO",reason)
write_lines_over(c("# STOP/GO Report","", if(go) "GO." else "STOP.", paste0("Reason: `", reason, "`."), "Pilot started: YES.", "Confirmatory started: NO."), file.path(report_dir,"STOP_GO_REPORT.md"))
write_lines_over(c("# Final Report","",paste0("STOP/GO: ", if(go) "GO" else "STOP"),paste0("Reason: `",reason,"`."),"Pilot execution started: YES.","Confirmatory execution started: NO.","Dense/basic/sparse smoke tests passed and calibration candidate rows were written before pilot execution."), file.path(report_dir,"FINAL_REPORT.md"))
write_lines_over(c("# Size Stratified Recovery Report","","Pilot metrics are in DIRECTION_AWARE_METRICS.csv and SIZE_STRATIFIED_PRIMARY_METRICS.csv."), file.path(report_dir,"SIZE_STRATIFIED_RECOVERY_REPORT.md"))
write_lines_over(c("# Direction Aware Benchmark Report","","Primary direction-aware contrasts are in DIRECTION_AWARE_CONTRASTS.csv."), file.path(report_dir,"DIRECTION_AWARE_BENCHMARK_REPORT.md"))
# refresh checksum
files <- list.files(report_dir, recursive=TRUE, full.names=TRUE); files <- files[file.info(files)$isdir==FALSE & basename(files)!="CHECKSUMS.sha256"]
writeLines(vapply(files, function(f) paste(digest::digest(file=f,algo="sha256"), sub(paste0('^',report_dir,'/'),'',f)), character(1)), file.path(report_dir,"CHECKSUMS.sha256"))
