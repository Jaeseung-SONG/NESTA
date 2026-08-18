#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if(length(args) >= 1) args[[1]] else file.path("tutorial", "output")
suppressPackageStartupMessages(library(data.table))
script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
required_columns <- c("SYMBOL", "TWAS.Z", "Mean_expression", "Initial.Heat",
                      "Final.Heat", "delta_NESTA", "Analysis_name")

read_and_check <- function(cell_type){
  path <- file.path(output_dir, paste0("cell_type_", cell_type, "_scores.tsv"))
  if(!file.exists(path)) stop("Missing tutorial output: ", path)
  tab <- data.table::fread(path)
  missing <- setdiff(required_columns, names(tab))
  if(length(missing)) stop(path, " is missing: ", paste(missing, collapse = ", "))
  if(any(!is.finite(tab$Final.Heat))) stop(path, " contains non-finite Final.Heat values")
  if(!isTRUE(all.equal(tab$delta_NESTA, tab$Final.Heat - tab$TWAS.Z,
                       tolerance = 1e-10))){
    stop(path, " has an incorrect delta_NESTA definition")
  }
  tab
}

a <- read_and_check("A")
b <- read_and_check("B")
top_genes <- function(x, n = 5) x$SYMBOL[order(-abs(x$Final.Heat), x$SYMBOL)][seq_len(n)]

if(identical(top_genes(a), top_genes(b))){
  stop("Cell type A and B produced identical top-five rankings; expected cell type-specific re-prioritization")
}

observed <- rbind(
  data.table(cell_type = "Cell_type_A", rank = seq_len(5), SYMBOL = top_genes(a)),
  data.table(cell_type = "Cell_type_B", rank = seq_len(5), SYMBOL = top_genes(b))
)
expected <- fread(file.path(script_dir, "expected", "top_genes.tsv"))
if(!identical(observed, expected)){
  stop("Top-gene ranking does not match tutorial/expected/top_genes.tsv")
}

cat("Tutorial validation passed\n")
