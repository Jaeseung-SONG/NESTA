# Lightweight NESTA interface checks for the final simulation study bundle.
# Full submitted implementation remains binding at /home/js/NESTA/Analysis/Nesta.R.
reject_nperm <- function(args) {
  if ("nperm" %in% names(args)) stop("Unsupported diffuStats argument `nperm`; use `n.perm`")
  TRUE
}

diffuse_checked <- function(graph, scores, method = "ber_p", n.perm = 300, seed = 9703, ...) {
  dots <- list(...)
  reject_nperm(dots)
  diffuStats::diffuse(graph = graph, scores = scores, method = method, n.perm = n.perm, seed = seed, ...)
}
