read_simple_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("[[:space:]]+#.*$", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  out <- list()
  for (line in lines) {
    if (!grepl(":", line, fixed = TRUE)) next
    key <- trimws(sub(":.*$", "", line))
    val <- trimws(sub("^[^:]+:", "", line))
    val <- sub("^['\"]", "", sub("['\"]$", "", val))
    if (!nzchar(val)) val <- NA_character_
    out[[key]] <- val
  }
  out
}

cfg_get <- function(cfg, key, required = TRUE, default = NA_character_) {
  val <- cfg[[key]]
  if (is.null(val) || length(val) == 0 || is.na(val) || !nzchar(val)) {
    if (required) stop("Missing required config key: ", key)
    return(default)
  }
  val
}

cfg_num <- function(cfg, key, default = NA_real_) as.numeric(cfg_get(cfg, key, required = FALSE, default = as.character(default)))
cfg_int <- function(cfg, key, default = NA_integer_) as.integer(cfg_get(cfg, key, required = FALSE, default = as.character(default)))
cfg_vec_num <- function(cfg, key) as.numeric(strsplit(cfg_get(cfg, key), ",")[[1]])

script_dir <- function() {
  file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}

ensure_output_tree <- function(root) {
  for (d in c("string_download", "nesta_string", "nesta_coexpression_reference",
              "tables", "reports", "archive_previous_invalid_run")) {
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  }
}

safe_slug <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}
