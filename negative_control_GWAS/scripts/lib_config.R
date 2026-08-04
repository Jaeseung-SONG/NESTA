read_simple_config <- function(path) {
  if (!file.exists(path)) stop("Config not found: ", path)
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

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}
