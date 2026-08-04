write_manifest <- function(root, report_dir) {
  files <- list.files(root, recursive=TRUE, full.names=TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  rel <- sub(paste0("^", root, "/"), "", files)
  sha <- vapply(files, function(f) digest::digest(file=f, algo="sha256"), character(1))
  tab <- data.frame(path=rel, size_bytes=file.info(files)$size, sha256=sha, stringsAsFactors=FALSE)
  write.csv(tab, file.path(report_dir, "GITHUB_READY_FILE_MANIFEST.csv"), row.names=FALSE)
  tab
}
write_checksums <- function(report_dir) {
  files <- list.files(report_dir, recursive=TRUE, full.names=TRUE)
  files <- files[file.info(files)$isdir == FALSE & basename(files) != "CHECKSUMS.sha256"]
  lines <- vapply(files, function(f) paste(digest::digest(file=f, algo="sha256"), sub(paste0("^", report_dir, "/"), "", f)), character(1))
  writeLines(lines, file.path(report_dir, "CHECKSUMS.sha256"))
}
