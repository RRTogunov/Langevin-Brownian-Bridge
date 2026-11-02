#' Path of one checkpoint
#' @param dir   checkpoint directory
#' @param label unique label of the fit (file-name safe)
#' @return file path "<dir>/fit_<label>.rds"
checkpoint_path <- function(dir, label) {
  file.path(dir, sprintf("fit_%s.rds", gsub("[^A-Za-z0-9._-]", "-", label)))
}

#' Write a checkpoint atomically (tmp file + rename; safe for array jobs)
#' @param record list to save
#' @param dir    checkpoint directory (created if missing)
#' @param label  unique label of the fit
#' @return invisible file path
write_checkpoint <- function(record, dir, label) {
  make_path(dir)
  fp  <- checkpoint_path(dir, label)
  tmp <- paste0(fp, ".tmp", Sys.getpid())
  saveRDS(record, tmp)
  file.rename(tmp, fp)
  invisible(fp)
}

#' Does a checkpoint exist?
#' @param dir checkpoint directory; @param label fit label
#' @return logical
checkpoint_exists <- function(dir, label) file.exists(checkpoint_path(dir, label))

#' Read all checkpoints in a directory
#' @param dir     checkpoint directory
#' @param pattern regex on file names
#' @return list of records (empty list if none)
collect_checkpoints <- function(dir, pattern = "^fit_.*\\.rds$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  lapply(files, readRDS)
}
