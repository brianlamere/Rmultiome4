#helper functions used by project, but not involved in the scientific workflow
update_provenance <- function(seurat_obj, milestone, params = NULL) {
  setwd(Rmultiome_files)
  git_hash <- tryCatch(system("git rev-parse HEAD", intern = TRUE),
                       error=function(e) NA)
  git_branch <- tryCatch(system("git rev-parse --abbrev-ref HEAD",
                                intern = TRUE), error=function(e) NA)
  git_tag <- tryCatch(system("git describe --tags --always",
                             intern = TRUE), error=function(e) NA)
  git_dirty <- tryCatch(length(system("git status --porcelain",
                                      intern = TRUE)) > 0, error=function(e) NA)
  info <- list(
    step = milestone,
    timestamp = Sys.time(),
    git_commit = git_hash,
    git_branch = git_branch,
    git_tag = git_tag,
    git_dirty = git_dirty,
    params = params
  )
  if (is.null(seurat_obj@misc$provenance)) {
    seurat_obj@misc$provenance <- list()
  }
  seurat_obj@misc$provenance[[length(seurat_obj@misc$provenance) + 1]] <- info
  seurat_obj
}

#this will print the provenance log.
print_provenance <- function(seurat_obj) {
  prov <- seurat_obj@misc$provenance
  if (is.null(prov) || length(prov) == 0) {
    cat("No provenance information found.\n")
    return(invisible(NULL))
  }
  for (i in seq_along(prov)) {
    cat(sprintf("Step %d: %s\n", i, prov[[i]]$step))
    cat(sprintf("  Timestamp: %s\n", prov[[i]]$timestamp))
    cat(sprintf("  Git commit: %s\n", prov[[i]]$git_commit))
    cat(sprintf("  Branch: %s | Tag: %s | Dirty: %s\n",
                prov[[i]]$git_branch, prov[[i]]$git_tag, prov[[i]]$git_dirty))
    if (!is.null(prov[[i]]$params)) {
      cat("  Params:\n")
      print(prov[[i]]$params)
    }
    cat("\n")
  }
}

get_rds_path <- function(sample, milestone) {
  file.path(rdsdir, paste0(sample, "_", milestone, ".rds"))
}

#' Update trimming_settings for a specific sample.
#'
#' @param settings_df The trimming_settings data.frame.
#' @param sample_name The sample name as stored in the 'sample' column.
#' @param ... Named arguments with fields to update (e.g., nCAlt = 30000).
#' @return Updated settings_df.
#' @examples
#' trimming_settings <- update_trimming_settings(trimming_settings, "Sample1", nCAlt = 25000, nCRlt = 7000)
update_trimming_settings <- function(settings_df, sample_name, ...) {
  # Check sample column exists
  if (!"sample" %in% colnames(settings_df)) stop("settings_df must have a 'sample' column.")
  args <- list(...)
  # Only allow updating existing columns (except sample)
  valid_cols <- setdiff(colnames(settings_df), "sample")
  wrong_cols <- setdiff(names(args), valid_cols)
  if (length(wrong_cols)) stop("Unknown columns in update: ", paste(wrong_cols, collapse=", "))
  # Find row index for this sample
  idx <- which(settings_df$sample == sample_name)
  # If not present, add a new row
  if (length(idx) == 0) {
    new_row <- as.list(rep(NA, length(colnames(settings_df))))
    names(new_row) <- colnames(settings_df)
    new_row$sample <- sample_name
    for (nm in names(args)) new_row[[nm]] <- args[[nm]]
    settings_df <- rbind(settings_df, as.data.frame(new_row, stringsAsFactors = FALSE))
    message(sprintf("Added new entry for sample '%s': %s", sample_name, 
                    paste(sprintf("%s=%s", names(args), args), collapse=", ")))
  } else {
    # Update in place
    old_vals <- as.list(settings_df[idx, names(args), drop = FALSE])
    for (nm in names(args)) settings_df[idx, nm] <- args[[nm]]
    message(sprintf("Updated sample '%s': %s", sample_name, 
                    paste(sprintf("%s: %s -> %s", names(args), old_vals, args), collapse=", ")))
  }
  rownames(settings_df) <- NULL
  return(settings_df)
}

#this just computes them mean percent of mitochondrial genes in an entire sample
meanMT <- function(samplename) {
  return(mean(samplename@meta.data[["percent.mt"]]))
}