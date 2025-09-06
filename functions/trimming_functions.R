#trimming functions

#' trimSample: Trims cells and features from a Seurat object using sample-specific thresholds.
#'
#' This version expects the Seurat object to have its sample identity in the @project.name slot.
#' It looks up trimming parameters in the global `trimming_settings` data.frame.
#'
#' @param seurat_obj A Seurat object (with @project.name set to sample name)
#' @param trimming_settings (optional) Data frame with sample-specific trimming parameters.
#'        If not supplied, will look for `trimming_settings` in the global environment.
#' @return A filtered/truncated Seurat object
#' @export
trimSample <- function(seurat_obj, trimming_settings = NULL) {
  # Get sample name from Seurat object
  sample_name <- seurat_obj@project.name
  # Fetch trimming_settings if not provided
  if (is.null(trimming_settings)) {
    if (!exists("trimming_settings", envir = .GlobalEnv)) {
      stop("trimming_settings not found in global environment, and not provided as argument.")
    }
    trimming_settings <- get("trimming_settings", envir = .GlobalEnv)
  }
  # Look up parameters for this sample
  params <- trimming_settings[trimming_settings$sample == sample_name, ]
  if (nrow(params) == 0) {
    stop(sprintf("No trimming settings found for sample '%s'.", sample_name))
  }
  # Extract thresholds, with names matching your trimming_settings
  nCAlt     <- params$nCAlt
  nCRlt     <- params$nCRlt
  nCAgt     <- params$nCAgt
  nCRgt     <- params$nCRgt
  nslt      <- params$nslt
  TSSesgt   <- params$TSSesgt
  
  # --- Example Trimming Logic ---
  # (Replace with your actual trimming/filtering code as appropriate)
  # Here we assume basic filtering on nCount_RNA, nFeature_RNA, etc.
  
  # Ensure required metadata columns exist
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% colnames(seurat_obj@meta.data))) {
    stop("Required metadata columns not found in Seurat object.")
  }
  
  # Example filtering (customize as needed for your assays and QC metrics)
  trimmed <- subset(
    seurat_obj,
    subset = nCount_RNA < nCAlt &
      nCount_RNA > nCAgt &
      nFeature_RNA < nCRlt &
      nFeature_RNA > nCRgt
    # Add further criteria for nslt, TSSesgt, etc., as needed
  )
  
  # Optionally, store trimming info in the misc slot for future provenance
  trimmed@misc$trimming <- list(
    sample = sample_name,
    nCAlt = nCAlt,
    nCRlt = nCRlt,
    nCAgt = nCAgt,
    nCRgt = nCRgt,
    nslt = nslt,
    TSSesgt = TSSesgt,
    timestamp = Sys.time()
  )
  
  return(trimmed)
}