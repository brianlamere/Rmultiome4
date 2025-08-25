#trimming functions

#' Subset a Seurat object by QC metrics
#'
#' @param samplename Seurat object to trim
#' @param nCAlt, nCRlt, nCAgt, nCRgt, TSSeslt, TSSesgt, nslt, nsgt, pMTlt Trimming thresholds
#' @return Trimmed Seurat object
trimSample <- function(samplename,
                       nCAlt = 50000,
                       nCAgt = 1000,
                       nCRlt = 20000,
                       nCRgt = 1000,
                       TSSeslt = 8,
                       TSSesgt = 2,
                       nslt = 2.5,
                       nsgt = 0.5,
                       pMTlt = 10) {
  req_cols <- c("nCount_ATAC", "nCount_RNA", "nucleosome_signal",
                "percent.mt", "TSS.enrichment")
  missing_cols <- setdiff(req_cols, colnames(samplename@meta.data))
  if (length(missing_cols) > 0) {
    warning("Missing required columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }
  myobj <- samplename # No need to copy if subsetting returns new object
  myobj_trim <- subset(
    x = myobj,
    subset = nCount_ATAC < nCAlt &
      nCount_RNA < nCRlt &
      nCount_ATAC > nCAgt &
      nCount_RNA > nCRgt &
      nucleosome_signal < nslt &
      nucleosome_signal > nsgt &
      percent.mt < pMTlt &
      TSS.enrichment < TSSeslt &
      TSS.enrichment > TSSesgt
  )
  n_initial <- ncol(myobj)
  n_trimmed <- ncol(myobj_trim)
  if (n_trimmed == 0) {
    warning("All cells were filtered out! Please check thresholds.")
  } else if (n_trimmed < 0.5 * n_initial) {
    warning("More than half of cells were filtered; check trimming parameters.")
  }
  
  oldname <- myobj@project.name
  myobj_trim@project.name <- paste(oldname, "trimmed", sep = "_")
  return(myobj_trim)
}