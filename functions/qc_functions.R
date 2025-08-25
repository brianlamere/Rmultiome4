#functions for QC of objects

#QC density plots against ATAC assay
QCDensA <- function(samplename) {
  DefaultAssay(samplename) <- "ATAC"
  DS1sn2 <- DensityScatter(samplename, x = 'nCount_ATAC', y = 'TSS.enrichment',
                           log_x = TRUE, quantiles = TRUE
  ) +
    ggtitle(paste("ATAC+TSS density plot for ", samplename@project.name))
  DS2sn2 <- DensityScatter(samplename, x = 'nCount_RNA', y = 'TSS.enrichment',
                           log_x = TRUE, quantiles = TRUE
  ) +
    ggtitle(paste("RNA+TSS density plot for ", samplename@project.name))
  DS3sn2 <- DensityScatter(samplename, x = 'nCount_ATAC', y = 'nCount_RNA',
                           log_x = TRUE, quantiles = TRUE
  ) +
    ggtitle(paste("ATAC+RNA density plot for ", samplename@project.name))
  print(DS1sn2)
  print(DS2sn2)
  print(DS3sn2)
}

#QC Vln plots against ATAC assay
QCVlnA <- function(samplename) {
  DefaultAssay(samplename) <- "ATAC"
  VPsn2 <- VlnPlot(
    object = samplename,
    features = c("nCount_ATAC", "nCount_RNA", "TSS.enrichment",
                 "nucleosome_signal"),
    ncol = 4,
    pt.size = 0) 
  print(VPsn2)
}

#QC density plots against RNA assay
QCVlnR <- function(samplename) {
  DefaultAssay(samplename) <- "RNA"
  VP2sn2 <- VlnPlot(samplename, features = "percent.mt") +
    ggtitle(paste("MT Density plot for ", samplename@project.name))
  print(VP2sn2)
  print(paste("The mean value of the MT percentages in this sample is ",
              meanMT(samplename)))
}