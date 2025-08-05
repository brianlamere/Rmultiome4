#Functions derived substantially from the following project:
# https://github.com/lzillich/CN_multiome_cocaine

library(BSgenome.Hsapiens.UCSC.hg38)

scratchdir1 = "/path/to/scratch"

#DoubletFinder path with Normalize to ScaleData.
predoubNormObj <- function(samplename,
                           FN_dims = 1:10,
                           umap_dims = 1:10,
                           FVF_nfeatures = 3000) {
  #consider adding flag that says SCT version was run and stopping
  myobj <- copy(samplename)
  DefaultAssay(myobj) <- "RNA"
  myobj <- NormalizeData(myobj)
  myobj <- FindVariableFeatures(myobj, selection.method = "vst",
                                nfeatures = FVF_nfeatures)
  myobj <- ScaleData(myobj)
  myobj <- RunPCA(myobj)
  myobj <- FindNeighbors(myobj, reduction = "pca", dims = FN_dims)
  myobj <- FindClusters(myobj, resolution = 0.2)
  myobj <- RunUMAP(myobj, dims = umap_dims)
  UMAPPlot(myobj, reduction = "umap")
  return(myobj)
}

#DoubletFinder path for pK indentification with no "ground truth"
predoubpKI <- function(samplename, sweep_PCs = 1:10) {
  DefaultAssay(samplename)<-"RNA"
  myobj <- copy(samplename)
  sweep.res.list_myobj <- paramSweep(myobj, PCs = sweep_PCs, sct = FALSE)
  sweep.stats_myobj <- summarizeSweep(sweep.res.list_myobj, GT = FALSE)
  bcmvn_myobj <- find.pK(sweep.stats_myobj)
  pk.factors <- bcmvn_myobj[which.max(bcmvn_myobj$BCmetric), "pK"]
  optimal.pk <- as.numeric(as.character(pk.factors))
  print(paste("Optimal pK value is:", optimal.pk))
  return(optimal.pk)
  return(pk.factors)
}

## Homotypic Doublet Proportion Estimate
HDPEobj <- function(samplename, optimal.pk, assumed_doublets = 0.05,
                    dFPCs = 1:10, dFpN = 0.25) {
  myobj <- copy(samplename)
  annotations <- myobj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  nExp_poi <- round(assumed_doublets * nrow(myobj@meta.data))  
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  myobj <- doubletFinder(myobj, PCs = dFPCs, pN = dFpN,
                         pK = optimal.pk, nExp = nExp_poi.adj,
                         sct = FALSE)
  return(myobj)
}

redoPeaks <- function(seuratCombinedObj, sample.list, peakwidthlt = 10000,
                      peakwidthgt = 20, scratchdir = scratchdir1) {
  samples <- sample.list
  seuCmb <- copy(seuratCombinedObj)
  ranges_list <- list()
  for (i in seq_along(samples)) {
    ranges_list[[i]] <- samples[[i]]@assays$ATAC@ranges
  }
  # Combine into a GRangesList, unlist, and reduce
  peaks <- reduce(unlist(as(ranges_list, "GRangesList")))
  #return(peaks)

  peakwidths <- width(peaks)
  peaks <- peaks[peakwidths < peakwidthlt & peakwidths > peakwidthgt]

  counts_atac_merged <- FeatureMatrix(seuCmb@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(seurat))

  seuCmb[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                              fragments = seuCmb@assays$ATAC@fragments,
                              annotation = seuCmb@assays$ATAC@annotation,
                              sep = c(":","-"),
                              genome = "hg38")

  options(future.globals.maxSize = 20000 * 1024^2)

  peaks <- CallPeaks(seuCmb,
                   assay="ATAC",
                   group.by="orig.ident",
                   macs2.path="/usr/bin/macs3",
                   outdir = scratchdir,
                   fragment.tempdir = scratchdir)

  counts_atac <- FeatureMatrix(seuCmb@assays$ATAC@fragments,
                             features = peaks,
                             cells = colnames(seuCmb))

  seuCmb[["ATAC"]] <- CreateChromatinAssay(counts_atac,
                            fragments = seuCmb@assays$ATAC@fragments,
                            annotation = seuCmb@assays$ATAC@annotation,
                            genome = 'hg38')
  
  standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
  
  idx_standard_chroms <- 
    which(as.character(seqnames(granges(seuCmb[["ATAC"]]))) %in% standard_chroms)
  
  seuCmb[["ATAC"]] <- 
    subset(seuCmb[["ATAC"]],
               features = rownames(seuCmb[["ATAC"]])[idx_standard_chroms])
  seqlevels(seuCmb[['ATAC']]@ranges) <- 
    intersect(seqlevels(granges(seuCmb[['ATAC']])),
              unique(seqnames(granges(seuCmb[['ATAC']]))))
  return(seuCmb)
}
