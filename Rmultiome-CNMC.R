#Functions derived substantially from the following project:
# https://github.com/lzillich/CN_multiome_cocaine

scratchdir1 = "/your/path/to/scratch"

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

#this comes from CN_multiome_cocaine/2_merge_call_peaks.Rmd largely, with
#some mods and flow changes, but the individual steps were substantially theirs
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

  #the "official" method for using macs3 is to just give seurat the macs3
  #path and claim it is macs2.
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
