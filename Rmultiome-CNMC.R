#Functions derived substantially from the following project:
# https://github.com/lzillich/CN_multiome_cocaine

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


postMergeATACProcObj1 <- function(samplename, FTFmin.cutoff = 50) {
  myobj <- copy(samplename)
  DefaultAssay(myobj) <- "ATAC"
  #what is method
  myobj <- RunTFIDF(myobj, method = 1)
  myobj <- FindTopFeatures(myobj, min.cutoff = FTFmin.cutoff)

  print("Now starting the RunSVD step.")
  myobj <- RunSVD(myobj, n = 50)
  print("Now starting the RunUMAP step.")
  myobj <- RunUMAP(myobj,
                   reduction = "lsi",
                   dims = 2:30,
                   reduction.name = "umap_atac",
                   reduction.key = "UMAPATAC_")
  return(myobj)
}

postMergeATACProcObj2 <- function(samplename) {
  myobj <- copy(samplename)
  #breaks here, crashes in next step.  "centering and scaling data matrix"
  #dies with "invalid class "chromatinassay" object:
  print("Now starting the ScaleData step.")
  myobj <- ScaleData(myobj, block.size=500)
  print("Now starting the RunHarmony step.")
  myobj <- RunHarmony(myobj,
                      group.by.vars = "orig.ident",
                      reduction.use = "lsi",
                      dims.use = 2:30,
                      max.iter.harmony = 50,
                      assay.use = "ATAC",
                      reduction.save = "harmony_atac")
  return(myobj)
}
