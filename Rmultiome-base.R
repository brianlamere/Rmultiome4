#base library for other parts of the toolset
#copyright under LGPL
#reference the git repo https://github.com/brianlamere/Rmultiome
#I sincerely hope others find continued value in anything here for
# their own projects!  Peace.

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)
library(ggplot2)
library(data.table)
library(DoubletFinder)
library(AnnotationHub)

#vignette references:
# https://stuartlab.org/signac/articles/pbmc_multiomic
# https://satijalab.org/seurat/articles/integration_introduction
# https://satijalab.org/seurat/articles/seurat5_integration_bridge

#Doublelet section came from here, but it currently gives inconsistent results
#so will remove unless bug is fixed https://github.com/lzillich/CN_multiome_cocaine

#same annotations are used for all samples.  I set a default value only so 
#the function window shows what annotations we're using, modify as needed
loadannotations <- function(ensdb = EnsDb.Hsapiens.v86) {
  annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  return(annotation)
}

#different method used by lzillich
loadannotations2 <- function(ensdb = EnsDb.Hsapiens.v86, style = "UCSC") {
  ah <- AnnotationHub()
  ensdbs <- query(ah, c("EnsDb.Hsapiens"))
  ensdb_id <- ensdbs$ah_id[grep(paste0(" 105 EnsDb"), ensdbs$title)]
  ensdb <- ensdbs[[ensdb_id]]
  seqlevelsStyle(ensdb) <- style
  annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
  genome(annotations) <- "hg38"
  return(annotations)
}

#I don't understand why so much is being done to the reference object, but
#that is what was in my inherited code
loadreference <- function(reference = "allen_m1c_2019_ssv4.rds") {
  refpath <- paste("your/path/to/references", reference, sep = "/")
  Reductions(reference)
  reference <- SCTransform(reference, verbose = FALSE)
  reference <- RunPCA(reference, npcs = 50)
  reference <- UpdateSeuratObject(reference)
}

#function for creating the base objects from which all the raw data comes
baseobjects <- function(samplename) {
  #now that we know the middle, we can construct full file names
  fullrna <- paste(sourcedir, samplename, h5filename, sep = "/")
  fullatac <- paste(sourcedir, samplename, atacfilename, sep = "/")
  counts <- Read10X_h5(filename = fullrna)
  rna_counts <- counts$`Gene Expression`
  atac_counts <- counts$Peaks
  baseSeuratObj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = samplename
  )
  baseSeuratObj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fullatac,
    annotation = annotations
  )
  #adding the percentage of mitochondrial genes
  DefaultAssay(baseSeuratObj) <- "RNA"
  baseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(baseSeuratObj,
                                                        pattern = "^MT-")
  #these worry me since they become outdated once object is modified
  DefaultAssay(baseSeuratObj) <- "ATAC"
  baseSeuratObj <- NucleosomeSignal(baseSeuratObj)
  baseSeuratObj <- TSSEnrichment(baseSeuratObj)
  return(baseSeuratObj)
}

#this just computes the mean percent of mitochondrial genes in an entire sample
#My goal was to have a longer comment than characters of code in the section
meanMT <- function(samplename) {
  return(mean(samplename@meta.data[["percent.mt"]]))
}

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

#function for trimming samples
#it is extrodinarily recommend you use QA to come up with your own values.
trimSample <- function(samplename,
                       nCAlt = 50000, # nCount_RNA is less than this
                       nCAgt = 1000,
                       nCRlt = 6500,
                       nCRgt = 1000,
                       #we were using nuclei, which have much lower TSS generally
                       TSSeslt = 8,
                       TSSesgt = 2,
                       nslt = 2.5,
                       nsgt = 0.5,
                       pMTlt = 10) {
  #you can only cut shorter, not longer.  leave the base object alone.  Don't
  #trim into the object you're trimming
  myobj <- copy(samplename)
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
  #I would rather use the "command" slot in Seurat object class, but am
  #concerned that no one else does...?
  oldname <- samplename@project.name
  myobj_trim@project.name <- paste(oldname, "trimmed", sep = "_")
  return(myobj_trim)
}

#DoubletFinder path with Normalize to ScaleData.
predoubNormObj <- function(samplename,
                           FN_dims = 1:10,
                           umap_dims = 1:10,
                           FVF_nfeatures = 3000) {
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


