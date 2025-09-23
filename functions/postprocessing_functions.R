# Contains rna and ATAC modality functions, and harmony/batch effects reductions

#' @param seurat_obj Seurat object.
#' @return Seurat object with FindVariableFeatures, ScaleData, and RunPCA
post_merge_rna <- function(pm_rna_obj) {
  DefaultAssay(pm_rna_obj) <- "RNA"
  pm_rna_obj <- FindVariableFeatures(pm_rna_obj, assay = "RNA")
  pm_rna_obj <- ScaleData(pm_rna_obj, assay = "RNA",
                         features = VariableFeatures(pm_rna_obj))
  pm_rna_obj <- RunPCA(pm_rna_obj, assay = "RNA",
                      features = VariableFeatures(pm_rna_obj))
  return(pm_rna_obj)
}

#' @param seurat_obj Seurat object.
#' @return Chromatin Assay with RunTFIDF, FindTopFeatures, and RunSVD
post_merge_atac <- function(pm_atac_obj) {
  DefaultAssay(pm_atac_obj) <- "ATAC"
  pm_atac_obj <- RunTFIDF(pm_atac_obj)
  pm_atac_obj <- FindTopFeatures(pm_atac_obj, min.cutoff = 'q0')
  pm_atac_obj <- RunSVD(pm_atac_obj)
}

harmonize_both <- function(harmony_obj, harmony_max_iter = 50,
                         harmony_project.dim = FALSE) {
  DefaultAssay(harmony_obj) <- "RNA"
  harmony_obj <- RunHarmony(
    harmony_obj,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    plot_convergence = TRUE,
    max_iter = harmony_max_iter,
    reduction.save = "harmony"
  )
  DefaultAssay(harmony_obj) <- "ATAC"
  harmony_obj <- RunHarmony(
    object = harmony_obj,
    group.by.vars = "orig.ident",
    reduction.use = "lsi",
    project.dim = harmony_project.dim
  )
  return(harmony_obj)
}

FMMN_task <- function(FMMN_obj, dims_pca = 2:30, dims_harmony = 2:30, knn = 40) {
  FMMN_obj <- FindMultiModalNeighbors(
    object = FMMN_obj,
    reduction.list = list("pca", "harmony"),
    k.nn = knn,
    dims.list = list(dims_pca, dims_harmony)
  )
  return(FMMN_obj)
}

cluster_data <- function(
    harmony_obj, 
    assay = "RNA", 
    alg = 3, 
    res = 0.5, 
    cluster_dims = 2:30, 
    singleton_handling = c("merge", "keep", "discard"),
    cluster_seed = 42
) {
  singleton_handling <- match.arg(singleton_handling)
  DefaultAssay(harmony_obj) <- assay
  
  # Determine group.singletons argument for FindClusters
  if (singleton_handling == "merge") {
    group_singletons <- TRUE
  } else {
    group_singletons <- FALSE
  }
  
  if (assay == "RNA" || assay == "wsnn") {
    graph.name <- "wsnn"
    reduction <- "weighted.nn"
    reduction.name <- "wnn.umap"
    reduction.key <- "wnnUMAP_"
    nn.name <- "weighted.nn"
  } else if (assay == "ATAC") {
    graph.name <- "ATAC_snn"  # This assumes you have built this graph with FindNeighbors(reduction = "lsi")
    reduction <- "lsi"
    reduction.name <- "atac.umap"
    reduction.key <- "ATACUMAP_"
    nn.name <- "ATAC_snn"
  } else {
    stop("Unsupported assay for clustering, use RNA or ATAC.")
  }
  
  
  # Clustering step
  harmony_obj <- FindClusters(
    harmony_obj,
    graph.name = graph.name,
    algorithm = alg,
    resolution = res,
    group.singletons = group_singletons,
    #I do not like needing the below, but I had subtypes of same type that
    #had slight instability, that I was merging anyway but would sometimes 
    #merge on their own, changing my number of clusters thus type assignments
    random.seed = cluster_seed
  )
  
  # If discarding singletons, subset them out
  if (singleton_handling == "discard") {
    if ("singleton" %in% levels(harmony_obj$seurat_clusters)) {
      singleton_cells <- WhichCells(harmony_obj, idents = "singleton")
      harmony_obj <- subset(harmony_obj, cells = setdiff(colnames(harmony_obj), singleton_cells))
      # Drop unused cluster level
      harmony_obj$seurat_clusters <- droplevels(harmony_obj$seurat_clusters)
    }
  }
  
  # UMAP step
  harmony_obj <- RunUMAP(
    harmony_obj,
    reduction = reduction,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    dims = cluster_dims,
    nn.name = nn.name
  )
  
  return(harmony_obj)
}


target_markers <- function(harmony_obj, numMarks = 5) {
  DefaultAssay(harmony_obj) <- "RNA"
  cluster_counts <- table(harmony_obj$seurat_clusters)
  non_singleton_clusters <- names(cluster_counts[cluster_counts > 1])
  
  # Subset object to exclude singleton clusters
  seurat_nonsingleton <- subset(harmony_obj, idents = non_singleton_clusters)
  
  # Find markers for the non-singleton clusters only
  markers_nonsingleton <- FindAllMarkers(seurat_nonsingleton, assay = "RNA")
  
  # For top 5 markers per cluster (adjust 'n' as needed)
  top_markers_nonsingleton <- markers_nonsingleton %>% group_by(cluster) %>% 
    top_n(n = 5, wt = avg_log2FC)
  #or with few/no singletons...
  #markers <- FindAllMarkers(merged_data, assay = "RNA", only.pos = TRUE)
  #top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  #print(top_markers)
  return(harmony_obj)
}
