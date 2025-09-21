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

cluster_data <- function(harmony_obj, alg = 3, res = 0.5, cluster_dims = 2:30) {
  #clustering. change resolution to increase/decrease number of clusters
  DefaultAssay(harmony_obj) <- "RNA"
  harmony_obj <- FindClusters(
    harmony_obj,
    graph.name = "wsnn",
    algorithm = alg,
    resolution = res,
    group.singletons = FALSE
  )

  harmony_obj <- RunUMAP(
    harmony_obj,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    dims = cluster_dims
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
