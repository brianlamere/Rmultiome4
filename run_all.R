source("/projects/opioid/Rmultiome/Rmultiome-main.R")
source("/projects/opioid/Rmultiome/settings.R")

standard_chroms <- paste0("chr", c(1:22, "X", "Y"))

missing <- setdiff(samplelist, trimming_settings$sample)
if (length(missing) > 0) {
  stop(
    "These samples are in samplelist but missing in trimming_settings: ",
    paste(missing, collapse = ", ")
  )
}

EnsDbAnnos <- loadannotations()

for (sample in samplelist) {
  cat(sprintf("\nProcessing sample: %s\n", sample))
  
  # STEP 1: RAW (import)
  base_path <- get_rds_path(sample, "base")
  if (!file.exists(base_path)) {
    base_obj <- base_object(sample)
    print("Adding chromosome mapping information to ATAC assay.")
    base_obj <- chromosome_mapping(base_obj, rna_annos = EnsDbAnnos)
    print("Removing non-standard chromosomes from ATAC and RNA.")
    base_obj <- remove_nonstandard_chromosomes(base_obj)
    base_obj <- update_provenance(base_obj, "raw_import")
    saveRDS(base_obj, base_path)
  } else {
    base_obj <- readRDS(base_path)
  }
  
  # STEP 2: 1D TRIM
  trimmed_path <- get_rds_path(sample, "trimmed")
  if (!file.exists(trimmed_path)) {
    trim_obj <- base_obj
    print("Calculating NucleosomeSignal and TSSEnrichment for ATAC data.")
    DefaultAssay(trim_obj) <- "ATAC"
    trim_obj <- NucleosomeSignal(trim_obj)
    trim_obj <- TSSEnrichment(trim_obj)
    print("Trimming based on QC")
    trim_obj <- trimSample(trim_obj)
    trim_obj <- update_provenance(trim_obj, "trim_data")
    saveRDS(trim_obj, trimmed_path)
  } else {
    trim_obj <- readRDS(trimmed_path)
  }
  
  # STEP 3: 2D TRIM
  kde_path <- get_rds_path(sample, "kdetrim")
  if (!file.exists(kde_path)) {
    kde_obj <- trim_obj
    print("Doing n-Dimensional KDE trimming.")
    kde_obj <- kdeTrimSample(kde_obj,
                             atac_percentile = 0.96,
                             rna_percentile = 0.96,
                             combine_method = "intersection" #or "union" if desired
                             )
    saveRDS(kde_obj, kde_path)
  } else {
    kde_obj <- readRDS(kde_path)
  }
    
  # STEP 4: preRNA
  preRNA_path <- get_rds_path(sample, "preRNA")
  if (!file.exists(preRNA_path)) {
    obj <- kde_obj
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = 2500)
    obj <- update_provenance(obj, "pre-merge_rna")
    saveRDS(obj, preRNA_path)
  } else {
    obj <- readRDS(preRNA_path)
  }
  
  # STEP 5: preATAC/pipeline1
  pipeline1_path <- get_rds_path(sample, "pipeline1")
  if (!file.exists(pipeline1_path)) {
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj)
    obj <- update_provenance(obj, "pre-merge_atac")
    saveRDS(obj, pipeline1_path)
  }
  cat(sprintf("Sample %s completed successfully.\n", sample))
}

merged_data <- merge_sample_objects(samplelist)

saveRDS(merged_data, "/projects/opioid/vault/merged_samples.rds")
#merged_data <- readRDS("/projects/opioid/vault/merged_samples.rds")

#post-merge RNA modality
merged_data <- post_merge_rna(merged_data)
saveRDS(merged_data, "/projects/opioid/vault/postRNA.rds")

#post-merge ATAC modality
merged_data <- post_merge_atac(merged_data)
saveRDS(merged_data, "/projects/opioid/vault/postATAC.rds")

#find the elbow to set dimensionality for Neighbors.
findElbow(merged_data)

#leaving this here in case you want to restart pre-harmony
#merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
DefaultAssay(merged_data) <- "RNA"
#time for harmony and FindMultiModalNeighbors
harmony_obj <- harmony_FMMN(merged_data, harmony_max_iter = 50,
                            dims_pca = 1:40, dims_harmony = 1:40)

saveRDS(harmony_obj, "/projects/opioid/vault/harmonized40.rds")

harmony_obj <- cluster_data(harmony_obj, alg = 3, res = 0.05, cluster_dims = 1:40)

#adding file name info that corresponds to the resolution used for FindClusters
saveRDS(harmony_obj, "/projects/opioid/vault/pre_mapping_05_40.rds")

premap_obj <- readRDS("/projects/opioid/vault/pre_mapping.rds")

DefaultAssay(premap_obj) <- "RNA"
DimPlot(harmony_obj, reduction = "wnn.umap", group.by = "orig.ident", raster = FALSE) +
  ggtitle("Algorith = SLM, Resolution = 0.05")
DimPlot(harmony_obj, reduction = "wnn.umap", group.by = "seurat_clusters", raster = FALSE) +
  ggtitle("Algorith = SLM, Resolution = 0.05")

VlnPlot(premap_obj, features = "ST18", group.by = "seurat_clusters", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#markers_nonsingleton4 <- target_markers(merged_data4)
#saveRDS(markers_nonsingleton4, "/projects/opioid/vault/markers4.rds")
#markers_nonsingleton4 <- readRDS("/projects/opioid/vault/markers4.rds")

