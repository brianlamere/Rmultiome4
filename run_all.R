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

for (sample in samplelist) {
  # STEP 1: RAW (import)
  base_obj <- base_object(sample)
  base_obj <- chromosome_mapping(base_obj)
  base_obj <- update_provenance(base_obj, "raw_import")
  saveRDS(base_obj, get_rds_path(sample, "base"))
  
  # STEP 2: QC
  #trimming based on nCount
  trimmed_object <- trimSample(base_object)
  
  #seurat_qc <- do_qc(raw_obj) # your QC function
  #saveRDS(seurat_qc, get_rds_path(sample, "qc"))
  
  # STEP 3: Filtering
  #seurat_filtered <- filter_cells(seurat_qc)
  #saveRDS(seurat_filtered, get_rds_path(sample, "filtered"))
  
  # STEP 4: Normalization
  #seurat_norm <- normalize_data(seurat_filtered)
  #saveRDS(seurat_norm, get_rds_path(sample, "norm"))
  
  # STEP 5: Clustering
  #seurat_clustered <- cluster_cells(seurat_norm)
  #saveRDS(seurat_clustered, get_rds_path(sample, "clustered"))
}

