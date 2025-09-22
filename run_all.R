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

#leaving this here in case you want to restart pre-harmony
#merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
DefaultAssay(merged_data) <- "RNA"

#At this point you stop being production, and move to doing a parameter sweep to
#find the right settings to use.  You will do most of the rest of these tasks
#many times if you do it correctly, since you will be using subjective and objective
#assessments of the different settings, including with cell typing, before moving
#back to differential analysis.  Don't return here until you know what to put on
#the next lines for dims, knn, and resolution

#time for harmony and FindMultiModalNeighbors
harmony_obj <- harmonize_both(merged_data, harmony_max_iter = 50)
saveRDS(harmony_obj, "/projects/opioid/vault/harmonized.rds")

rm(harmony_obj)
rm(premap_obj)

harmony_obj <- FMMN_task(harmony_obj, dims_pca = 2:40, dims_harmony = 2:40, knn = 40)
premap_obj <- cluster_data(harmony_obj, alg = 3, res = 0.04,
                           cluster_dims = 2:40, cluster_seed = 1984)
DimPlot(premap_obj,label=T, raster=FALSE)

#adding file name info that corresponds to the resolution used for FindClusters
saveRDS(premap_obj, "/projects/opioid/vault/pre_mapping_dim240_knn40_res0.04_seed1984.rds")

DimPlot(premap_obj,label=T, raster=FALSE)
#premap_obj <- readRDS("/projects/opioid/vault/pre_mapping_dim240_knn40_res0.05.rds")

labeled_obj <- premap_obj
#Note: these assignments are ONLY TRUE if 2:40/40/0.05 is used and only for this data!
cluster_to_celltype <- c(
  "0" = "Oligodendrocyte",
  "1" = "Glutamatergic neurons",
  "2" = "Astrocyte",
  "3" = "Microglia",
  "4" = "GABAergic neurons",
  "5" = "GABAergic neurons",
  "6" = "Oligodendrocyte precursor cells",
  "7" = "Endothelial",
  "8" = "Mature Neurons",
  "9" = "Oligodendrocyte",
  "10" = "Glutamatergic neurons",
  "11" = "Dopaminergic neurons",
  "12" = ""
)

labeled_obj$celltypes <- cluster_to_celltype[as.character(labeled_obj$seurat_clusters)]

sample_metadata <- data.frame(
  orig.ident = c("LG300", "LG301", "LG22", "LG25", "LG38", "LG05", "LG26", "LG31", "LG08", "LG23", "LG33"),
  HIV_status = c("No HIV", "No HIV", "HIV+", "HIV+", "HIV+", "HIV+", "HIV+", "HIV+", "HIV+", "HIV+", "HIV+"),
  Opioid_exposure = c("Low", "Low", "Low", "Low", "Low", "Acute", "Acute", "Acute", "Chronic", "Chronic", "Chronic")
)

meta <- labeled_obj@meta.data %>%
  left_join(sample_metadata, by = "orig.ident")

labeled_obj@meta.data <- meta

#to check
DimPlot(labeled_obj,label=T, raster=FALSE)
DimPlot(labeled_obj, group.by = "celltypes", label = TRUE, raster=FALSE)
table(harmony_240_40_0.05$HIV_status, harmony_240_40_0.05$Opioid_exposure)