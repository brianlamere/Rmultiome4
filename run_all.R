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
  
  # STEP 2: TRIM
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
  
  # STEP 3: preRNA
  preRNA_path <- get_rds_path(sample, "preRNA")
  if (!file.exists(preRNA_path)) {
    obj <- trim_obj
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = 2500)
    obj <- update_provenance(obj, "pre-merge_rna")
    saveRDS(obj, preRNA_path)
  } else {
    obj <- readRDS(preRNA_path)
  }
  
  # STEP 4: preATAC/pipeline1
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

merged_data <- readRDS("/projects/opioid/vault/merged_samples.rds")

# Collapse layers to get unified counts and data
merged_data[["RNA"]] <- JoinLayers(merged_data[["RNA"]])
#merged_data[["ATAC"]] <- JoinLayers(merged_data[["ATAC"]])

# Get list of layers for ATAC
atac_layers <- Layers(merged_data[["ATAC"]])

# Extract count matrices from each layer and combine
counts_list <- lapply(atac_layers, function(lyr) {
  GetAssayData(merged_data[["ATAC"]], layer = lyr, slot = "counts")
})

# Combine matrices by columns (cells)
combined_counts <- do.call(cbind, counts_list)

# Create new ChromatinAssay with the collapsed matrix
new_atac <- CreateChromatinAssay(counts = combined_counts,
                                 genome = merged_data[["ATAC"]]@genome)
merged_data[["ATAC"]] <- new_atac

#post-merge RNA modality
merged_data <- FindVariableFeatures(merged_data, assay = "RNA")
print(length(VariableFeatures(merged_data)))
merged_data <- ScaleData(merged_data, assay = "RNA",
                         features = VariableFeatures(merged_data))
merged_data <- RunPCA(merged_data, assay = "RNA",
                      features = VariableFeatures(merged_data))
saveRDS(postRNA, "/projects/opioid/vault/postRNA.rds")

#post-merge ATAC modality
postATAC <- 
