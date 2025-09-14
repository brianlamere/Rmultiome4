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

#post-merge RNA modality
DefaultAssay(merged_data) <- "RNA"
merged_data <- FindVariableFeatures(merged_data, assay = "RNA")
merged_data <- ScaleData(merged_data, assay = "RNA",
                         features = VariableFeatures(merged_data))
merged_data <- RunPCA(merged_data, assay = "RNA",
                      features = VariableFeatures(merged_data))
saveRDS(merged_data, "/projects/opioid/vault/postRNA.rds")

#post-merge ATAC modality
DefaultAssay(merged_data) <- "ATAC"
merged_data <- RunTFIDF(merged_data)
merged_data <- FindTopFeatures(merged_data, min.cutoff = 'q0')
merged_data <- RunSVD(merged_data)
saveRDS(merged_data, "/projects/opioid/vault/postATAC.rds")
merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")

#time for harmony
DefaultAssay(merged_data) <- "RNA"
merged_data <- RunHarmony(
  merged_data,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  max_iter = 50,
  reduction.save = "harmony"
)

DefaultAssay(merged_data) <- "ATAC"
merged_data <- RunHarmony(
  object = merged_data,
  group.by.vars = "orig.ident",
  reduction.use = "lsi",
  project.dim = FALSE
)

merged_data <- FindMultiModalNeighbors(
  object = merged_data,
  reduction.list = list("pca", "harmony"),
  dims.list = list(1:50, 1:50)
)

merged_data <- FindClusters(
  merged_data,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.7
)

merged_data <- RunUMAP(
  merged_data,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

saveRDS(merged_data, "/projects/opioid/vault/pre_mapping.rds")

# unless something above has been changed, the below plot now shows pretty good
#integration, with donors spread across the umap
DimPlot(merged_data, reduction = "wnn.umap", group.by = "orig.ident")

#now we need to do the work of mapping cell types so we can then cluster by
#cell types.

#compare these results to PanglaoDB, CellMarker, etc
markers <- FindAllMarkers(merged_data, assay = "RNA", only.pos = TRUE)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(top_markers)

saveRDS(markers, "/projects/opioid/vault/markers.rds")

library(SeuratDisk)
SaveH5Seurat(merged_data, filename = "merged_data.h5Seurat")
Convert("merged_data.h5Seurat", dest = "h5ad")



#Link peaks to nearby genes, or use motif analysis (e.g. Signac::FindMotifs()) 
# to infer cell type-specific TFs.
#atac_markers <- FindAllMarkers(merged_data, assay = "ATAC", only.pos = TRUE)

