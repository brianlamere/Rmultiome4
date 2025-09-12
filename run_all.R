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

pipeline1_step <- "base"
continue_pipeline <- TRUE

for (sample in samplelist) {
  # STEP 1: RAW (import)
  if (pipeline1_step == "base") {
    # Recompute base_obj from scratch
    base_obj <- base_object(sample)
    print("Adding chromosome mapping information to ATAC assay.")
    base_obj <- chromosome_mapping(base_obj, rna_annos = EnsDbAnnos)
    print("Removing non-standard chromosomes from ATAC and RNA.")
    base_obj <- remove_nonstandard_chromosomes(base_obj)
    base_obj <- update_provenance(base_obj, "raw_import")
    saveRDS(base_obj, get_rds_path(sample, "base"))
    if (continue_pipeline == TRUE) {
      pipeline1_step <- "trim"
    }
  } else if (pipeline1_step == "trim") {
    trim_obj <- readRDS(get_rds_path(sample, "base"))
    print("Calculating NucleosomeSignal and TSSEnrichment for ATAC data.")
    DefaultAssay(trim_obj) <- "ATAC"
    trim_obj <- NucleosomeSignal(trim_obj)
    trim_Obj <- TSSEnrichment(trim_obj)
    print("Trimming based on QC")
    trim_obj <- trimSample(trim_obj)
    trim_obj <- update_provenance(trim_obj, "trim_data")
    saveRDS(trim_obj, get_rds_path(sample, "trimmed"))
    if (continue_pipeline == TRUE) {
      pipeline1_step <- "preRNA"
    }
  } else if (pipeline1_step == "preRNA") {
    # Normalize, find variable features
    obj <- readRDS(get_rds_path(sample, "trimmed"))
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = FVF_nfeatures)
    obj <- update_provenance(obj, "pre-merge_rna")
    saveRDS(obj, get_rds_path(sample, "preRNA"))
    if (continue_pipeline == TRUE) {
      pipeline1_step <- "preATAC"
    }
  } else if (pipeline1_step == "preATAC") {
    obj <- readRDS(get_rds_path(sample, "preRNA"))
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj)
    obj <- update_provenance(obj, "pre-merge_atac")
    # Save as pipeline1 to signal completion of pre-merge steps
    saveRDS(obj, get_rds_path(sample, "pipeline1"))
  }
}  

pipeline2_step == "stop"

{
  if (pipeline2_step == "merge") {
    merge_sample_objects(samplelist)
  } else if (pipeline2_step == "postRNA") {
    #
  } else if (pipeline2_step == "postATAC") {
    #
  } else if (pipeline2_step == "harmony") {
    #
  } else if (pipeline2_step == "integrate") {
    #
  }
}
