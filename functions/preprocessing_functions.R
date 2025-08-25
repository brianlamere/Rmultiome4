#Functions that are part of preprocessing; nothing here should, by definition
#process the seurat5 object, merely populate it with the data to be processed.

#we're using v86, though v114 is available it is from this month.
load_annotations <- function(ensdb = EnsDb.Hsapiens.v114) {
  annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  return(annotation)
}

#' Create the initial base seurat object
#'
#' @param samplename raw file from Read10X_h5
#' @return base Seurat5 object with percent.mt added
base_object <- function(samplename) {
  fullrna <- paste(rawdatadir, samplename, h5filename, sep = "/")
  fullatac <- paste(rawdatadir, samplename, atacfilename, sep = "/")
  counts <- Read10X_h5(filename = fullrna)
  rna_counts <- counts$`Gene Expression`
  atac_counts <- counts$Peaks
  print("Creating the RNA assay for the Seurat object...")
  baseSeuratObj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = samplename
  )
  print("Adding he ATAC assay for the Seurat object...")
  baseSeuratObj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fullatac,
    annotation = EnsDbAnnos
  )
  print("Calculating a slot for percent.mt for downstream QC")
  DefaultAssay(baseSeuratObj) <- "RNA"
  baseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(baseSeuratObj,
                                                        pattern = "^MT-")     
  return(baseSeuratObj)
}

#' Add chromosome mapping to RNA assay meta.features of a Seurat object
#'
#' @param seurat_obj A Seurat object with an RNA assay
#' @param assay Name of the RNA assay (default is "RNA")
#' @param edb EnsDb annotation database (default is EnsDb.Hsapiens.v86)
#' @return Seurat object with chromosome info added to @meta.features
#' @examples
#' seurat_obj <- chromosome_mapping(seurat_obj)
chromosome_mapping <- function(seurat_obj, assay = "RNA", edb = NULL, warn_threshold = 16000) {
  if (is.null(edb)) {
    if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      stop("Please install EnsDb.Hsapiens.v86: BiocManager::install('EnsDb.Hsapiens.v86')")
    }
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  }
  if (!requireNamespace("ensembldb", quietly = TRUE)) {
    stop("Please install ensembldb: BiocManager::install('ensembldb')")
  }
  if (!requireNamespace("AnnotationFilter", quietly = TRUE)) {
    stop("Please install AnnotationFilter: BiocManager::install('AnnotationFilter')")
  }
  gene_symbols <- rownames(seurat_obj[[assay]])
  filter_obj <- AnnotationFilter::GeneNameFilter(gene_symbols)
  chrom_map <- ensembldb::genes(
    edb,
    filter = filter_obj,
    return.type = "DataFrame"
  )[, c("gene_id", "gene_name", "seq_name")]
  chrom_vec <- chrom_map$seq_name[match(gene_symbols, chrom_map$gene_name)]
  names(chrom_vec) <- gene_symbols
  
  n_mapped <- sum(!is.na(chrom_vec))
  n_total <- length(chrom_vec)
  
  # Only warn if fewer than warn_threshold features are mapped
  if (n_mapped < warn_threshold) {
    warning(sprintf("chromosome_mapping: Only %d/%d features could be mapped to a chromosome (expected 18,000-24,000).", n_mapped, n_total))
  } else {
    message(sprintf("chromosome_mapping: %d/%d features mapped to a chromosome.", n_mapped, n_total))
  }
  
  # For Seurat v5: update feature.info in @misc
  feature_info <- seurat_obj[[assay]]@misc$feature.info
  feature_info$chromosome <- chrom_vec
  seurat_obj[[assay]]@misc$feature.info <- feature_info
  seurat_obj
}

#' Remove features not on standard chromosomes from RNA and ATAC assays in a Seurat object.
#'
#' @param seurat_obj A Seurat object with RNA and/or ATAC assays.
#' @param standard_chroms Character vector of chromosomes to keep (default: chr1-22, X, Y).
#' @return The filtered Seurat object with only standard chromosome features.
remove_nonstandard_chromosomes <- function(seurat_obj, standard_chroms = paste0("chr", c(1:22, "X", "Y"))) {
  # RNA
  if (!is.null(seurat_obj[["RNA"]])) {
    rna_features <- rownames(seurat_obj[["RNA"]])
    rna_chroms <- seurat_obj[["RNA"]]@meta.features$chromosome
    keep_rna <- rna_features[rna_chroms %in% standard_chroms]
    message(
      sprintf(
        "RNA: Keeping %d of %d features on standard chromosomes.",
        length(keep_rna), length(rna_features)
      )
    )
    seurat_obj[["RNA"]] <- subset(seurat_obj[["RNA"]], features = keep_rna)
  }
  # ATAC
  if (!is.null(seurat_obj[["ATAC"]])) {
    atac_features <- rownames(seurat_obj[["ATAC"]])
    atac_chroms <- sub(":.*", "", atac_features)
    keep_atac <- atac_features[atac_chroms %in% standard_chroms]
    message(
      sprintf(
        "ATAC: Keeping %d of %d features on standard chromosomes.",
        length(keep_atac), length(atac_features)
      )
    )
    seurat_obj[["ATAC"]] <- subset(seurat_obj[["ATAC"]], features = keep_atac)
  }
  return(seurat_obj)
}

