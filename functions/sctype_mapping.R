source("/projects/opioid/Rmultiome/Rmultiome-main.R")
source("/projects/opioid/Rmultiome/settings.R")
merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
DefaultAssay(merged_data) <- "RNA"
harmony_240_40_0.05 <- harmonize_both(merged_data, harmony_max_iter = 50)


harmony_240_40_0.05 <- FMMN_task(harmony_240_40_0.05, dims_pca = 2:40, dims_harmony = 2:40, knn = 40)
harmony_240_40_0.05 <- cluster_data(harmony_240_40_0.05, alg = 3, res = 0.05, cluster_dims = 2:40)
saveRDS(harmony_240_40_0.05, "/projects/opioid/vault/pre_mapping_240_40_0.05.rds")


library("HGNChelper")
library("openxlsx")

source("/usr/local/src/sc-type-1.0/R/gene_sets_prepare.R")
source("/usr/local/src/sc-type-1.0/R/sctype_score_.R")

#db_ <- gene_sets_prepare("/usr/local/src/sc-type-1.0/ScTypeDB_full.xlsx", "Brain")
db_ = "/usr/local/src/sc-type-1.0/ScTypeDB_full.xlsx"
tissue <- "Brain"

gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(harmony_240_40_0.05[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(harmony_240_40_0.05[["RNA"]]$scale.data) else
  as.matrix(harmony_240_40_0.05[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. harmony_240_40_0.05[["RNA"]]$counts or harmony_240_40_0.05[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(harmony_240_40_0.05@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(harmony_240_40_0.05@meta.data[harmony_240_40_0.05@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(harmony_240_40_0.05@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

#available as a wrapper
library("HGNChelper")
library("openxlsx")
source("/usr/local/src/sc-type-1.0/R/auto_detect_tissue_type.R")
source("/usr/local/src/sc-type-1.0/R/gene_sets_prepare.R")
source("/usr/local/src/sc-type-1.0/R/sctype_score_.R")
source("/usr/local/src/sc-type-1.0/R/sctype_wrapper.R")
sctype_obj2 <- run_sctype(harmony_240_40_0.05,known_tissue_type="Brain",
                   custom_marker_file="/usr/local/src/sc-type-1.0/ScTypeDB_full.xlsx",
                   name="sctype_classification",plot=TRUE)

