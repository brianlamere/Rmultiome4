source("/projects/opioid/Rmultiome/Rmultiome-main.R")
source("/projects/opioid/Rmultiome/settings.R")
merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
DefaultAssay(merged_data) <- "RNA"
premap_obj <- harmonize_both(merged_data, harmony_max_iter = 50)


premap_obj <- FMMN_task(premap_obj, dims_pca = 2:40, dims_harmony = 2:40, knn = 40)
premap_obj <- cluster_data(premap_obj, alg = 3, res = 0.05, cluster_dims = 2:40)
saveRDS(premap_obj, "/projects/opioid/vault/pre_mapping_240_40_0.05.rds")
h240k40r004 <- readRDS("/projects/opioid/vault/pre_mapping_dim240_knn40_res0.04.rds")

library("HGNChelper")
library("openxlsx")

source("/usr/local/src/sc-type-1.0/R/gene_sets_prepare.R")
source("/usr/local/src/sc-type-1.0/R/sctype_score_.R")

#db_ <- gene_sets_prepare("/usr/local/src/sc-type-1.0/ScTypeDB_full.xlsx", "Brain")
db_ = "/usr/local/src/sc-type-1.0/ScTypeDB_full.xlsx"
tissue <- "Brain"

gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(premap_obj[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(premap_obj[["RNA"]]$scale.data) else
  as.matrix(premap_obj[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. premap_obj[["RNA"]]$counts or premap_obj[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(premap_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(premap_obj@meta.data[premap_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(premap_obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3],)

#available as a wrapper
library("HGNChelper")
library("openxlsx")
source("/usr/local/src/sc-type-1.0/R/auto_detect_tissue_type.R")
source("/usr/local/src/sc-type-1.0/R/gene_sets_prepare.R")
source("/usr/local/src/sc-type-1.0/R/sctype_score_.R")
source("/usr/local/src/sc-type-1.0/R/sctype_wrapper.R")
sctype_obj3 <- run_sctype(premap_obj,known_tissue_type="Brain",
                   custom_marker_file="/usr/local/src/sc-type-1.0/ScTypeDB_full.xlsx",
                   name="sctype_classification",plot=TRUE)

