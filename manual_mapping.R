#manual mapping
FeaturePlot(merged_data, features = c("ST18", "MOBP", "TMEM144", "CTNNA3", "KCNH8"))
VlnPlot(merged_data, features = c("ST18", "MOBP", "TMEM144", "CTNNA3", "KCNH8"), group.by = "seurat_clusters")
VlnPlot(merged_data, features = c("ST18"), group.by = "seurat_clusters")