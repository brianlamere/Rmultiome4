#manual mapping
FeaturePlot(merged_data4, features = c("ST18", "MOBP", "TMEM144",
                                       "CTNNA3", "KCNH8")) +
  ggtitle("Feature Algorith = SLM, resolution = 0.4")
FeaturePlot(merged_data, features = c("CBLN2", "MEG3", "CUX2", "HS6ST3", "CA10"))
VlnPlot(merged_data, features = c("ST18", "MOBP", "TMEM144", "CTNNA3", "KCNH8"), group.by = "seurat_clusters")
VlnPlot(merged_data, features = c("ST18"), group.by = "seurat_clusters")
VlnPlot(merged_data, features = "ST18", group.by = "seurat_clusters", pt.size = 0)
