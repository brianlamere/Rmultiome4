run_DiffExpress_and_export <- function(
    seurat_obj,
    celltype_col = "celltypes",
    celltype,
    group_col = "group",
    ident.1,
    ident.2,
    output_prefix = "DiffExpress_results"
) {
  # Subset Seurat object by cell type
  cell_subset <- subset(seurat_obj, subset = !!as.name(celltype_col) == celltype)
  DefaultAssay(obj_assigned) <- "RNA"
  # Run DE
  markers <- FindMarkers(
    cell_subset,
    ident.1 = ident.1,
    ident.2 = ident.2,
    assay = "RNA",
    group.by = group_col
  )
  
  # Compose output filename
  fname <- sprintf("%s_%s_%s_vs_%s.csv", output_prefix, celltype, ident.1, ident.2)
  
  # Write full results to CSV
  write.csv(markers, fname)
  
  # Return the DE result (optional)
  invisible(markers)
}

run_DiffAccess_and_export <- function(
    seurat_obj,
    celltype_col = "celltypes",
    celltype,
    group_col = "group",
    ident.1,
    ident.2,
    output_prefix = "DiffAccess_results"
) {
  # Subset Seurat object by cell type
  cell_subset <- subset(seurat_obj, subset = !!as.name(celltype_col) == celltype)
  DefaultAssay(obj_assigned) <- "RNA"
  # Run DE
  markers <- FindMarkers(
    cell_subset,
    ident.1 = ident.1,
    ident.2 = ident.2,
    group.by = group_col,
    assay = "ATAC",
    test.use = "LR"
  )
  
  # Compose output filename
  fname <- sprintf("%s_%s_%s_vs_%s.csv", output_prefix, celltype, ident.1, ident.2)
  
  # Write full results to CSV
  write.csv(markers, fname)
  
  # Return the DE result (optional)
  invisible(markers)
}
