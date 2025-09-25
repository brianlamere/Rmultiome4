library(qlcMatrix)
library(BSgenome.Hsapiens.UCSC.hg38)

# 1. Compute gene activity matrix if you want to visualize/compare gene activity (optional)
# note, don't bother unless you want to compare inferred gene activity to
# measured RNA, to check concordance. 

DefaultAssay(tagged_obj) <- "ATAC"
activity_mat <- GeneActivity(tagged_obj)
tagged_obj[["ACTIVITY"]] <- CreateAssayObject(counts = activity_mat)

# 2. (Optional) Normalize gene activity for visualization/comparison
tagged_obj <- NormalizeData(tagged_obj, assay = "ACTIVITY")

tagged_obj <- RegionStats(tagged_obj,
                          genome = BSgenome.Hsapiens.UCSC.hg38,
                          assay = "ATAC")

# 3. Set default assay back to ATAC for LinkPeaks
DefaultAssay(tagged_obj) <- "ATAC"

# 4. Run LinkPeaks
tagged_obj <- LinkPeaks(
  object = tagged_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",      # or "ACTIVITY" if you want to use gene activity
  expression.slot = "data",      # "data" is log-normalized by default
  distance = 10e3                # or whatever window you want (10kb is default)
)