# Edit this file to match your local environment.
#the settings in this file are very specific to installation and project!

#vault of rds files saved at different milestones
rdsdir <- "/projects/opioid/vault"

# Project raw data directory
rawdatadir <- "/projects/opioid/rawdata"

#project source code directory
Rmultiome_files <- "/projects/opioid/Rmultiome"

# Main 10X h5 file (RNA)
h5filename <- "filtered_feature_bc_matrix.h5"

# ATAC fragments file name
atacfilename <- "atac_fragments.tsv.gz"

# Reference files directory
referencedir <- "/projects/opioid/references"

# Reference annotation file
referencefile <- "allen_m1c_2019_ssv4.rds"

# Scratch space for temporary files
scratchdir1 <- "/projects/opioid/scratch"

#list of samples
samplelist <- c("LG05", "LG08", "LG22", "LG23", "LG25", "LG26", "LG31", 
                "LG33", "LG38", "LG300", "LG301")

trimming_settings <- data.frame(
  # Instead of applying 1D filters on nCount or nFeature values individually,
  # I used 2D density plots of nCount vs TSS.enrichment and nFeature vs 
  # TSS.enrichment.  Cells were retained if they fell within regions of the plot
  # with density > 1, regardless of their individual nCount or nFeature values.
  
  #target technical limits:
  #     percent.mt 0 < X < 5
  #     nucleosome_signal 0.5 < X < 2.5
  #     TSS.enrichment 4 < X < 10
  #     nCount_RNA 200 < X < 20000 - lower than cells
  #     nFeature_RNA 100 < X < 5000 - lower than cells
  #     nCount_ATAC 500 < X < 50000 - lower than cells
  #     nFeature_ATAC 100 < X < 25000 - lower than cells
  #sample names
  sample = c(
    "LG05", "LG08", "LG22", "LG23",
    "LG25", "LG26", "LG31", "LG33",
    "LG38", "LG300", "LG301"
  ),
  #nCount ATAC, less than
  max_nCount_ATAC = c(
    40000, 20000, 20000, 60000,
    15000, 40000, 30000, 60000,
    45000, 50000, 60000
  ),
  #nCount RNA, less than
  max_nCount_RNA = c(
    25000, 15000, 25000, 15000,
    20000, 15000, 15000, 40000,
    45000, 45000, 50000
  ),
  #nCount ATAC, greater than
  min_nCount_ATAC = c(
    1000, 500, 400, 3000,
    2750, 1200, 3000, 2100,
    1400, 2100, 1000
  ),
  #nCount RNA, greater than
  min_nCount_RNA = c(
    200, 200, 200, 100,
    100, 100, 100, 100,
    800, 100, 300
  ),
  #nucleosome signal, less than
  max_nss = c(
    1.5, 1.5, 2.5, 1.2,
    1.5, 2.25, 1.75, 1.5,
    1.75, 1.5, 2
  ),
  min_nss = c(
    0.4, 0.4, 0.4, 0.4,
    0.4, 0.4, 0.4, 0.4,
    0.4, 0.4, 0.4
  ),
  #percent of mitochondrial content.  Ideal is 0 for nuclei
  #we are setting a hard limit of 10, but lower if the 90% percentile is less
  max_percentMT = c(
    8, 6, 6, 6,
    8, 8.5, 6, 6,
    2.5, 4.5, 4.5
    ),
  #nucleasome signal, less than
  max_TSS = c(
    10, 9, 9, 9,
    8, 8, 8, 9,
    7.5, 8, 10
    ),
  #nucleasome signal, greater than
  min_TSS = c(
    2.5, 2.5, 1, 2.25,
    2.75, 1, 1, 2.5,
    1, 1, 1
    )
)
