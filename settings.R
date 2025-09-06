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

#milestones that will trigger new save files, etc
milestones <- c("raw", "qc", "filtered", "norm", "clustered")

trimming_settings <- data.frame(
  #sample names
  sample = c(
    "LG05", "LG08", "LG22", "LG23", "LG25", "LG26",
    "LG31", "LG33", "LG38", "LG300", "LG301"
  ),
  #nCount ATAC, less than
  nCAlt = c(
    40000, 30000, 20000, 60000, 60000, 40000,
    50000, 60000, 60000, 50000, 60000
  ),
  #nCount RNA, less than
  nCRlt = c(
    30000, 25000, 30000, 15000, 20000, 20000,
    20000, 40000, 70000, 70000, 50000
  ),
  #nCount ATAC, greater than
  nCAgt = c(
    1900, 900, 200, 2800, 2700, 1200,
    3400, 1200, 500, 400, 600
  ),
  #nCount RNA, greater than
  nCRgt = rep(1000, 11),
  #nucleosome signal, less than
  nslt = c(
    2, 2, 2, 2, 2, 3,
    2, 2, 2.2, 2, 2
  ),
  #nucleasome signal, greater than
  TSSesgt = rep(1, 11)
)
