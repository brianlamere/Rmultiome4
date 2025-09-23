library(SeuratObject, lib.loc = "/projects/Seurat4")
library(Seurat, lib.loc = "/projects/Seurat4")
library(Signac, lib.loc = "/projects/Seurat4")
library(ggplot2, lib.loc = "/projects/Seurat4")
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
library(SeuratDisk)
library(dplyr)
#library(ggplot2)
library(data.table)
library(harmony, lib.loc = "/projects/Seurat4")
library(future)
library(MASS)
source("/projects/opioid/Rmultiome/functions/preprocessing_functions.R")
source("/projects/opioid/Rmultiome/functions/postprocessing_functions.R")
source("/projects/opioid/Rmultiome/functions/trimming_functions.R")
source("/projects/opioid/Rmultiome/functions/helper_functions.R")
source("/projects/opioid/Rmultiome/functions/DE_functions.R")

if (FALSE) {
  #dir.create("/projects/Seurat4", showWarnings = FALSE)
  install.packages("remotes") # if not already installed
  remotes::install_github("satijalab/seurat@v4.3.0",
                          lib = "/projects/Seurat4")
  remotes::install_github("satijalab/seurat-object@v4.1.3",
                          lib = "/projects/Seurat4")
  remotes::install_github("mojaveazure/seurat-disk",
                          lib = "/projects/Seurat4")
  install.packages("/projects/scratch/signac-1.11.0.tar.gz",
                   repos = NULL, type = "source",
                   lib = "/projects/Seurat4")
  install.packages("/projects/scratch/harmony_1.1.0.tar.gz",
                   repos = NULL, type = "source",
                   lib = "/projects/Seurat4")
  BiocManager::install("HGNChelper")
}
