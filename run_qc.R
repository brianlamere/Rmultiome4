#run qc tests
source("/projects/opioid/Rmultiome/Rmultiome-main.R")
source("/projects/opioid/Rmultiome/functions/qc_functions.R")
source("/projects/opioid/Rmultiome/functions/helper_functions.R")
source("/projects/opioid/Rmultiome/settings.R")
source("/projects/opioid/Rmultiome/functions/trimming_functions.R")


#when settings are changed in settings.R, you'll need to save the file the re-source it

rm(list = ls())
rm(list= c("qcSample", "trimmed_sample"))

sample = "LG05"

qcSample <- readRDS(get_rds_path(sample, "base"))

DefaultAssay(qcSample) <- "ATAC"
qcSample <- NucleosomeSignal(qcSample)
qcSample <- TSSEnrichment(qcSample)

#QC Vln plots against ATAC assay
QCVlnA(qcSample)

QCVlnR(qcSample)

#QC density plots against ATAC assay
QCDensity_ATAC(qcSample)
QCDensity_RNA(qcSample)

source("/projects/opioid/Rmultiome/settings.R")
trimmed_sample <- trimSample(qcSample)

QCDensity_ATAC(trimmed_sample)
QCDensity_RNA(trimmed_sample)
QCVlnA(trimmed_sample)
QCVlnR(trimmed_sample)

plotKDETrim(trimmed_sample, "nCount_ATAC", "TSS.enrichment", kde_percentile = 0.95, 
            retained_cells = colnames(your_trimmed_obj), main = "KDE Trim ATAC")
