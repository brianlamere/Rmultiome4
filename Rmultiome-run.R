#functions get defined once and fixed when needed, but these
#are usually moving processes.  This is where I actually run
#the tool, but for publishing I will only include a few
#example lines.

source("your/path/here/Rmultiome-base.R")
#this file is for mostly static settings, such as where files are found
source("/your/path/here/Rmultiome-settings.R")

#if you're creating base objects, you'll need to load annotations once
annotations <- loadannotations2()

#now we create the base objects using the baseobjects function above
sample1 <- baseobjects("sample1")
sample2 <- baseobjects("sample2")
sample3 <- baseobjects("sample3")

#run these and other QA checks for every sample
QCVlnA((sample1))
QCVlnR(sample1)
QCDensA(sample1)

#create subset trimmed objects
sample1_trimmed <- trimSample(sample1, nCAlt = 50000, nCRlt = 40000, nCAgt = 1900,
                           nCRgt = 1000, nslt = 2, TSSesgt = 1)
sample2_trimmed <- trimSample(sample2, nCAlt = 30000, nCRlt = 25000, nCAgt = 900,
                           nCRgt = 1000, nslt = 2, TSSesgt = 1)
sample3_trimmed <- trimSample(sample3, nCAlt = 20000, nCRlt = 50000, nCAgt = 200,
                           nCRgt = 1000, nslt = 2.5, TSSesgt = 1)

#Change parameters on trimming lines until you're happy with the trimmed results
#when done this way versus way I've typically seen it, it will use more ram BUT
#will only take seconds to see results of new QA.
QCVlnA((sample1_trimmed))
QCVlnR(sample1_trimmed)
QCDensA(sample1_trimmed)

#prepare objects for further processing
sample1_pdP <- predoubProcObj(sample1_trimmed)
sample2_pdP <- predoubProcObj(sample2_trimmed)
sample3_pdP <- predoubProcObj(sample3_trimmed)

#compute pKI
sample1_pKI <- predoubpKI(sample1_pdN)
sample2_pKI <- predoubpKI(sample2_pdN)
sample3_pKI <- predoubpKI(sample3_pdN)

#create HDPE object
sample1_HDPE <- HDPEobj(sample1_pdN, optimal.pk = sample1_pKI)
sample2_HDPE <- HDPEobj(sample2_pdN, optimal.pk = sample2_pKI)
sample3_HDPE <- HDPEobj(sample3_pdN, optimal.pk = sample3_pKI)

#link peaks between RNA and ATAC
sample1_LP <- linkPeaks(sample1_HDPE)
sample2_LP <- linkPeaks(sample2_HDPE)
sample3_LP <- linkPeaks(sample3_HDPE)

#Merge objects
samplelist <- list(sample1_LP, sample2_LP, sample3_LP)
seurat <- merge(x = samplelist[[1]], y = samplelist[-1])

seurat@meta.data$donor[seurat@meta.data$orig.ident=="sample1"] <- "sample1"
seurat@meta.data$donor[seurat@meta.data$orig.ident=="sample2"] <- "sample2"
seurat@meta.data$donor[seurat@meta.data$orig.ident=="sample3"] <- "sample3"

seuratRNA <- postMergeRNAProcObj(seurat)

#broken into subsections during debugging, likely will re-merge later
seurat2 <- redoPeaks1(seurat, samplelist)
seurat3 <- redoPeaks2(seurat)

seurat4 <- postMergeATACProcObj1(seurat3)
#This is where it programmically fails, but it is due to problems in the data in an earlier unknown step
seurat5 <- postMergeATACProcObj2(seurat4)
