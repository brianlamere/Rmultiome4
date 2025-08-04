#functions get defined once and fixed when needed, but these
#are usually moving processes.  This is where I actually run
#the tool, but for publishing I will only include a few
#example lines.

source("your/path/here/Rmultiome-base.R")

sourcedir = "your/path/here/source"
h5filename = "filtered_feature_bc_matrix.h5"
atacfilename = "atac_fragments.tsv.gz"

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
sample1_pdN <- predoubNormObj(sample1_trimmed)
sample2_pdN <- predoubNormObj(sample2_trimmed)
sample3_pdN <- predoubNormObj(sample3_trimmed)

#compute pKI
sample1_pKI <- predoubpKI(sample1_pdN)
sample2_pKI <- predoubpKI(sample2_pdN)
sample3_pKI <- predoubpKI(sample3_pdN)

#create HDPE object
sample1_HDPE <- HDPEobj(sample1_pdN, optimal.pk = sample1_pKI)
sample2_HDPE <- HDPEobj(sample2_pdN, optimal.pk = sample2_pKI)
sample3_HDPE <- HDPEobj(sample3_pdN, optimal.pk = sample3_pKI)
