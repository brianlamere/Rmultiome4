# Rmultiome
A small tool set that was designed for 10X multiome data from brain nuclei

I was running into some serious issues while trying to use the 
[Satijalab atacseq integration vignette](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)
but found a similar project from [lzillich](https://github.com/lzillich/CN_multiome_cocaine)

Things which were substantially taken from lzillich I put in the
Rmultiome-CNMC.R file, after making them functions/etc.  There were
things I had already done that were pulled from Satijalab or the scientist
who handed this off to me with some code they had written, those are in the
Rmultiome-base.R file, along with my own steps and changing the workflow.
The intent is that one then starts their own project file where they use
the library of functions; I'll make an example ofthat process here in the 
future.
