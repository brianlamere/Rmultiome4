# Rmultiome
A small tool set that is being designed for 10X multiome data from brain nuclei

Rebased after learning many of the vignettes that seem to apply, don't really
for my particular data set (very high depth, brain nuclei, thus very strong 
batch effects but no need to do things like redoing peaks, etc).

A project goal is being able to do most or all of the steps on my local fairly
substantial PC (192GB RAM, i9-13900k CPU, 8TB nvme raid) and anything that
might need more than that, I have a x2eidn AWS instance with 512GB RAM.

## Configuration

Project settings are stored in `settings.R`.  
Edit this file to match your local paths and data files before running any
scripts.

**Example:**
```r
sourcedir <- "/path/to/data"
h5filename <- "rna_data.h5"
atacfilename <- "atac_data.h5"
referencedir <- "/path/to/data"
referencefile <- "referencefile.rds"
scratchdir1 <- "/path/to/scratch/dir"
```

## Configuration

Edit `settings.R` to specify the path to your Cell Ranger output file.
By default, this is `"filtered_feature_bc_matrix.h5"`, but you may need to
change this if your file is named differently 
(e.g., "DonorA_filtered_feature_bc_matrix.h5").

### RDS Milestone Files

Milestone files are saved as `rds/{sample}_{milestone}.rds`.
This allows you to revert or restart analysis at any stage.

Example:
- `sampleA_qc.rds`
- `sampleA_filtered.rds`
- `sampleA_norm.rds`