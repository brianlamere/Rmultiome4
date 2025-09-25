
#this isn't a set of functions.  It should be and will be, but currently (alpha) isn't

# Load libraries
library(Signac)
library(Seurat)

# Directory containing your CSVs
da_dir <- "/projects/opioid/DiffAccess/"

# List CSV files
csv_files <- list.files(da_dir, pattern = "*.csv", full.names = TRUE)

# Get peak-gene links from the Seurat object
peak_gene_links <- Links(tagged_obj) # returns a GRanges object

# Create a results list to store linked genes per comparison
linked_genes_list <- list()

for (csv in csv_files) {
  # Extract celltype and comparison from filename
  fname <- basename(csv)
  # Remove extension, parse celltype and comparison
  name <- gsub("^DiffAccess_results_|\\.csv$", "", fname)
  
  # Load DA results
  da_res <- read.csv(csv, row.names = 1)
  
  # Get significant DA peaks
  sig_peaks <- rownames(da_res[da_res$p_val_adj < 0.05, ])
  
  # Subset links by peak column
  da_links <- links[links$peak %in% sig_peaks]
  
  # Get unique linked genes
  linked_genes <- unique(da_links$gene)
  
  # Save results
  linked_genes_list[[name]] <- linked_genes
  
  # Optionally: write all peak-gene links for these DA peaks to CSV
  out_csv <- file.path(da_dir, paste0("peak_gene_links_", name, ".csv"))
  write.csv(as.data.frame(da_links), out_csv, row.names = FALSE)
}

# Example: access linked genes for a specific comparison
# linked_genes_list[["Astrocyte_Low_vs_Acute"]]
