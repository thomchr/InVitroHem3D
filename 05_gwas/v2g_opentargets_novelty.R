source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
library(dplyr)
library(arrow)
library(nVennR)

setwd(file.path(DATA_DIR, "v2g"))
# -----------------------------------------------------------------
# 1. Load Your High-Confidence Data
# -----------------------------------------------------------------
my_v2g_data <- readRDS("allv2g.rds")
core = readLines("intersect1908genes.txt")
unique_genes <- unique(my_v2g_data$ensembl_gene_id[which(my_v2g_data$gene_name %in% core)])
# -----------------------------------------------------------------
# 2. Open Connection to Public Database
# -----------------------------------------------------------------
# (This doesn't load data into RAM yet, just creates the connection)
ot_v2g_dataset <- open_dataset(file.path(DATA_DIR, "V2G/open_targets_v2g"))

output_file <- "Novel_SNP_to_Gene_Pairs.core1908.csv"

# If the script was run previously, delete the old file so we start fresh
if (file.exists(output_file)) {
  file.remove(output_file)
}

print(paste("Starting loop: Processing", length(unique_genes), "genes to find novel pairs..."))

# -----------------------------------------------------------------
# 3. Process Gene-by-Gene and Append to File
# -----------------------------------------------------------------
for (i in seq_along(unique_genes)) {
  
  current_gene <- unique_genes[i]
  
  # A. Isolate YOUR data for the current gene
  my_gene_data <- my_v2g_data %>% filter(ensembl_gene_id == current_gene) %>%
    select(ensembl_gene_id, chr, variant_id,variant_pos) %>%
    distinct()
  
  # B. Extract PUBLIC data ONLY for the current gene into RAM
  public_snps <- ot_v2g_dataset %>%
    filter(gene_id == current_gene) %>%
    select(gene_id, chr_id, position) %>%
    collect()
  
  # C. Check for matches and isolate novelty
  if (nrow(public_snps) > 0) {
    # Format public data to match "chrX" style
    public_snps <- public_snps %>%
      mutate(chr_id = paste0("chr", chr_id)) %>%
      distinct()
    
    # Keep only the rows in your data that do NOT match the public data
    novel_pairs <- anti_join(
      my_gene_data,
      public_snps,
      by = c(
        "ensembl_gene_id" = "gene_id", 
        "chr" = "chr_id",    # REPLACE with your SNP chromosome column
        "variant_pos" = "position"   # REPLACE with your SNP position column
      )
    )
  } else {
    # If Open Targets has NO SNPs for this gene, ALL of your SNPs are novel!
    novel_pairs <- my_gene_data
  }
  
  # D. Append results to the final file
  if (nrow(novel_pairs) > 0) {
    write.table(
      novel_pairs, 
      file = output_file, 
      sep = ",", 
      row.names = FALSE, 
      # col.names only writes headers if the file doesn't exist yet (i.e., the first write)
      col.names = !file.exists(output_file), 
      append = TRUE
    )
  }
  
  # Optional: Print progress every 100 genes so you know it's working
  if (i %% 100 == 0) {
    print(paste("Processed", i, "out of", length(unique_genes), "genes..."))
  }
}

print("Finished! All novel pairs have been successfully saved.")