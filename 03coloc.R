#if needed, install packages
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("coloc", quietly = TRUE)) install.packages("coloc")
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")

#load packages
library(data.table)
library(dplyr)
library(coloc)
library(argparse)

#set up argparse
parser <- ArgumentParser()
parser$add_argument("--phecode", help="all of us phenotype ID")
parser$add_argument("--oid", help="mesa phenotype ID")
parser$add_argument("--db_pop", help="mesa population ID")
parser$add_argument("--gwas_pop", help="all of us population ID")

args <- parser$parse_args()

#find bucket
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')

#read in MESA pQTL table
name_of_pqtl_file <- paste0("mesa_", args$db_pop, "_cis_snps.txt")
pqtl_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_pqtl_file, ".gz .")
system(pqtl_command, intern=TRUE)
unzip_command <- paste0("gunzip ", name_of_pqtl_file, ".gz")
system(unzip_command, intern=TRUE)
qtl_data <- fread(name_of_pqtl_file, header=TRUE)
remove_command <- paste0("rm /home/jupyter/", name_of_pqtl_file)
system(remove_command, intern=TRUE)

#read in AoU GWAS table
name_of_gwas_file <- paste0(args$gwas_pop, "_coloc_mesa_", args$phecode, ".tsv")
gwas_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_gwas_file, " .")
system(gwas_command, intern=TRUE)
gwas_data <- fread(name_of_gwas_file, header=TRUE)

#extract phenotype
phenotype <- args$oid
#cat("Processing" , phenotype, "for MESA", args$db_pop, "+ AoU", args$gwas_pop, "\n")

#filter QTL data for current phenotype
qtl_subset <- qtl_data[qtl_data$phenotype_id == phenotype, ]

#find common SNPs for each phenotype
common_variants <- intersect(qtl_subset$variant_id, gwas_data$ID)

if (length(common_variants) > 0) {
  #filter to common SNPs
  qtl_coloc <- qtl_subset[qtl_subset$variant_id %in% common_variants, ]
  gwas_coloc <- gwas_data[gwas_data$ID %in% common_variants, ]
  
  #merge tables
  merged_data <- inner_join(gwas_coloc, qtl_coloc, by = c("ID" = "variant_id"))
  head(merged_data)
  
  pre_filter <- (nrow(merged_data))
  
  #remove duplicate SNPs
  duplicate_snps <- merged_data$ID[duplicated(merged_data$ID)]
  if (length(duplicate_snps) > 0) {
    cat("Found", length(duplicate_snps), "duplicate SNPs. Removing duplicates...\n")
    #keep most significant
    merged_data <- merged_data %>%
      group_by(ID) %>%
      slice_min(pval_nominal, n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  
  post_filter <-(nrow(merged_data))
  
  cat("Pre-filter SNP count: ", pre_filter, " Post-filter SNP count: ", post_filter, "\n")
  
  #mesa pop sizes
  if (args$db_pop == "META"){
    mesa_pop_size <- 2953
  } else if (args$db_pop == "EUR"){
    mesa_pop_size <- 1270
  } else if (args$db_pop == "AFR"){
    mesa_pop_size <- 675
  } else if (args$db_pop == "AMR"){
    mesa_pop_size <- 642
  }
  
  #all of us pop sizes
  aou_pop <- args$gwas_pop
  phecode <- args$phecode
  phe_counts <- list(
    BI_164    = list(META = c(case=34238, control=156661),
                     EUR  = c(case=19393, control= 88565),
                     AFR  = c(case= 8322, control= 34053),
                     AMR  = c(case= 5676, control= 27743)),
    CV_401.1  = list(META = c(case=69064, control=124692),
                     EUR  = c(case=41257, control= 68448),
                     AFR  = c(case=16709, control= 25976),
                     AMR  = c(case= 9692, control= 24388)),
    EM_202.2  = list(META = c(case=31368, control=168259),
                     EUR  = c(case=15439, control= 97907),
                     AFR  = c(case= 8475, control= 35468),
                     AMR  = c(case= 6643, control= 28246)),
    EM_239    = list(META = c(case=65107, control=128113),
                     EUR  = c(case=44523, control= 64192),
                     AFR  = c(case=10458, control= 32846),
                     AMR  = c(case= 8446, control= 25580)),
    GU_582.2  = list(META = c(case=13723, control=188331),
                     EUR  = c(case= 8021, control=106189),
                     AFR  = c(case= 3505, control= 41195),
                     AMR  = c(case= 1924, control= 33696)),
    MS_708    = list(META = c(case=40336, control=150312),
                     EUR  = c(case=27618, control= 79091),
                     AFR  = c(case= 7987, control= 34736),
                     AMR  = c(case= 4140, control= 29873)),
    NS_333.1  = list(META = c(case=26263, control=173183),
                     EUR  = c(case=18066, control= 94271),
                     AFR  = c(case= 4765, control= 39718),
                     AMR  = c(case= 2928, control= 32283)),
    RE_475    = list(META = c(case=24791, control=172922),
                     EUR  = c(case=13976, control= 97730),
                     AFR  = c(case= 6049, control= 37697),
                     AMR  = c(case= 4256, control= 30614))
  )
  
  #assign outputs
  cc <- phe_counts[[phecode]][[aou_pop]]
  aou_cases    <- as.integer(cc["case"])
  aou_pop_size <- as.integer(cc["case"] + cc["control"])
  
  #prepare datasets
  dataset1 <- list(
    beta = merged_data$slope,
    varbeta = merged_data$slope_se^2,
    snp = merged_data$ID,
    type = "quant",
    N = mesa_pop_size,
    MAF = merged_data$af
  )
  
  dataset2 <- list(
    beta = merged_data$BETA,
    varbeta = merged_data$SE^2,
    snp = merged_data$ID,
    type = "cc",
    N = aou_pop_size,
    s= aou_cases,
    sdY = 1
  )
  
  #run coloc
  result <- coloc.abf(dataset1, dataset2)
  
  #view summary
  cat("Results for", phenotype, ":\n")
  print(result$summary)
  cat("\n")
  outfile <- "/home/jupyter/coloc_output.txt"
  
  if (nrow(common) > 0) {
    
    result <- coloc.abf(dataset1, dataset2)
    
    row <- data.frame(
      phecode = args$phecode,
      olink_id = args$oid,
      mesa_pop = args$db_pop,
      aou_pop = args$gwas_pop,
      t(result$summary),
      check.names = FALSE
    )
    
    write.table(
      row,
      file = outfile,
      sep = "\t",
      row.names = FALSE,
      col.names = !file.exists(outfile),
      append = TRUE,
      quote = FALSE
    )
    
  } else {
    cat("No common variants found for", phenotype, "\n")
  }
}
