sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)

# input data
f <- snakemake@input[["units"]]
d <- read.csv(f)

cohort = snakemake@params[["cohort"]]
clonalSNVs_dir = snakemake@params[["clonalSNVs_dir"]]

# filter for target cohort
d = d |> filter(sample_type %in% cohort) 

# read data and calculate mean AF
d[c("TF", "TF_uCI", "TF_lCI",                     
    "n_pass", "n_obs")] <- NA

for (row in 1:nrow(d)){ #
    pt = d[row, "sample_ID"][1]
    fd = d[row, "cfDNA_ID"][1]
       
    AF_table <- read.csv(paste0(clonalSNVs_dir, "results/", pt, "/", fd, "/SNVs_in_cfDNA.txt"), sep = "\t")
    
    n_pass <- nrow(AF_table)
    n_obs <- sum(AF_table$cfDNA_n_alt_total > 0)
    AF_table <- AF_table |> mutate(AF = cfDNA_n_alt_total/cfDNA_DP_total)
    mean_AF <- mean(AF_table$AF)
    
    d[row, c("TF", "TF_uCI", "TF_lCI")] = c(mean_AF, NA, NA)
    d[row, c("n_pass", "n_obs")] = c(n_pass, n_obs)
}


write.csv(d, snakemake@output[["results"]], row.names = FALSE)

sink()