sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(stringr)

# input data
estimates = snakemake@input[["estimates"]]

# combine count-filtered UT sets TF estimates
res = NULL
for (file in estimates){
    print(file)
    cutoff = str_split(str_split(file, "/", 9)[[1]][7], "ct", 2)[[1]][2]
    d = read.csv(file) # read in TF estimate information

    file_prefix = paste(str_split(file, "/", 9)[[1]][1:3], collapse="/") # find correct path
    UTmdata <- read.table(paste0(file_prefix, "/unique_tumor_count_filter_test/UT_mdata_minct", cutoff, ".txt"), header = TRUE) # UT set size
    
    row <- tibble(tf_mean = d$tf_mean[1], 
                  tf_lower_CI = d$tf_lower_CI[1], 
                  tf_upper_CI = d$tf_upper_CI[1], 
                  cutoff = cutoff, 
                  nUT = UTmdata$nUT_final[1])
    res <- rbind(res, row)    
}

# write output
write.csv(res, snakemake@output[["combined_estimates"]], row.names = FALSE)

sink()