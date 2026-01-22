sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(stringr)

# input parameters
max_cutoff <- snakemake@params[["max_cutoff"]]

# input data
res <- NULL
input_file <- snakemake@input[["data"]]
gc_data <- read.table(input_file, header = TRUE, sep = "\t")
print("Head of data")
print(head(gc_data, n = 2))



for (gc in unique(gc_data$GC)){
    min_btwn_peaks <- 2
    mean_peak <- 1
    
    d = gc_data |> filter(GC == gc) |> select(count, n_kmers)
    d$n_kmers <- as.numeric(d$n_kmers)
    d$count <- as.numeric(d$count)

    # find minimum position (count with the smallest number of k-mers) between noise and signal distributions (lower cutoff)
    d_minimal <- d %>% filter(between(count, 1, max_cutoff))
    d_minimal$count <- as.numeric(d_minimal$count)
    min_btwn_peaks <- round(optimize(approxfun(d_minimal$count,d_minimal$n_kmers),interval=c(1,10))$minimum)
    print("min_btwn_peaks")
    print(min_btwn_peaks)
    # find the mode of the signal distribution (distribution with the higher mode) 
    d_minimal <- d_minimal %>% filter(count >= min_btwn_peaks)
    mean_peak <- d_minimal$count[which(d_minimal$n_kmers == max(d_minimal$n_kmers))]

    # if there are more than one count with the largest number of k-mers, choose the smaller count
    if (length(mean_peak) > 1){
        mean_peak <- mean_peak[1]
    }
    print("mean_peak")
    print(mean_peak)
    if (length(mean_peak) == 0){
        mean_peak = 1
    }
    # if difference between minimum position (count with smallest number of k-mers between the two distributions)
    # and mode of the signal distribution is smaller than 2:
    # it is likely that the current k-mer count distribution is a uniformly decreasing distribution and there is no second (higher) mode
    # in this case, set the minimum position to 2
    if (mean_peak - min_btwn_peaks < 2){
        min_btwn_peaks <- 2
    }

    # filter the count distribution based on lower and upper cutoffs set above and calculate mean and standard deviation (based on the signal distribution)
    cfDNA_mean_df <- d %>% filter(between(count, min_btwn_peaks, max_cutoff))
    filtered_n <- sum(cfDNA_mean_df$n_kmers)
    mean_count <- sum(cfDNA_mean_df$count*cfDNA_mean_df$n_kmers)/filtered_n
    cfDNA_mean_df <- cfDNA_mean_df %>% mutate(var_add = (count - mean_count)**2)
    var_count <- sum(cfDNA_mean_df$var_add*cfDNA_mean_df$n_kmers)/(filtered_n-1)
    sd_count <- sqrt(var_count)
    
    # save result
    gc_res <- tibble(gc_content = gc,
                     mean = mean_count, 
                     var = var_count, 
                     sd = sd_count, 
                     min_btwn_peaks = min_btwn_peaks)
    res <- rbind(res, gc_res)
}

# write output
write.csv(res, snakemake@output[["cfDNA_mean_gc_strat"]], row.names = FALSE)

sink()