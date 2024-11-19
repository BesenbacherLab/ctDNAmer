sink(snakemake@log[[1]], append=TRUE)

library(tidyverse)
library(stringr)


max_count <- snakemake@params[["max_count"]]
max_cutoff <- snakemake@params[["max_cutoff"]]
    
res <- NULL
input_file <- snakemake@input[["data"]]
gc_data <- read.table(input_file, header = TRUE, sep = "\t")
print("Head of data")
print(head(gc_data, n = 2))

for (gc in unique(gc_data$GC)){
    
    d = gc_data |> filter(GC == gc) |> select(count, n_kmers)
    d$n_kmers <- as.numeric(d$n_kmers)
    d$count <- as.numeric(d$count)

    # find min between noise and signal peak (lower cutoff), find the peak of signal
    d_minimal <- d %>% filter(between(count, 1, max_count))
    d_minimal$count <- as.numeric(d_minimal$count)
    min_btwn_peaks <- round(optimize(approxfun(d_minimal$count,d_minimal$n_kmers),interval=c(1,10))$minimum)
    d_minimal <- d_minimal %>% filter(count >= min_btwn_peaks)
    mean_peak <- d_minimal$count[which(d_minimal$n_kmers == max(d_minimal$n_kmers))]

    if (length(mean_peak) > 1){
        mean_peak <- mean_peak[1]
    }
    
    # if difference between minimum between peaks and mean peak is smaller than 2, it is likely that we have a uniformly decreasing distribution and no second peak, set min btwn peaks to 2 and calculate peak as the mean count
    if (mean_peak - min_btwn_peaks < 2){
        cfDNA_mean_df <- d %>% filter(between(count, 1, max_cutoff))
        filtered_n <- sum(cfDNA_mean_df$n_kmers)
        mean_peak <- round(sum(cfDNA_mean_df$count*cfDNA_mean_df$n_kmers)/filtered_n)
        min_btwn_peaks <- 2
    }

    # filter based on lower and upper cutoffs set above and calculate mean and standard deviation
    cfDNA_mean_df <- d %>% filter(between(count, min_btwn_peaks, max_cutoff))
    filtered_n <- sum(cfDNA_mean_df$n_kmers)
    mean_count <- sum(cfDNA_mean_df$count*cfDNA_mean_df$n_kmers)/filtered_n
    cfDNA_mean_df <- cfDNA_mean_df %>% mutate(var_add = (count - mean_count)**2)
    var_count <- sum(cfDNA_mean_df$var_add*cfDNA_mean_df$n_kmers)/(filtered_n-1)
    sd_count <- sqrt(var_count)
        
    gc_res <- tibble(gc_content = gc,
                     mean = mean_count, 
                     var = var_count, 
                     sd = sd_count, 
                     min_btwn_peaks = min_btwn_peaks)
    res <- rbind(res, gc_res)
}


write.csv(res, snakemake@output[["cfDNA_mean_gc_strat"]], row.names = FALSE)

sink()