sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)

# input parameters
max_count = snakemake@params[["max_count"]]

# input data
f <- snakemake@input[["data"]]
d <- read.table(f)
colnames(d) <- c("count", "n_kmers")
d$count <- as.numeric(d$count)

# find minimum position (count with the smallest number of k-mers) between noise and signal distributions (lower cutoff)
d_minimal <- d %>% filter(between(count, 1, max_count))
min_btwn_peaks <- round(optimize(approxfun(d_minimal$count,d_minimal$n_kmers),interval=c(1,10))$minimum)
# find the mode of the signal distribution (distribution with the higher mode) 
d_minimal <- d_minimal %>% filter(count >= min_btwn_peaks)
mean_peak <- d_minimal$count[which(d_minimal$n_kmers == max(d_minimal$n_kmers))]

# find the upper cutoff (min count which holds less than 0.5% of the remaining k-mers (after removing low-count k-mers)
d_low_rm <- d %>% filter(count >= min_btwn_peaks)
d_sum_n <- sum(d_low_rm$n_kmers)
d_low_rm <- d_low_rm %>% filter(count >= mean_peak)
max_cutoff <- min(d_low_rm$count[which(d_low_rm$n_kmers <= ((d_sum_n*0.5)/100))])

# if difference between minimum position (count with smallest number of k-mers between the two distributions)
# and mode of the signal distribution is smaller than 2:
# it is likely that the current k-mer count distribution is a uniformly decreasing distribution and there is no second (higher) mode
# in this case, set the minimum position to 2
if (mean_peak - min_btwn_peaks < 2){
    print(paste0("Special case: distribution seems to be uniformly decreasing, setting min_btwn_peaks to 2"))
    min_btwn_peaks <- 2
}

# filter the count distribution based on lower and upper cutoffs set above and calculate mean and standard deviation (based on the signal distribution)
mean_df <- d %>% filter(between(count, min_btwn_peaks, max_cutoff))
filtered_n <- sum(mean_df$n_kmers)
mean_count <- sum(mean_df$count*mean_df$n_kmers)/filtered_n
mean_df <- mean_df %>% mutate(var_add = (count - mean_count)**2)
var_count <- sum(mean_df$var_add*mean_df$n_kmers)/(filtered_n-1)
sd_count <- sqrt(var_count)
    
# write output
res <- tibble(mean = mean_count, 
              var = var_count, 
              sd = sd_count, 
              min_btwn_peaks = min_btwn_peaks, 
              max_cutoff = max_cutoff)
write.csv(res, snakemake@output[["mean_count"]], row.names = FALSE)

sink()