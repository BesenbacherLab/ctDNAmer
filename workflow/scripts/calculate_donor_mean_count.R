sink(snakemake@log[[1]], append=TRUE)

library(tidyverse)


max_count = snakemake@params[["max_count"]]


f <- snakemake@input[["count_sum_filt"]]
d <- read.table(f)
colnames(d) <- c("count", "n_kmers")
d$count <- as.numeric(d$count)

# find min between noise and signal peak (lower cutoff), find the peak of signal
d_minimal <- d %>% filter(between(count, 1, max_count))
min_btwn_peaks <- round(optimize(approxfun(d_minimal$count,d_minimal$n_kmers),interval=c(1,10))$minimum)
d_minimal <- d_minimal %>% filter(count >= min_btwn_peaks)
mean_peak <- d_minimal$count[which(d_minimal$n_kmers == max(d_minimal$n_kmers))]

# find upper cutoff (min count which holds less than 0.5% of k-mers)
d_low_rm <- d %>% filter(count >= min_btwn_peaks)
d_sum_n <- sum(d_low_rm$n_kmers)
d_low_rm <- d_low_rm %>% filter(count >= mean_peak)
max_cutoff <- min(d_low_rm$count[which(d_low_rm$n_kmers <= ((d_sum_n*0.5)/100))])

# if difference between minimum between peaks and mean peak is smaller than 2, it is likely that we have a uniformly decreasing distribution and no second peak, set min btwn peaks to 2 and calculate peak as the mean count
if (mean_peak - min_btwn_peaks < 2){
    print(paste0("Special case: distribution seems to be uniformly decreasing, setting min_btwn_peaks to 2"))
    min_btwn_peaks <- 2
}

# filter based on lower and upper cutoffs set above and calculate mean and standard deviation
mean_df <- d %>% filter(between(count, min_btwn_peaks, max_cutoff))
filtered_n <- sum(mean_df$n_kmers)
mean_count <- sum(mean_df$count*mean_df$n_kmers)/filtered_n
mean_df <- mean_df %>% mutate(var_add = (count - mean_count)**2)
var_count <- sum(mean_df$var_add*mean_df$n_kmers)/(filtered_n-1)
sd_count <- sqrt(var_count)
    

res <- tibble(mean = mean_count, var = var_count, sd = sd_count, min_btwn_peaks = min_btwn_peaks, max_cutoff = max_cutoff)
write.csv(res, snakemake@output[["mean_count"]], row.names = FALSE)

sink()