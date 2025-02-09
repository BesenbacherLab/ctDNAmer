sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)

# input data
f = snakemake@input[["data"]]

d <- read.table(f)
colnames(d) <- c("count", "n_kmers")
d$count <- as.numeric(d$count)
d <- d[1:50, ] # remove k-mers with counts below 50

t_mean_df <- d %>% filter(count > 1) # remove k-mers observed once
filtered_n <- sum(t_mean_df$n_kmers)
mean_count <- sum(t_mean_df$count*t_mean_df$n_kmers)/filtered_n # calculate mean
t_mean_df <- t_mean_df %>% mutate(var_add = (count - mean_count)**2)
var_count <- sum(t_mean_df$var_add*t_mean_df$n_kmers)/(filtered_n-1) # calculate variance
sd_count <- sqrt(var_count) # calculate standard deviation

res <- tibble(mean_count = mean_count, 
              var_count = var_count,
              sd_count = sd_count)
write.csv(res, snakemake@output[["tumor_signal_mean"]], row.names = FALSE)

sink()