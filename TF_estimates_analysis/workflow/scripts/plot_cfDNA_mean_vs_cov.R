sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)

res_dir = snakemake@params[["res_dir"]]
cohort = snakemake@params[["cohort"]]
res_subfd = snakemake@params[["res_subfd"]]

# read in data
f <- snakemake@input[["samples"]]
d = read.table(f, header = TRUE)
d <- d |> rowwise() |> 
    mutate(sample_type = ifelse(grepl("frfr", tumor), "frfr", "ffpe")) |> 
    ungroup() |> 
    select(-tumor, -germline)

f <- snakemake@input[["units"]]
d_cfDNA = read.table(f, header = TRUE)
d_cfDNA <- left_join(d_cfDNA, d, by = "sample_ID")
d_cfDNA <- d_cfDNA |> filter(sample_type == cohort) |> select(-cfDNA)

d_cfDNA <- d_cfDNA |>
  mutate(tmp_chunks = stringr::str_split(cfDNA_ID, stringr::fixed("_"),  n = 3)) |>
  mutate(cfDNA_ID_s = map_chr(tmp_chunks, 1),
         sub_value = map_chr(tmp_chunks, 2), 
         sub_val = map_chr(tmp_chunks, 3)) |>
  select(-c(tmp_chunks, sub_value, sub_val))

# aggregate k-mer mean count estimates
d_cfDNA[c("mean_count", "sd_count")] <- NA
for (i in 1:nrow(d_cfDNA)){ 
    dat <- read.table(paste0(res_dir, d_cfDNA$sample_ID[i], "/", d_cfDNA$cfDNA_ID[i], "/", res_subfd, "cfDNA_iGL_int_cfDNAc_mean.csv"), sep = ",", header = TRUE)
    d_cfDNA[i, c("mean_count", "sd_count")] = c(dat$mean[1], dat$sd[1])
}

f <- snakemake@input[["wgs_cov"]]
covs <- read.csv(f, sep = "\t")
covs <- covs |> select(SAMPLE, MEAN_COVERAGE, SD_COVERAGE, MEDIAN_COVERAGE) |> 
    rowwise() |>
    mutate(cfDNA_ID_s = str_split_fixed(SAMPLE, "_", 2)[, 1]) |> 
    ungroup() |>
    select(-SAMPLE)
d_cfDNA <- left_join(d_cfDNA, covs, by = "cfDNA_ID_s") |> select(-cfDNA_ID_s)



options(repr.plot.width=7, repr.plot.height=7)
p_mean <- ggplot(d_cfDNA |> filter(mean_count < 100)) + 
    geom_point(aes(x = mean_count, y = MEAN_COVERAGE)) + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    annotate(geom = "text", x = 25, y = 60, label = paste0("Pearson correlation: ", round(cor(d_cfDNA$mean_count, d_cfDNA$MEAN_COVERAGE), 3)), size = 6) + 
    geom_abline(intercept = 0, slope = 1) + 
    #ggtitle() + 
    xlab("cfDNA k-mers mean count") + 
    ylab("cfDNA sample mean coverage") 

p_sd <- ggplot(d_cfDNA |> filter(mean_count < 100)) + 
    geom_point(aes(x = sd_count, y = SD_COVERAGE)) + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    geom_abline(intercept = 0, slope = 1) + 
    annotate(geom = "text", x = 10, y = 23, label = paste0("Pearson correlation: ", round(cor(d_cfDNA$sd_count, d_cfDNA$SD_COVERAGE), 3)), size = 6) + 
    xlab("cfDNA k-mers standard deviation") + 
    ylab("cfDNA sample coverage standard deviation")

sink()

png(filename = snakemake@output[["cfDNA_mean_vs_cov"]])
print(p_mean)
dev.off()

png(filename = snakemake@output[["cfDNA_sd_vs_cov"]])
print(p_sd)
dev.off()