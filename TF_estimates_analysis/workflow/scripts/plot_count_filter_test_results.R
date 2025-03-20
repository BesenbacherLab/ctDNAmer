sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)

res_dir <- snakemake@params[["res_dir"]]

# input data
f <- snakemake@input[["cft_samples"]]
d <- read.table(f, header = TRUE)
d <- d |> 
    mutate(label_n = seq(1, nrow(d), 1), pt_label = paste0("cft_pt", label_n)) |> 
    select(-label_n, -tumor, -germline)

f <- snakemake@input[["cft_units"]]
d_cfDNA = read.table(f, header = TRUE)
d_cfDNA <- left_join(d_cfDNA, d, by = "sample_ID")

d_cfDNA <- d_cfDNA |> group_by(sample_ID) |> 
    filter(timepoint %in% c(max(timepoint), min(timepoint))) |> 
    rowwise() |>
    mutate(sample_type = ifelse(timepoint <= 0, "ctDNA-positive", "ctDNA-negative")) |>
    ungroup() |> 
    select(-cfDNA)


# aggregate count filter test results from each patient and sample
res <- NULL
for (row in 1:nrow(d_cfDNA)){
    pt <- d_cfDNA$sample_ID[row]
    fd <- d_cfDNA$cfDNA_ID[row]
    sample_type <- d_cfDNA$sample_type[row]
    pt_label <- d_cfDNA$pt_label[row]

    if (sample_type == "ctDNA-positive"){
        f <- paste0(res_dir, "results/patients/", pt, "/", fd, "/unique_tumor_count_filter_test/tf_estimation/preop_combined_tf_estimates.csv")
    } else {
        f <- paste0(res_dir, "results/patients/", pt, "/", fd, "/unique_tumor_count_filter_test/tf_estimation/postop_combined_tf_estimates.csv")
    }
    
    res_fd <- read.csv(f, header = TRUE)
    min_ci <- min(res_fd$cutoff)
    res_fd <- res_fd |> rowwise() |> mutate(cutoff_type = ifelse(cutoff == min_ci, "baseline", "count-filtered")) |> ungroup()
    res_fd <- res_fd |> mutate(pt = pt, sample_type = sample_type, 
                               pt_label = paste0(pt_label, "    (BL c_t = ", min_ci, ")"), 
                               sample_type_pt = paste(pt, sample_type, "_"))
    res <- rbind(res, res_fd)
    
}

# plot
colors <- c("ctDNA-positive" = "indianred", "ctDNA-negative" = "steelblue")
shapes <- c("baseline" = 8, "count-filtered" = 20)

options(repr.plot.width=20, repr.plot.height=8)
p <- ggplot(res) + 
    geom_point(aes(x = nUT, y = tf_mean, color = sample_type, shape = cutoff_type), size = 3) + 
    geom_line(aes(x = nUT, y = tf_mean, color = sample_type, group = sample_type_pt)) + 
    geom_ribbon(aes(x = nUT, y = tf_mean, ymin = tf_lower_CI, ymax = tf_upper_CI, fill = sample_type, group = sample_type_pt), alpha=0.5) +
    theme_bw() + 
    geom_vline(aes(xintercept = 20000), linetype = 2) + 
    facet_wrap(~pt_label, nrow = 3, labeller=label_value) + 
    scale_x_continuous(breaks = seq(0, 300000, 50000)) + 
    scale_shape_manual(name = "UT set", values = shapes) + 
    scale_color_manual(name = "cfDNA sample", values = colors) + 
    scale_fill_manual(name = "cfDNA sample", values = colors) + 
    theme(text = element_text(size = 16)) + 
    xlab("Number of UT k-mers") + 
    ylab("TF")

sink()

png(filename = snakemake@output[["cft_res"]], width = 1200, height = 500)
print(p)
dev.off()