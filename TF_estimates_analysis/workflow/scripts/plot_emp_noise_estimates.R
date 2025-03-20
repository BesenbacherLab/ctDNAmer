sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)

res_dir <- snakemake@params[["res_dir"]]
cohort <- snakemake@params[["cohort"]]

# input data
f <- snakemake@input[["samples"]]
d <- read.table(f, header = TRUE)
d <- d |> rowwise() |> 
    mutate(sample_type = ifelse(grepl("frfr", tumor), "frfr", "ffpe")) |> 
    ungroup() |> 
    select(-tumor, -germline)

# aggregate empirical noise distribution estimates
d[c("noise_mean_m30", "noise_var_m30")] <- NA
for (row in 1:nrow(d)){
    pt <- d$sample_ID[row]

    n_est <- read.table(paste0(res_dir, "results/patients/", pt, "/empirical_noise/estimates.csv"), sep = ",", header = TRUE)
    mean <- n_est$mean_mu[1]
    phi <- n_est$mean_phi[1]

    d[row, "noise_mean_m30"] = mean*30 # assuming a constant kmer mean count of 30
    d[row, "noise_var_m30"] = (mean*30 + ((mean*30)**2/phi))
}


#plotting
if (cohort %in% c("frfr", "ffpe")){
    options(repr.plot.width=7, repr.plot.height=6)
    p <- ggplot(d |> filter(sample_type == cohort)) + 
        geom_point(aes(x = noise_mean_m30, y = noise_var_m30)) + 
        theme_bw() + 
        theme(text = element_text(size = 18), 
              plot.margin = margin(t = 0.1, r = 0.5, b = 0.1, l = 0.5, unit = "cm")) + 
        geom_abline(linetype = 2) + 
        xlab("Mean") + 
        ylab("Variance") 
    
} else {
    options(repr.plot.width=5, repr.plot.height=5)
    colors = c("frfr" = "darkseagreen", "ffpe" = "indianred")
    p <- ggplot(d) + 
        geom_point(aes(x = noise_mean_m30, y = noise_var_m30, color = sample_type)) + 
        theme_bw() + 
        theme(text = element_text(size = 18), 
              legend.position = "bottom") + 
        scale_color_manual(name = "Sample type", values = colors) + 
        geom_abline(linetype = 2) + 
        xlab("Mean") + 
        ylab("Variance") 
        
}


sink()

png(filename = snakemake@output[["emp_noise_estimates"]])
print(p)
dev.off()