sink(snakemake@log[[1]], append=TRUE)  # logging

# packages
library(tidyverse)
library(stringr)

# input data and parameters
tracked_SNVs_fnames <- snakemake@input[["SNVs"]] # input tracked SNVs data

pt_id <- snakemake@params[["pt_id"]]
cfDNA_timepoints <- snakemake@params[["cfDNA_timepoints"]]

# calculate mean AF
res_m <- NULL
res_dist <- NULL

for (i in 1:length(tracked_SNVs_fnames)){
    fname = tracked_SNVs_fnames[i] # file name (SNV cfDNA information)
    print(paste0("Working on file: ", fname))
    
    # find cfDNA id from file path
    path_split = str_split(fname, "/")[[1]]
    for (j in 1:length(path_split)){
        sub_path = path_split[j]
        if (sub_path == pt_id){
            cfDNA_id = path_split[j+1]
        }
    }
    # get cfDNA smple timepoint
    timepoint = as.numeric(cfDNA_timepoints[cfDNA_id])

    AF_table <- read.csv(fname, sep = "\t") 
    n_pass <- nrow(AF_table)
    n_obs <- sum(AF_table$cfDNA_n_alt_total > 0)

    # calculate allele frequency
    AF_table <- AF_table %>% 
        select(chrom, start, end, ref, alt, cfDNA_n_alt_total, cfDNA_DP_total) %>% 
        mutate(cfDNA_id = cfDNA_id, 
               timepoint = as.numeric(timepoint),
               AF = cfDNA_n_alt_total/cfDNA_DP_total, 
               label = paste0("Timepoint: ", timepoint, 
               "\nn mean AF: ", round(mean(AF), 4), 
               "\nn pass var: ", n_pass,
               "\nn obs var: ", n_obs))
    
    # calculate mean allele frequency
    mean_AF <- mean(AF_table$AF)
    
    mean_AF_df <- tibble(pt_id = pt_id, cfDNA_id = cfDNA_id, timepoint = as.numeric(timepoint),
                        mean_AF = mean_AF, n_pass = n_pass, n_obs = n_obs, i = i)
    
    # bind to results DF's
    res_dist = rbind(res_dist, AF_table)
    res_m <- rbind(res_m, mean_AF_df)
}

# plot AF distribution
p_dist <- ggplot(res_dist %>% filter(AF > 0)) + 
    geom_histogram(aes(x = AF), binwidth = 0.005, fill = "grey", color = "darkgrey") + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    facet_wrap(~label, nrow = 1) + 
    ylab("Count")  + 
    xlab("Allele frequency")


# plot mean AF
colors1 <- c("Surgery/First treatment" = "firebrick4")
p_meanAF <- ggplot(data = res_m) + 
    # Clinical data
    geom_vline(data = NULL, aes(xintercept=0, color = "Surgery/First treatment"), size = 1.1) +          # Surgery
    scale_color_manual(name="Clinical timepoints", values=colors1) +

    # TF estimates - clonal mutations 
    geom_line(data = res_m, aes(x = timepoint, y = mean_AF), size = 0.7) + 
    geom_point(data = res_m, aes(x = timepoint, y = mean_AF), size = 3) +
    theme_bw() +
    theme(text = element_text(size = 16),legend.box.margin = margin(10, 0, 0, 0)) + #
    theme(legend.position = "bottom",  legend.box = "vertical", legend.spacing.y = unit(0.05, 'cm')) + # 
    xlab("Sample time point (days since first treatment)") +
    ylab("mean AF")


# write results
ggsave(snakemake@output[["mean_VAF_dist"]], p_dist, create.dir = TRUE, device='png', width = 20, height = 3, dpi=100)
ggsave(snakemake@output[["mean_VAF_timeline"]], p_meanAF, create.dir = TRUE,  device='png', width = 10, height = 5, dpi=150)
write.table(res_m, snakemake@output[["mean_VAF"]], row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)

sink()