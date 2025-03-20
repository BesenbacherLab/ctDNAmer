sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(scales)

cohort <- snakemake@params[["cohort"]]

# input data
f <- snakemake@input[["ctDNAmers_res"]]
d <- read.csv(f)
if (cohort %in% c("frfr", "ffpe")){
    d <- d |> filter(sample_type == cohort)
}

# format data to categories
d_sum <- NULL
for (pt in unique(d$sample_ID)){ #
    d_pt <- d %>% filter(sample_ID == pt) %>% arrange(timepoint)
    pre_op <- d_pt$TF[1]

    if (is.na(d_pt$adjuvant_chemo_start_days[1])){
        after_treatment <- d_pt$TF[2]
    } else {
        after_treatment <- d_pt$TF[min(which(d_pt$timepoint > d_pt$adjuvant_chemo_end_days[1]))]
    }
    
    time_to_r <- d_pt$time_to_relapse_days[1]
    if (is.na(time_to_r)){
        relapse = "no recurrence"
        pre_relapse = NA
        last_sample = d_pt$TF[which.max(d_pt$timepoint)]
    } else {
        relapse = "recurrence"
        d_pt_pre_r <- d_pt[d_pt$timepoint <= time_to_r, ]
        pre_relapse = d_pt_pre_r$TF[which.max(d_pt_pre_r$timepoint)]
        last_sample = NA
    }
    row <- tibble(pt = pt, 
                  relapse = relapse, 
                  preop_TF = pre_op, 
                  postop_TF = after_treatment, 
                  pre_relapse_TF = pre_relapse, 
                  last_sample_TF = last_sample)
    d_sum <- rbind(d_sum, row)
}
d_sum_long <- pivot_longer(d_sum, cols = c(-pt, -relapse), names_to = "sample", values_to = "TF_estimate") |> 
    rowwise() |> 
    mutate(c_label = factor(ifelse(sample == "preop_TF", "preop", 
                                   ifelse(sample == "postop_TF", "after initial \ntreatment", 
                                          ifelse(sample == "pre_relapse_TF", "before clinical \nrecurrence", "last follow-up \ncfDNA"))), 
                            levels=c("preop", "after initial \ntreatment", "before clinical \nrecurrence", "last follow-up \ncfDNA"))) |>
    ungroup() |> 
    filter(!(is.na(TF_estimate)))


# plotting
colors <- c("recurrence" = "indianred", "no recurrence" = "steelblue")

options(repr.plot.width=10, repr.plot.height=3.5)
p <- ggplot(d_sum_long) + 
    geom_boxplot(aes(x = c_label, y = TF_estimate, fill = relapse), alpha = 0.5, outlier.shape = NA) +  
    geom_line(aes(x = c_label, y = TF_estimate, group = pt, color = relapse), linetype = 2) +  
    theme_bw() + 
    scale_fill_manual(values = colors) + 
    scale_color_manual(values = colors, name = NULL) + 
    theme(text = element_text(size = 16)) + 
    labs(fill = "Recurrence \nstatus") + 
    facet_wrap(~relapse, scales = "free_x") + 
    guides(fill = guide_legend(override.aes = list(color = NA)), 
           color = "none", 
           shape = "none")  + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    ylab("log10(TF)") + 
    xlab("")

sink()

png(filename = snakemake@output[["TF_sum"]], width = 800, height = 300)
print(p)
dev.off()