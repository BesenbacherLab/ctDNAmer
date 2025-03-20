sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)


cohort <- snakemake@params[["cohort"]]

# input data
f <- snakemake@input[["ctDNAmers_res"]]
d <- read.csv(f)
if (cohort %in% c("frfr", "ffpe")){
    d <- d |> filter(sample_type == cohort)
}

# calculate and plot lead time
d_lead <- d |>
    filter(relapse == 1, timepoint > 0, ctDNA_detected == "pos") |>
    select(sample_ID, timepoint, time_to_relapse_days, ctDNA_detected) |>
    group_by(sample_ID) |>
    summarize(kmers_relapse = min(timepoint), time_to_relapse_days = time_to_relapse_days[1]) |>
    mutate(lead_time_days = time_to_relapse_days - kmers_relapse, 
           lead_time_months = lead_time_days/365*12, 
           kmers_relapse = kmers_relapse,
           kmers_relapse_months = kmers_relapse/365*12,
           time_to_relapse_months = time_to_relapse_days/365*12) |>
    arrange(lead_time_months)  

options(repr.plot.width=8, repr.plot.height=5)
colors <- c("ctDNA-mers" = "darkgreen", "imaging" = "firebrick")

d_lead$sample_ID <- reorder(d_lead$sample_ID, d_lead$lead_time_months)
median_lead <- median(d_lead$lead_time_months)

# plotting
#options(repr.plot.width=6, repr.plot.height=2.5) 
options(repr.plot.width=7, repr.plot.height=3.5)
p <- ggplot(d_lead, aes(x=sample_ID, y = lead_time_months)) + 
    geom_segment(aes(x=as.numeric(sample_ID)+.001, xend=as.numeric(sample_ID)+.001,
                     y = kmers_relapse_months, yend = time_to_relapse_months), size = 0.7) + 
    geom_point(aes(x=as.numeric(sample_ID)+.001, y = kmers_relapse_months, color = "ctDNA-mers"), size = 3) +
    geom_point(aes(x=as.numeric(sample_ID)+.001, y = time_to_relapse_months, color = "imaging"), size = 3) +
    theme_bw() + 
    xlab("Patient") + 
    ylab("Time until recurrence (months)") + 
    theme(axis.text.y = element_blank(), 
          text = element_text(size = 16), 
          legend.spacing.y = unit(0.1, 'cm'), 
          legend.spacing.x = unit(0.1, 'cm'), 
          legend.position = "right") + 
    labs(color = "Recurrence \ndetection") +
    scale_color_manual(values = colors) + 
    coord_flip() + 
    annotate(geom = "text", y = 25, x = as.numeric(d_lead$sample_ID[4])+.001, label = paste0("Median lead time: \n", round(median_lead, 1), " months"), size = 5) 


sink()

png(filename = snakemake@output[["lead_time"]], width = 480, height = 240)
print(p)
dev.off()