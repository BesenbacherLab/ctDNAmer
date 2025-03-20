sink(snakemake@log[[1]], append=TRUE) # logging

# packages
install.packages("ggmagnify", repos = c("https://hughjonesd.r-universe.dev", 
                 "https://cloud.r-project.org"))


library(tidyverse)
library(yardstick)
library(ggmagnify)


# input data
f <- snakemake@input[["ctDNA_mers_TF"]]
d_kmers <- read.csv(f)

f <- snakemake@input[["ctDNA_mers_det_stats"]]
d_kmers_cp <- read.csv(f)
kmers_cp <- d_kmers_cp$optimal_cutpoint[1]
auc_kmers <- d_kmers_cp$AUC[1]


f <- snakemake@input[["clonal_SNV_AF"]]
d_snvs <- read.csv(f)

f <- snakemake@input[["clonal_SNV_det_stats"]]
d_snvs_cp <- read.csv(f)
snvs_cp <- d_snvs_cp$optimal_cutpoint[1]
auc_snvs <- d_snvs_cp$AUC[1]

############## combined ROC curves ##############
d_roc_kmers <- d_kmers |> filter(ctDNA_status != "unknown") |>  filter(!is.na(TF)) |> arrange(TF)
d_roc_kmers <- d_roc_kmers |> mutate(value = TF, class = as.factor(ifelse(ctDNA_status == "pos", "Class1", "Class2")))
d_roc_kmers <- data.frame(d_roc_kmers |> select(class, value))
colnames(d_roc_kmers) <- c("y", "y_score")
roc_kmers <- roc_curve(d_roc_kmers, y, y_score)


d_roc_snvs <- d_snvs |> filter(ctDNA_status != "unknown") |>  filter(!is.na(TF)) |> arrange(TF)
d_roc_snvs <- d_roc_snvs |> mutate(value = TF, class = as.factor(ifelse(ctDNA_status == "pos", "Class1", "Class2")))
d_roc_snvs <- data.frame(d_roc_snvs |> select(class, value))
colnames(d_roc_snvs) <- c("y", "y_score")
roc_snvs <- roc_curve(d_roc_snvs, y, y_score)

colors <- c("ctDNA-mers" = "steelblue", "clonal SNVs" = "indianred")

options(repr.plot.width=10, repr.plot.height=8)
p_ROC <- ggplot() + 
    geom_path(data = roc_kmers, aes(x = 1-specificity, y = sensitivity, color = "ctDNA-mers"), linewidth = 1) + 
    geom_path(data = roc_snvs, aes(x = 1-specificity, y = sensitivity, color = "clonal SNVs"), linewidth = 1) + 
    theme_bw() + 
    theme(text = element_text(size = 25))  + 
    geom_abline(intercept = 0, slope = 1) + 
    annotate(geom = "point", 
             x = 1-roc_kmers[which(roc_kmers$.threshold == kmers_cp), "specificity"][[1]], 
             y = roc_kmers[which(roc_kmers$.threshold == kmers_cp), "sensitivity"][[1]], size = 4) + 
     annotate(geom = "segment", 
             x = 1-roc_kmers[which(roc_kmers$.threshold == kmers_cp), "specificity"][[1]], 
             y = roc_kmers[which(roc_kmers$.threshold == kmers_cp), "sensitivity"][[1]], 
             xend = 1-roc_kmers[which(roc_kmers$.threshold == kmers_cp), "specificity"][[1]], 
             yend = 1-roc_kmers[which(roc_kmers$.threshold == kmers_cp), "specificity"][[1]], color = "steelblue", linetype = 2, linewidth = 1.3) + 

    annotate(geom = "point", 
             x = 1-roc_snvs[which(roc_snvs$.threshold == snvs_cp), "specificity"][[1]], 
             y = roc_snvs[which(roc_snvs$.threshold == snvs_cp), "sensitivity"][[1]], size = 4) + 
     annotate(geom = "segment", 
             x = 1-roc_snvs[which(roc_snvs$.threshold == snvs_cp), "specificity"][[1]], 
             y = roc_snvs[which(roc_snvs$.threshold == snvs_cp), "sensitivity"][[1]], 
             xend = 1-roc_snvs[which(roc_snvs$.threshold == snvs_cp), "specificity"][[1]], 
             yend = 1-roc_snvs[which(roc_snvs$.threshold == snvs_cp), "specificity"][[1]], color = "indianred", linetype = 2, linewidth = 1.3) + 
    annotate(geom = "text", 
             x = 0.75, 
             y = 0.3, 
             label = paste0("clonal SNVs AUC: ", round(auc_snvs, 3)), 
             size = 7, color = "indianred") + 
    annotate(geom = "text", 
             x = 0.75, 
             y = 0.35, 
             label = paste0("ctDNA-mers AUC: ", round(auc_kmers, 3)), 
             size = 7, color = "steelblue") + 
    scale_color_manual(name = "Method", values = colors)


############## Correlation plot ##############

# prepare data
d_kmers_sub <- d_kmers |> select(-sample_type, -start_of_last_intervention, -cfDNA_mean, -wt, -tphi, -wgl, -glphi, -wn)|> 
    rename(ctDNA_detected_kmers = ctDNA_detected, 
           TF_kmers = TF, 
           TF_lCI_kmers = TF_lCI,
           TF_uCI_kmers = TF_uCI)

d_snvs_sub <- d_snvs |> select(-sample_type, -start_of_last_intervention, -n_pass, -n_obs) |> 
    rename(ctDNA_detected_snvs = ctDNA_detected,
           TF_snvs = TF, 
           TF_lCI_snvs = TF_lCI,
           TF_uCI_snvs = TF_uCI)

d_ff <- left_join(d_kmers_sub, d_snvs_sub, by = c("sample_ID", "cfDNA_ID", "timepoint", "relapse", "time_to_relapse_days", "adjuvant_chemo_start_days", "adjuvant_chemo_end_days", "ctDNA_status"))
d_ff_sub <- d_ff |> 
    filter(ctDNA_status %in% c("pos", "neg")) |> 
    rowwise() |> 
    mutate(ctDNA_status = ifelse(ctDNA_status == "pos", "positive", "negative")) |>
    ungroup()

cor_TF_AF <- cor(d_ff_sub$TF_snvs, d_ff_sub$TF_kmers, use = "complete.obs")


# plot
colors <- c("positive" = 'indianred', "negative" = 'steelblue')
p <- ggplot(d_ff_sub) + 
    geom_point(aes(x = TF_kmers, y = TF_snvs, color = ctDNA_status), size = 3) + 
    geom_abline(intercept = 0, slope = 1, color = "darkgrey", linetype = 2, linewidth = 1) + 
    theme_bw() + 
    theme(text = element_text(size = 25)) +
    scale_color_manual(name = "ctDNA", values=colors) +
    xlab("ctDNA-mers TF") + 
    ylab("clonal SNVs mean AF") + 
    annotate(geom = "text", 
             x = 0.15, #0.015, 
             y = 0.38, #0.095, 
             label = paste0("Pearson cor.: ", round(cor_TF_AF, 3)), size = 7) +
    aes(xmax=1) + 
    aes(ymax=0.4)

# add magnification
from <- list(0, 0.02, 0, 0.02)
to <- list(0.7, 1, 0.02, 0.22)
p_cor <- p + geom_magnify(from = from, to = to)

sink()

png(filename = snakemake@output[["ROC_curve"]], width = 700, height = 500)
print(p_ROC)
dev.off()

png(filename = snakemake@output[["correlation_plot"]], width = 700, height = 500)
print(p_cor)
dev.off()
