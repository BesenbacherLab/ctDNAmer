sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(cutpointr)
library(yardstick)


# input data
f <- snakemake@input[["units"]]
d <- read.csv(f)

# params
colors <- c("pos" = "indianred", "neg" = "steelblue", "pos., TF > 0.001" = "indianred1", "neg., TF > 0.001" = "steelblue1")
pl_cutoff <- 0.001
label_pos = 'pos., TF > 0.001'
label_neg = 'neg., TF > 0.001'



d_obs <- d |> filter(ctDNA_status != "unknown") |>  filter(!is.na(TF)) |> arrange(TF)
# maximize Youden's J
opt_cut <- cutpointr(data = d_obs, 
                     x = TF, 
                     class = ctDNA_status, 
                     direction = ">=", # larger value for positive class
                     pos_class = "pos",
                     neg_class = "neg", 
                     method = maximize_metric, 
                     metric = youden)

sum_opt_cut <- summary(opt_cut)
sum_cutpointr <- sum_opt_cut$cutpointr
cutpoint = sum_cutpointr[[1]][c("optimal_cutpoint")][[1]]
print(paste0("Number of observations: ", sum_opt_cut$n_obs, ", n pos: ", sum_opt_cut$n_pos, ", n neg: ", sum_opt_cut$n_neg))
print(paste0("Cutpoint chosen based on maximizing Youden's J: ", round(cutpoint, 8)))

############# ROC curve #############
d_roc <- d_obs
d_roc <- d_roc |> mutate(value = TF, class = as.factor(ifelse(ctDNA_status == "pos", "Class1", "Class2")))
d_roc <- data.frame(d_roc |> select(class, value))
colnames(d_roc) <- c("y", "y_score")
roc <- roc_curve(d_roc, y, y_score)

p1 <- ggplot(roc) +
    geom_path(aes(x = 1-specificity, y = sensitivity), linewidth = 0.5) + 
    theme_bw() + 
    theme(text = element_text(size = 20))  + 
    geom_abline(intercept = 0, slope = 1) + 
    annotate(geom = "point", 
             x = 1-roc[which(roc$.threshold == cutpoint), "specificity"][[1]], 
             y = roc[which(roc$.threshold == cutpoint), "sensitivity"][[1]], size = 4) + 
     annotate(geom = "segment", 
             x = 1-roc[which(roc$.threshold == cutpoint), "specificity"][[1]], 
             y = roc[which(roc$.threshold == cutpoint), "sensitivity"][[1]], 
             xend = 1-roc[which(roc$.threshold == cutpoint), "specificity"][[1]], 
             yend = 1-roc[which(roc$.threshold == cutpoint), "specificity"][[1]], color = "indianred", linetype = 2, linewidth = 1.3) + 
    annotate(geom = "text", 
             x = 0.75, 
             y = 0.3, 
             label = paste0("AUC: ", round(sum_cutpointr[[1]][c("AUC")]$AUC, 3), 
                            "\nmax. Youden: ", round(sum_cutpointr[[1]][c("youden")]$youden, 2)), 
             size = 6, color = "black")


############# plotting the TF values #############
d_obs_bin <- d_obs |> 
    rowwise() |> 
    mutate(TF_bin = ifelse(TF > pl_cutoff, pl_cutoff, TF), 
           ctDNA_status_l = factor(ctDNA_status, 
                                   levels = c("pos", "neg"), 
                                   labels = c("ctDNA positive (preop & recur. pt postop until last treatment)", "ctDNA negative (non-recur. pt postop after treatment)")),
           ctDNA_status_bin = ifelse(TF > pl_cutoff, ifelse(ctDNA_status == "pos", label_pos, label_neg), ctDNA_status)) |> 
    ungroup() |> 
    arrange(TF)
d_obs_bin$cfDNA_ID_f <- factor(d_obs_bin$cfDNA_ID, levels=d_obs_bin$cfDNA_ID)  

p2 <- ggplot(d_obs_bin) +  
    geom_col(aes(x = ctDNA_status, 
                 y = TF_bin, 
                 group = cfDNA_ID_f, 
                 fill = ctDNA_status_bin),
                 size=0.2,
                 color = "lightgrey", position = position_dodge(width = 1.1)) +
    geom_hline(aes(yintercept = cutpoint), linetype = 2, linewidth = 0.7, color = "black") + 
    scale_fill_manual(name = "Sample ctDNA status", values = colors) + 
    theme_bw() + 
    theme(text = element_text(size = 20), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          legend.position = "bottom", 
          legend.spacing.y = unit(2.0, 'cm'),
          legend.text = element_text(margin = margin(t = 5, l = 5, b = 5))) + 
    xlab("Sample") + 
    ylab("TF") + 
    facet_wrap(~ctDNA_status_l, scales = "free_x", labeller=label_value, nrow = 1) + 
  guides(fill = guide_legend(byrow = TRUE)) 


d_f_nozeros <- d |> rowwise() |> mutate(ctDNA_detected = ifelse(TF >= cutpoint, "pos", "neg")) |> ungroup()
d_f <- d |> rowwise() |> mutate(ctDNA_detected = ifelse(TF >= cutpoint, "pos", "neg"), 
                                        TF = ifelse(ctDNA_detected == "neg", 0, TF), 
                                        TF_uCI = ifelse(ctDNA_detected == "neg", 0, TF_uCI), 
                                        TF_lCI = ifelse(ctDNA_detected == "neg", 0, TF_lCI)) |> ungroup()

det_stats <- sum_cutpointr[[1]][c("optimal_cutpoint", "youden", "acc", "sensitivity", "specificity", "AUC")] |> as.data.frame()
conf_m <- sum_opt_cut$confusion_matrix[[1]]

write.csv(d_f, snakemake@output[["results"]], row.names = FALSE)
write.csv(d_f_nozeros, snakemake@output[["results_noTF_zero"]], row.names = FALSE)
write.csv(det_stats, snakemake@output[["detection_stats"]], row.names = FALSE)
write.csv(conf_m, snakemake@output[["conf_matrix"]], row.names = FALSE)

sink()

png(filename = snakemake@output[["TF_ROC"]])
print(p1)
dev.off()

png(filename = snakemake@output[["TF_barplot"]], width = 1200, height = 400)
print(p2)
dev.off()