sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(ggnewscale)


# input data
f <- snakemake@input[["ctDNA_mers_TF"]]
d_kmers <- read.csv(f)

f <- snakemake@input[["clonal_SNV_AF"]]
d_snvs <- read.csv(f)

# Prepare data
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


# set colors
colors1 <- c("surgery" = "black", "relapse" = "steelblue4")
colors2 <- c("No" = 'steelblue4', "Yes" = 'firebrick4')
colors3 <- c("ctDNA-mers TF" = "black", "clonal SNVs mean AF" = "darkgrey") #, "clonal SNVs median" = "lightgrey"
colors_imaging <- c("No abnormality" = "palegreen4", 'Suspect_finding' = "lightpink1",
                    "Recurrence" = "firebrick1",
                    "Regression" = "darkolivegreen3",
                    "Progression" = "firebrick",
                    "Other cancer" = "maroon",
                    "Stable disease" = "khaki1")
colors_intervention <- c("adjuvant \nchemotherapy" = "cadetblue4", 
                         "intervention \nchemo- or radiation therapy" = "coral1", 
                         "intervention \n(surgery, RFA, MWA, radiotherapy)" = "coral1")

# get imaging info
f <- snakemake@input[["imagings"]]
imagings <- read.csv(f)
imagings <- imagings |> 
    select(patient_id, imaging_result, time_since_op_days) |> 
    mutate(imaging_res_label = ifelse(imaging_result == 0, "No abnormality",
                                      ifelse(imaging_result == 1, "Recurrence", 
                                             ifelse(imaging_result == 2, "Stable disease", 
                                                   ifelse(imaging_result == 3, "Regression", 
                                                         ifelse(imaging_result == 4, "Progression", 
                                                               ifelse(imaging_result == 5, "Suspect_finding", "Other cancer"))))))) |>
    rename(sample_ID = patient_id)


# get recurrence treatment info
f <- snakemake@input[["interventions"]]
interventions = read.csv(f)
interventions <- interventions |>
    select(patient_id, intervention_start_days, intervention_end_days, intervention_type) |> 
    rename(sample_ID = patient_id)

intervention_chemo = interventions |> filter(intervention_type %in% c(4, 5, 12)) |> select(sample_ID, intervention_start_days, intervention_end_days)
intervention_other = interventions |> filter(!intervention_type %in% c(4, 5, 12)) |> select(sample_ID, intervention_start_days, intervention_end_days)



########### Recurring patient timelines ###########

# data prep
data_ff <- d_ff |> filter(relapse == 1) |> arrange(sample_ID)
df_pt_nrs = tibble(sample_ID = unique(data_ff$sample_ID), pt_nr = seq(1, length(unique(data_ff$sample_ID)), 1))
data_ff <- data_ff |> mutate(ctDNA_detected_kmers = ifelse(ctDNA_detected_kmers == "neg", "No", "Yes"), 
                             ctDNA_detected_snvs = ifelse(ctDNA_detected_snvs == "neg", "No", "Yes"))


data_ff <- left_join(data_ff, df_pt_nrs, by = "sample_ID")
data_ff <- data_ff |> rowwise() |> mutate(pt_label = paste0("pt_rec_", pt_nr)) |> ungroup()
data_ff_label <- data_ff |> select(sample_ID, pt_label) |> distinct()
d_minmax <- data_ff |> group_by(sample_ID) |> mutate(TF_min = min(TF_lCI_kmers, TF_snvs), TF_max = max(TF_uCI_kmers, TF_snvs)) |> 
    select(sample_ID, TF_min, TF_max) |> distinct() |> ungroup()

ad_chem_ff <- data_ff |> select(sample_ID, adjuvant_chemo_start_days, adjuvant_chemo_end_days) |> distinct()
ad_chem_ff <- left_join(ad_chem_ff, d_minmax, by = "sample_ID") 
ad_chem_ff <- left_join(ad_chem_ff, data_ff_label, by = "sample_ID")
ad_chem_ff <- ad_chem_ff |> select(-sample_ID)

imaging_ff <- imagings |> filter(sample_ID %in% unique(data_ff$sample_ID)) 
imaging_ff <- left_join(imaging_ff, data_ff_label, by = "sample_ID")
imaging_ff <- imaging_ff |> select(-sample_ID)


intervention_chemo_ff <- left_join(intervention_chemo |> filter(sample_ID %in% unique(data_ff$sample_ID)) |> distinct(), d_minmax, by = "sample_ID")
intervention_chemo_ff <- left_join(intervention_chemo_ff, data_ff_label, by = "sample_ID")
intervention_chemo_ff <- intervention_chemo_ff |> select(-sample_ID)

intervention_other_ff <- left_join(intervention_other |> filter(sample_ID %in% unique(data_ff$sample_ID)) |> distinct(), d_minmax, by = "sample_ID")
intervention_other_ff <- left_join(intervention_other_ff, data_ff_label, by = "sample_ID")
intervention_other_ff <- intervention_other_ff |> select(-sample_ID)

# plot
options(repr.plot.width=15, repr.plot.height=15)
p_recur <- ggplot(data = data_ff) + 

    # Clinical data
    geom_vline(data = NULL, aes(xintercept=0, color = "surgery"), linewidth = 1.1, linetype = 8) +      # Surgery
    scale_color_manual(name="Clinical timepoints", values=colors1) + # 
    new_scale_color() +

    # TF estimates
    geom_ribbon(aes(x = timepoint, ymin = TF_lCI_kmers, ymax = TF_uCI_kmers), alpha = 0.2, color="grey") + 
    geom_point(aes(x = timepoint, y = TF_kmers, color = ctDNA_detected_kmers), size = 2) + #
    scale_color_manual(name = "ctDNA detected", values=colors2) + 
    new_scale_color() +

    # clonal SNVs mean AF
    
    geom_point(aes(x = timepoint, y = TF_snvs, color = ctDNA_detected_snvs), size = 2) + #
    scale_color_manual(name = "ctDNA detected", values=colors2) + 
    new_scale_color() +

    geom_line(aes(x = timepoint, y = TF_kmers, color = "ctDNA-mers TF"), size = 0.7) + 
    geom_line(aes(x = timepoint, y = TF_snvs, color = "clonal SNVs mean AF"), size = 0.7) + 
    scale_color_manual(name = "TF estimation", values=colors3, limits = c("ctDNA-mers TF", "clonal SNVs mean AF")) + #, "clonal SNVs median"
    new_scale_color() +

    #imaging
    geom_vline(data = imaging_ff, aes(xintercept = time_since_op_days, color = imaging_res_label), linewidth = 1) + 
    scale_color_manual(name = "Imaging", values=colors_imaging) +
    new_scale_color() +

    

    # adjuvant chemo
    geom_rect(data = ad_chem_ff,
              aes(xmin = adjuvant_chemo_start_days, xmax = adjuvant_chemo_end_days, ymin = TF_min, ymax = TF_max,
                  fill = "adjuvant \nchemotherapy"), alpha = 0.3) + 

    # intervention administered once
      geom_rect(data = intervention_other_ff,                                             # Intervention
               aes(xmin = intervention_start_days-10, xmax = intervention_start_days+10, ymin = TF_min, ymax = TF_max,
                   fill = "intervention \n(surgery, RFA, MWA, radiotherapy)"), alpha = 0.3) +
    
    # intervention chemo
    geom_rect(data = intervention_chemo_ff,
             aes(xmin = intervention_start_days, xmax = intervention_end_days, ymin = TF_min, ymax = TF_max,
                 fill = "intervention \nchemo- or radiation therapy"), alpha = 0.3) + 
    scale_fill_manual(name = "Treatment", values = colors_intervention) + 
    
    theme_bw() +
    theme(text = element_text(size = 11), 
          legend.box.margin = margin(6, 6, 6, 6), #116, 6, 6, 6
          legend.position = "right", 
          legend.box = "vertical",
          legend.spacing.y = unit(-0.1, 'cm')) + 
    xlab("Sample time point (days since surgery)") +
    facet_wrap(~pt_label, ncol = 2, scales = "free_y") + 
    ylab("")


########### Non-recurring patient timelines ###########

# data prep
data_ff <- d_ff |> filter(relapse == 0) |> arrange(sample_ID)
df_pt_nrs = tibble(sample_ID = unique(data_ff$sample_ID), pt_nr = seq(1, length(unique(data_ff$sample_ID)), 1))
data_ff <- data_ff |> mutate(ctDNA_detected_kmers = ifelse(ctDNA_detected_kmers == "neg", "No", "Yes"), 
                             ctDNA_detected_snvs = ifelse(ctDNA_detected_snvs == "neg", "No", "Yes"))

data_ff <- left_join(data_ff, df_pt_nrs, by = "sample_ID")
data_ff <- data_ff |> rowwise() |> mutate(pt_label = paste0("pt_no_rec_", pt_nr)) |> ungroup()
data_ff_label <- data_ff |> select(sample_ID, pt_label) |> distinct()
d_minmax <- data_ff |> group_by(sample_ID) |> mutate(TF_min = min(TF_lCI_kmers, TF_snvs), TF_max = max(TF_uCI_kmers, TF_snvs)) |> 
    select(sample_ID, TF_min, TF_max) |> distinct() |> ungroup()

ad_chem_ff <- data_ff |> select(sample_ID, adjuvant_chemo_start_days, adjuvant_chemo_end_days) |> distinct()
ad_chem_ff <- left_join(ad_chem_ff, d_minmax, by = "sample_ID") 
ad_chem_ff <- left_join(ad_chem_ff, data_ff_label, by = "sample_ID")
ad_chem_ff <- ad_chem_ff |> select(-sample_ID)

imaging_ff <- imagings |> filter(sample_ID %in% unique(data_ff$sample_ID)) 
imaging_ff <- left_join(imaging_ff, data_ff_label, by = "sample_ID")
imaging_ff <- imaging_ff |> select(-sample_ID)


intervention_chemo_ff <- left_join(intervention_chemo |> filter(sample_ID %in% unique(data_ff$sample_ID)) |> distinct(), d_minmax, by = "sample_ID")
intervention_chemo_ff <- left_join(intervention_chemo_ff, data_ff_label, by = "sample_ID")
intervention_chemo_ff <- intervention_chemo_ff |> select(-sample_ID)

intervention_other_ff <- left_join(intervention_other |> filter(sample_ID %in% unique(data_ff$sample_ID)) |> distinct(), d_minmax, by = "sample_ID")
intervention_other_ff <- left_join(intervention_other_ff, data_ff_label, by = "sample_ID")
intervention_other_ff <- intervention_other_ff |> select(-sample_ID)

# plot
options(repr.plot.width=20, repr.plot.height=25)
p_nonrecur <- ggplot(data = data_ff) + 

    # Clinical data
    geom_vline(data = NULL, aes(xintercept=0, color = "surgery"), linewidth = 1.1, linetype = 8) +      # Surgery
    scale_color_manual(name="Clinical timepoints", values=colors1) + # 
    new_scale_color() +

    # TF estimates
    geom_ribbon(aes(x = timepoint, ymin = TF_lCI_kmers, ymax = TF_uCI_kmers), alpha = 0.2, color="grey") + 
    geom_point(aes(x = timepoint, y = TF_kmers, color = ctDNA_detected_kmers), size = 2) + #
    scale_color_manual(name = "ctDNA detected", values=colors2) + 
    new_scale_color() +

    # clonal SNVs mean AF
    
    geom_point(aes(x = timepoint, y = TF_snvs, color = ctDNA_detected_snvs), size = 2) + #
    scale_color_manual(name = "ctDNA detected", values=colors2) + 
    new_scale_color() +

    geom_line(aes(x = timepoint, y = TF_kmers, color = "ctDNA-mers TF"), size = 0.7) + 
    geom_line(aes(x = timepoint, y = TF_snvs, color = "clonal SNVs mean AF"), size = 0.7) + 
    scale_color_manual(name = "TF estimation", values=colors3, limits = c("ctDNA-mers TF", "clonal SNVs mean AF")) + 
    new_scale_color() +

    #imaging
    geom_vline(data = imaging_ff, aes(xintercept = time_since_op_days, color = imaging_res_label), linewidth = 1) + 
    scale_color_manual(name = "Imaging", values=colors_imaging) +
    new_scale_color() +

    # adjuvant chemo
    geom_rect(data = ad_chem_ff,
              aes(xmin = adjuvant_chemo_start_days, xmax = adjuvant_chemo_end_days, ymin = TF_min, ymax = TF_max,
                  fill = "adjuvant \nchemotherapy"), alpha = 0.3) + 

    # intervention administered once
      geom_rect(data = intervention_other_ff,                                             # Intervention
               aes(xmin = intervention_start_days-3, xmax = intervention_start_days+3, ymin = TF_min, ymax = TF_max,
                   fill = "intervention \n(surgery, RFA, MWA, radiotherapy)"), alpha = 0.3) +
    
    # intervention chemo
    geom_rect(data = intervention_chemo_ff,
             aes(xmin = intervention_start_days, xmax = intervention_end_days, ymin = TF_min, ymax = TF_max,
                 fill = "intervention \nchemo- or radiation therapy"), alpha = 0.3) + 
    scale_fill_manual(name = "Treatment", values = colors_intervention) + 
    
    theme_bw() +
    theme(text = element_text(size = 11), 
          legend.box.margin = margin(6, 6, 6, 6), #116, 6, 6, 6
          legend.position = "right", 
          legend.box = "vertical",
          legend.spacing.y = unit(-0.1, 'cm')) + 
    xlab("Sample time point (days since surgery)") +
    facet_wrap(~pt_label, ncol = 3, scales = "free_y") + 
    ylab("")

sink()

png(filename = snakemake@output[["recur_pt_timelines"]], width = 1000, height = 1000)
print(p_recur)
dev.off()


png(filename = snakemake@output[["no_recur_pt_timelines"]], width = 1200, height = 2000)
print(p_nonrecur)
dev.off()
