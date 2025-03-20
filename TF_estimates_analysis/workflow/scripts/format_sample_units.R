sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)


# input data
f <- snakemake@input[["samples"]]
d = read.table(f, header = TRUE)
d <- d |> rowwise() |> mutate(sample_type = ifelse(grepl("frfr", tumor), "frfr", "ffpe")) |> ungroup() |> 
    select(sample_ID, sample_type)


f <- snakemake@input[["units"]]
d_cfDNA = read.table(f, header = TRUE)
d_cfDNA <- left_join(d_cfDNA, d, by = "sample_ID")
d_cfDNA <- d_cfDNA |> select(-cfDNA)


# add clinical relapse and adjucvant chemo
f <- snakemake@input[["clin_data"]] 
clin_data <- read.csv(f)
clin_data <- clin_data |> 
    select(patient_id, is_relapse..0.no.1.yes., time_to_relapse_days, adjuvant_chemo_start_days, adjuvant_chemo_end_days) |> 
    rename(relapse = is_relapse..0.no.1.yes., sample_ID = patient_id)
d_cfDNA <-left_join(d_cfDNA, clin_data, by="sample_ID")


# get recurrence treatment info
f <- snakemake@input[["interventions"]]
d_intv = read.csv(f)
d_intv <- d_intv |>
    rename(sample_ID = patient_id) |>
    select(sample_ID, intervention_start_days, intervention_end_days, intervention_type) |> 
    filter(sample_ID %in% unique(d_cfDNA$sample_ID))

d_intv_last_start = d_intv |>
    group_by(sample_ID) |>
    filter(intervention_start_days == max(intervention_start_days)) |>
    select(-c(intervention_end_days, intervention_type)) |>
    rename(start_of_last_intervention = intervention_start_days) |> distinct()
d_cfDNA <- left_join(d_cfDNA, d_intv_last_start, by="sample_ID")


# add sample ground-truth labels
d_labels <- NULL
for (pt in unique(d_cfDNA$sample_ID)){

    d_pt <- d_cfDNA |> filter(sample_ID == pt) |> arrange(timepoint)
    relapse_pt <- d_pt$relapse[1]
    if (relapse_pt == 1){
        max_r_known_pt = max(d_pt$start_of_last_intervention)
        ad_chem_end_pt = d_pt$adjuvant_chemo_end_days[1]
        last_data_p <- d_pt$timepoint[nrow(d_pt)]
        if (is.na(max_r_known_pt)){
            max_r_known_pt = last_data_p
        }
    } else {
        ad_chem_end_pt = d_pt$adjuvant_chemo_end_days[1]
    }
    for (row in 1:nrow(d_pt)){
        samp_time <- d_pt$timepoint[row]

        if (samp_time <= 0){
            ctDNA_i = "pos"
        } else {
            if (relapse_pt == 0) {
                if (is.na(ad_chem_end_pt)){
                    ctDNA_i = "neg"
                } else {
                    if (samp_time > ad_chem_end_pt){
                        ctDNA_i = "neg"
                    } else {
                        ctDNA_i = "unknown"
                    }    
                }
            } else {
                if (is.na(ad_chem_end_pt) | samp_time > ad_chem_end_pt){
                    if (samp_time <= max_r_known_pt){
                        ctDNA_i = "pos"
                    } else {
                        ctDNA_i = "unknown"
                    }
                } else {
                    ctDNA_i = "unknown"
                }    
            }
        }
        d_pt_row <- d_pt[row, ] |> mutate(ctDNA_status = ctDNA_i)
        d_labels <- rbind(d_labels, d_pt_row)
    }
}

paste0("Number of cfDNA samples: ", dim(d_labels)[1])
paste0("Positive samples: ", sum(d_labels$ctDNA_status == "pos", na.rm = T))
paste0("Negative samples: ", sum(d_labels$ctDNA_status == "neg", na.rm = T))
paste0("Positive and negative samples: ", sum(d_labels$ctDNA_status == "pos", na.rm = T) + sum(d_labels$ctDNA_status == "neg", na.rm = T))
paste0("NA status samples: ", sum(d_labels$ctDNA_status == "unknown", na.rm = T))


write.csv(d_labels, snakemake@output[["formatted_units"]], row.names = FALSE)

sink()