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

f <- snakemake@input[["donors"]]
d_donor = read.table(f, header = TRUE)

f <- snakemake@input[["units"]]
d_cfDNA = read.table(f, header = TRUE)
d_cfDNA <- left_join(d_cfDNA, d, by = "sample_ID") |> select(-cfDNA)


f <- snakemake@input[["clin_data"]]
clin_data <- read.csv(f)
clin_data <- clin_data |> 
    select(patient_id, is_relapse..0.no.1.yes., adjuvant_chemo_start_days, adjuvant_chemo_end_days) |> 
    rename(relapse = is_relapse..0.no.1.yes., sample_ID = patient_id)
d_cfDNA <-left_join(d_cfDNA, clin_data, by="sample_ID")

# filter for the target cohort
if (cohort != "both"){
    d_cfDNA <- d_cfDNA |> filter(sample_type == cohort)
}

# select non-recurring patients samples that have been obtained after treatment end
d_cfDNA_nr <- d_cfDNA |>
    filter(relapse == 0) |> 
    group_by(sample_ID) |> 
    mutate(adjuvant_chem_end = ifelse(is.na(adjuvant_chemo_end_days), 0, adjuvant_chemo_end_days)) |> 
    filter(timepoint > adjuvant_chem_end) |> ungroup() |> select(-adjuvant_chem_end)

gc_lower <- 20
gc_upper <- 80
options(dplyr.summarise.inform = FALSE)

# calculate mean number of observed UT k-mers in non-recurring patients ctDNA-negative samples
res_signal <- NULL
for (i in 1:nrow(d_cfDNA_nr)){ #
    pt = d_cfDNA_nr$sample_ID[i]
    fd = d_cfDNA_nr$cfDNA_ID[i]
    s_type_i = d_cfDNA_nr$sample_type[i]
    n_s <- nrow(d_cfDNA_nr |> filter(sample_ID == pt))
    
    kmer_data <- read.table(paste0(res_dir, "results/patients/", pt, "/", fd, "/UT_cfDNA_annotation.txt"), header = FALSE)
    colnames(kmer_data) <- c("kmer", "tumor", "gc", "cfDNA")
    kmer_data <- kmer_data |>
        filter(cfDNA > 0, cfDNA <= 100) |>
        group_by(cfDNA) |>
        dplyr::summarize(n = n()) |>
        mutate(pt = pt, 
               cfDNA_f = fd,
               n_samples = n_s,
               sample_type = s_type_i)
    res_signal <- rbind(res_signal, kmer_data)   
}
res_signal_sum <- res_signal |> group_by(pt, cfDNA, sample_type) |> 
    summarise(mean_n_cfDNA = sum(n)/n_samples) |> distinct()
res_signal <- left_join(res_signal, res_signal_sum, by = c("pt", "cfDNA", "sample_type")) |> 
    rowwise() |>
    mutate(var_add = (n - mean_n_cfDNA)**2) |>
    ungroup() |> 
    group_by(pt, cfDNA, sample_type) |> 
    summarise(n_samples = unique(n_samples),
              var = sum(var_add)/(n_samples-1), 
              mean_n_cfDNA = unique(mean_n_cfDNA), 
              sd_signal = sqrt(var)) |>
    ungroup()


# calculate mean number of observed UT k-mers in healthy donors
res_noise <- NULL
for (pt in unique(d_cfDNA_nr$sample_ID)){
    UT_mdata <- read.table(paste0(res_dir, "results/patients/", pt, "/unique_tumor/UT_n_and_ci.txt"), header = TRUE)
    min_Tcount <- UT_mdata$final_lower_cutoff[1]
    max_Tcount <- UT_mdata$final_upper_cutoff[1]
    s_type_i <- d_cfDNA_nr$sample_type[which(d_cfDNA_nr$sample_ID == pt)][1]

    res_pt <- NULL
    for (donor in d_donor$donor_ID){
        
        f <- paste0(res_dir, "results/patients/", pt, "/empirical_noise/", donor, "/combined_int_gc_content_count_table.txt")
        d_pt <- read.table(f, header = TRUE)
        colnames(d_pt) <- c("GC", "UT_count", "cfDNA", "n_kmers")
        d_pt <- d_pt |> 
            filter(between(GC, gc_lower, gc_upper)) |>
            filter(between(UT_count, min_Tcount, max_Tcount)) |> 
            group_by(UT_count, cfDNA) |>
            summarise(n_kmers = sum(n_kmers, na.rm = T)) |> 
            mutate(donor = donor, pt = pt, sample_type = s_type_i) |> 
            ungroup()
        res_pt <- rbind(res_pt, d_pt)
    }
    res_pt_sum <- res_pt |> 
        group_by(donor, cfDNA, pt, sample_type) |> 
        summarise(n_kmers = sum(n_kmers)) |>
        ungroup() |> 
        group_by(cfDNA, pt, sample_type) |>
        summarise(n_total = sum(n_kmers, na.rm = T), 
                  mean_n_donor = n_total/nrow(d_donor)) |>
        ungroup() 
    res_pt <- left_join(res_pt, res_pt_sum, by = c("pt", "cfDNA", "sample_type")) |> 
        rowwise() |>
        mutate(var_add = (n_kmers - mean_n_donor)**2) |>
        ungroup() |> 
        group_by(pt, cfDNA, sample_type) |> 
        summarise(var = sum(var_add)/(nrow(d_donor)-1), 
                  mean_n_donor = unique(mean_n_donor), 
                  sd_noise = sqrt(var)) |> 
        ungroup()
    res_noise <- rbind(res_noise, res_pt)
}

# join data sets
res <- full_join(res_signal |> select(pt, cfDNA, sample_type, mean_n_cfDNA, sd_signal), 
                 res_noise |> select(pt, cfDNA, sample_type, mean_n_donor, sd_noise), by = c("pt", "cfDNA", "sample_type"))

# calculate correlation betwen matched ctDNA-negative samples and unmatched healthy cfDNA samples
cor_table <- res |> 
    filter(cfDNA <= 4) |> 
    filter(!(is.na(mean_n_donor))) |>
    filter(!(is.na(mean_n_cfDNA))) |>
    group_by(cfDNA, sample_type) |> 
    dplyr::summarize(cor = cor(mean_n_cfDNA, mean_n_donor, method = c("pearson"))) |> 
    ungroup()

# plotting
if (cohort %in% c("frfr", "ffpe")){
    options(repr.plot.width=16, repr.plot.height=4)
    
    res <- left_join(res, cor_table, by = c("cfDNA", "sample_type")) |> 
        rowwise() |>
        mutate(label = paste0("cfDNA count: ", cfDNA, "; cor: ", round(cor, 3))) |>
        ungroup()
    
    p <- ggplot(res %>% filter(cfDNA <= 4) |> filter(sample_type == cohort)) + 
        geom_point(aes(x = mean_n_donor, y = mean_n_cfDNA)) + 
        facet_wrap(~label, scales = "free", labeller=label_value, nrow = 1) + 
        geom_abline(linetype = 2) + 
        theme_bw() + 
        theme(text = element_text(size = 18)) + 
        ylab("mean number of k-mers\n in postop cfDNA") + 
        xlab("mean number of k-mers in healthy donor cfDNA")
} else {
    res <- left_join(res, cor_table, by = c("cfDNA", "sample_type")) |> 
        rowwise() |>
        mutate(label = paste0("cfDNA count: ", cfDNA)) |>
        ungroup()
    
    colors = c("frfr" = "darkseagreen", "ffpe" = "indianred")
    options(repr.plot.width=18, repr.plot.height=5)
    p <- ggplot(res %>% filter(cfDNA <= 4)) + 
        geom_point(aes(x = mean_n_donor, y = mean_n_cfDNA, color = sample_type)) + 
        facet_wrap(~label, scales = "free", labeller=label_value, nrow = 1) + 
        scale_color_manual(name = "Sample type", values = colors) + 
        geom_abline(linetype = 2) + 
        theme_bw() + 
        theme(text = element_text(size = 18), 
              legend.position = "bottom") + 
        ylab("mean number of k-mers\n in postop cfDNA") + 
        xlab("mean number of k-mers in healthy donor cfDNA") 
        
}


sink()

png(filename = snakemake@output[["emp_noise_data"]], width = 1200, height = 300)
print(p)
dev.off()