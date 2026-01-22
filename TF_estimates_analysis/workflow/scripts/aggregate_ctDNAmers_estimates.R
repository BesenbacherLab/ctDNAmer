sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)


# input data
f <- snakemake@input[["units"]]
d <- read.csv(f)

cohort = snakemake@params[["cohort"]]
TF_subf = snakemake@params[["TF_subfd"]]
res_subf = snakemake@params[["res_subfd"]]
ctDNA_mers_dir = snakemake@params[["ctDNA_mers_dir"]]

preop_wn <- NULL
d_BL <- d |> filter(timepoint <= 0)
for (row in 1:nrow(d_BL)){
    pt = d_BL[row, "sample_ID"]
    fd = d_BL[row, "cfDNA_ID"]
    
    est3 <- read.csv(paste0(ctDNA_mers_dir, pt, "/", res_subf, fd, "/", TF_subf, "preop_tf_estimate.csv"))
    preop_wn <- rbind(preop_wn, tibble(sample_ID = pt, 
                                       wgl = est3$wgl_mean, 
                                       glphi = est3$glphi_mean))
}


# get TF estimates
d[c("cfDNA_mean", "TF", "TF_uCI", "TF_lCI", "wt", "tphi", "wgl", "glphi", "wn")] <- NA

for (row in 1:nrow(d)){
    pt = d[row, "sample_ID"][1]
    fd = d[row, "cfDNA_ID"][1]
    timepoint = d[row, "timepoint"][1]
    preop_wn_pt <- preop_wn |> filter(sample_ID == pt)

    cfDNA_mean <- read.csv(paste0("./results/patients/", pt, "/", fd, "/cfDNA_iGL_int_cfDNAc_mean.csv"), header = TRUE)
    d[row, "cfDNA_mean"] = cfDNA_mean$mean[1]

    f_est <- ifelse(timepoint <= 0, "preop_tf_estimate.csv", "postop_tf_estimate.csv")
    est3 <- read.csv(paste0(ctDNA_mers_dir, pt, "/", res_subf, fd, "/", TF_subf, f_est))

    d[row, c("TF", "TF_uCI", "TF_lCI")] = c(est3$tf_mean, est3$tf_upper_CI, est3$tf_lower_CI)
    d[row, c("wgl", "glphi")] = c(preop_wn_pt$wgl, preop_wn_pt$glphi)
    
    d[row, c("wt", "tphi", "wn")] = c(est3$wt_mean, est3$tphi_mean, 1-est3$wt_mean-preop_wn_pt$wgl)

}

# exclude samples with cfDNA mean below 10 
ex_cfDNA <- d |> group_by(sample_ID) |> filter(cfDNA_mean < 10) |> ungroup()
ex_cfDNA <- ex_cfDNA$cfDNA_ID
print(paste0("Samples with cfDNA mean below 10 (excluded): ", ex_cfDNA))

nrow(d)
d <- d |> filter(cfDNA_mean >= 10)
nrow(d)

# filter for target cohort
d_co = d |> filter(sample_type %in% cohort) 

# summary of samples left
paste0("Number of cfDNA samples: ", dim(d_co)[1])
paste0("Positive samples: ", sum(d_co$ctDNA_status == "pos", na.rm = T))
paste0("Negative samples: ", sum(d_co$ctDNA_status == "neg", na.rm = T))
paste0("Positive and negative samples: ", sum(d_co$ctDNA_status == "pos", na.rm = T) + sum(d_co$ctDNA_status == "neg", na.rm = T))
paste0("NA status samples: ", sum(d_co$ctDNA_status == "unknown", na.rm = T))

write.csv(d_co, snakemake@output[["results"]], row.names = FALSE)

# write excluded cfDNA ID's
d_ex_cfDNA <- tibble(cfDNA_ID = ex_cfDNA)
write.csv(d_ex_cfDNA, snakemake@output[["ex_cfDNA"]], row.names = FALSE)

sink()