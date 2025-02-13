sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(coda)
library(rstan)
rstan_options(auto_write = TRUE)

######################### Input #########################

#### input models
mod = snakemake@input[["mod"]]

#### input parameters

## pre-treatment/baseline (ctDNA-positive) cfDNA sample parameter estimates
preop_estimates <- read.csv(snakemake@input[["preop_est"]])

## tumor component
TF_prior_beta_b <- as.integer(snakemake@params[["TF_prior_beta_b"]])
print(paste0("Tumor fraction prior Beta distribution b parameter (sets the prior sample size): ", TF_prior_beta_b))

wt_prior_n <- 100
print(paste0("Tumor component weight sample size: ", wt_prior_n))

wt_prior_n_scale <- 1/1
print(paste0("Tumor component weight sample size scale (used in case the prior sample size is made dependent on input data size with wt_prior_n = n_UT_scale): ", wt_prior_n_scale))

t_phi_lb_scale <- as.numeric(snakemake@params[["t_phi_lb_scale"]])
print(paste0("Tumor component phi parameter lower bound scaling factor: ", t_phi_lb_scale))

wT_mean = preop_estimates$wt_mean[1]
print(paste0("Tumor component weight prior mean (baseline cfDNA sample weight estimate): ", wT_mean))

wT_lb = wT_mean/2
print(paste0("Tumor component weight lower bound (mean between the prior weight mean and 0): ", wT_lb))

## germline component
gl_comp_mean <- 1/2
print(paste0("Germline component mean scaling factor: ", gl_comp_mean))

# germline component weight and variance scaling parameter baseline/pre-treatment (ctDNA-positive sample) estimates
w_gl <- preop_estimates$wgl_mean[1]
print(paste0("Germline component weight: ", w_gl))

gl_phi <- preop_estimates$glphi_mean[1]
print(paste0("Germline component phi parameter: ", gl_phi))

# noise component
f <- snakemake@input[["noise_rate"]]
noise_est <- read.csv(f)
noise_mu = noise_est$mean_mu[1]
noise_phi = noise_est$mean_phi[1]
print(paste0("Noise rate estimates, mu and phi: ", noise_mu, ", ", noise_phi))


#### input data

## kmer counts
f <- snakemake@input[["kmer_data"]]
data = read.table(f)
colnames(data) <- c("kmer", "tumor", "gc_content", "cfDNA")
n_UT <- nrow(data)
print(paste0("Total number of UT k-mers: ", n_UT))

## set tumor component weight prior sample sample size (scaled accordingly, if dependent on the input data size)
if (wt_prior_n == "n_UT"){
    wt_prior_n = n_UT
} else if (wt_prior_n == "n_UT_scale"){
    wt_prior_n = round(n_UT*wt_prior_n_scale)
}
print(paste0("Tumor component weight prior distribution sample size: ", wt_prior_n))
n_wt_p = wt_prior_n

## cfDNA mean count
f <- snakemake@input[["cfDNA_mean"]]
cfDNA_mean_df <- read.table(f, header = T, sep = ",") 
cfDNA_mean  <- cfDNA_mean_df$mean[1]

######################### 3-component mixture model #########################

# combine input data and parameters
data_list = list(n = n_UT,
                 n_wt_p = n_wt_p,
                 noise_mu = noise_mu, 
                 noise_phi = noise_phi,
                 cfDNA_mean = cfDNA_mean, 
                 c_cfDNA = data$cfDNA,
                 max_count_cfDNA_p1 = max(data$cfDNA) + 1,
                 t_phi_lb_scale = t_phi_lb_scale,
                 wT_mean = wT_mean,
                 wT_lb = wT_lb, 
                 w_gl = w_gl, 
                 gl_phi = gl_phi, 
                 TF_prior_beta_b = TF_prior_beta_b, 
                 gl_comp_mean = gl_comp_mean, 
                 t_phi_a = 1, 
                 t_phi_b = 1)

# set initial values
initf1 <- function() {list("TF" = 1e-5,
                           "t_phi" = runif(1, (cfDNA_mean*1e-5)**2/((cfDNA_mean*1e-5*t_phi_lb_scale) - (cfDNA_mean*1e-5)), 100),
                           "w_t" = runif(1, wT_lb, 1))} 


# posterior sampling
set.seed(1)
start_time <- Sys.time()
out <- rstan::stan(file = mod, 
        data = data_list,
        chains = 4,
        include = TRUE, pars = c("TF", "w_t", "t_phi", "gl_assign", "t_assign", "n_assign"), 
        iter = 2500,
        warmup = 500,
        cores = 4, 
        seed = 1, 
        init = initf1)

end_time <- Sys.time()
print(paste0("Time used for burn in and posterior sampling: ", end_time - start_time))

print("Model summary")
sum_mod <- summary(out, pars = c("TF", "w_t", "t_phi", "lp__"))
sum_mod <- sum_mod$summary
print(sum_mod)

TF_n_eff = sum_mod["TF", "n_eff"]
tphi_n_eff = sum_mod["t_phi", "n_eff"]

# handle convergence issues by setting an upper bound to t_phi
if (TF_n_eff < 1000 | tphi_n_eff < 1000){
    print("Effective sample size of TF and/or tphi too low, rerunning model with more informative t_phi prior to ensure convergence")
    print(paste0("TF effective sample size: ", TF_n_eff))
    print(paste0("t_phi effective sample size: ", tphi_n_eff))

    # combine input data and parameters
    data_list = list(n = n_UT,
                    n_wt_p = n_wt_p,
                    noise_mu = noise_mu, 
                    noise_phi = noise_phi,
                    cfDNA_mean = cfDNA_mean, 
                    c_cfDNA = data$cfDNA,
                    max_count_cfDNA_p1 = max(data$cfDNA) + 1,
                    t_phi_lb_scale = t_phi_lb_scale, 
                    wT_mean = wT_mean,
                    wT_lb = wT_lb, 
                    w_gl = w_gl, 
                    gl_phi = gl_phi, 
                    TF_prior_beta_b = TF_prior_beta_b, 
                    gl_comp_mean = gl_comp_mean, 
                    t_phi_a = 1, 
                    t_phi_b = 100)

    # set initial values
    initf1 <- function() {list("TF" = 1e-5,
                            "t_phi" = runif(1, (cfDNA_mean*1e-5)**2/((cfDNA_mean*1e-5*t_phi_lb_scale) - (cfDNA_mean*1e-5)), 100),
                            "w_t" = runif(1, wT_lb, 1))} 


    # posterior sampling
    set.seed(1)
    start_time <- Sys.time()
    out <- rstan::stan(file = mod, 
            data = data_list,
            chains = 4,
            include = TRUE, pars = c("TF", "w_t", "t_phi", "gl_assign", "t_assign", "n_assign"), 
            iter = 2500,
            warmup = 500,
            cores = 4, 
            seed = 1, 
            init = initf1)

    end_time <- Sys.time()
    print(paste0("Time used for burn in and posterior sampling (model with a strict t_phi upper bound): ", end_time - start_time))

    print("Model summary")
    sum_mod <- summary(out, pars = c("TF", "w_t", "t_phi", "lp__"))
    sum_mod <- sum_mod$summary
    print(sum_mod)

}

# post processing
print("Extracting MCMC samples")
list_of_draws <- extract(out, pars = c("gl_assign", "t_assign", "n_assign"))
flush.console()

print("Model summary statistics for all parameters")
mean_estimate_tf <- sum_mod["TF", "mean"][[1]]
lower_CI_tf <- sum_mod["TF", "2.5%"][[1]]
upper_CI_tf <- sum_mod["TF", "97.5%"][[1]]
n_eff_tf <- sum_mod["TF", "n_eff"][[1]]
r_hat_tf <- sum_mod["TF", "Rhat"][[1]]
print("TF: Mean, 2.5% and 97.5% quantiles, n_eff, r_hat")
print(c(mean_estimate_tf, lower_CI_tf, upper_CI_tf, n_eff_tf, r_hat_tf))

mean_estimate_wt <- sum_mod["w_t", "mean"][[1]]
lower_CI_wt <- sum_mod["w_t", "2.5%"][[1]]
upper_CI_wt <- sum_mod["w_t", "97.5%"][[1]]
n_eff_wt <- sum_mod["w_t", "n_eff"][[1]]
r_hat_wt <- sum_mod["w_t", "Rhat"][[1]]
print("w_t: Mean, 2.5% and 97.5% quantiles, n_eff, r_hat")
print(c(mean_estimate_wt, lower_CI_wt, upper_CI_wt, n_eff_wt, r_hat_wt))


mean_estimate_t_phi <- sum_mod["t_phi", "mean"][[1]]
lower_CI_t_phi <- sum_mod["t_phi", "2.5%"][[1]]
upper_CI_t_phi <- sum_mod["t_phi", "97.5%"][[1]]
n_eff_t_phi <- sum_mod["t_phi", "n_eff"][[1]]
r_hat_t_phi <- sum_mod["t_phi", "Rhat"][[1]]
print("t_phi: Mean, 2.5% and 97.5% quantiles, n_eff, r_hat")
print(c(mean_estimate_t_phi, lower_CI_t_phi, upper_CI_t_phi, n_eff_t_phi, r_hat_t_phi))

# save component assignments
unique_cfDNA_val = sort(unique(data$cfDNA))
cfDNA_assignments_res <- NULL
for (i in 1:length(unique_cfDNA_val)){
    c_cfDNA = unique_cfDNA_val[i]
    n_cfDNA = sum(data$cfDNA == c_cfDNA)
    t_assign = list_of_draws$t_assign[ ,c_cfDNA+1]
    gl_assign = list_of_draws$gl_assign[ ,c_cfDNA+1]
    n_assign = list_of_draws$n_assign[ ,c_cfDNA+1]
    res_i = tibble(count = c_cfDNA, 
                   n_kmers = n_cfDNA,
                   tumor = round(mean(t_assign)), 
                   germline = round(mean(gl_assign)), 
                   noise = round(mean(n_assign)))
    cfDNA_assignments_res <- rbind(cfDNA_assignments_res, res_i)
}
print("Mean component assignments, df head")
print(head(cfDNA_assignments_res))

######################### Write results #########################
res_3cat <- tibble(tf_mean = mean_estimate_tf, tf_lower_CI = lower_CI_tf, tf_upper_CI = upper_CI_tf, tf_n_eff = n_eff_tf, tf_r_hat = r_hat_tf,
                   wt_mean = mean_estimate_wt, wt_lower_CI = lower_CI_wt, wt_upper_CI = upper_CI_wt, wt_n_eff = n_eff_wt, wt_r_hat = r_hat_wt,
                   tphi_mean = mean_estimate_t_phi, tphi_lower_CI = lower_CI_t_phi, tphi_upper_CI = upper_CI_t_phi, tphi_n_eff = n_eff_t_phi, tphi_r_hat = r_hat_t_phi)
write.csv(res_3cat, snakemake@output[["estimates"]], row.names = FALSE)
write.csv(cfDNA_assignments_res, snakemake@output[["cfDNA_assignments_df"]], row.names = FALSE)
sink()

png(filename = snakemake@output[["traceplot"]])
rstan::traceplot(out, inc_warmup = FALSE, nrow = 4, pars = c("TF", "w_t", "t_phi", "lp__"))
dev.off()

png(filename = snakemake@output[["autocorr"]])
stan_ac(out, pars = c("TF", "w_t", "t_phi"), nrow = 3)
dev.off()

png(filename = snakemake@output[["density"]])
stan_hist(out, nrow = 3, pars = c("TF", "w_t", "t_phi"))
dev.off()

sink(snakemake@output[["summary_txt"]])
print(sum_mod)
sink()