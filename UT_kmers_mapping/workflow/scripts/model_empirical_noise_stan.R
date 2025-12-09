sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(coda)
library(rstan)
rstan_options(auto_write = TRUE)

# input parameters
donor_list <- snakemake@params[["donor_list"]]
gc_lower <- as.integer(snakemake@params[["gc_lower"]])
gc_upper <- as.integer(snakemake@params[["gc_upper"]])
print(paste0("GC content cutoffs: ", gc_lower, ", ", gc_upper))

# input data
UT_mdata = read.table(snakemake@input[["UT_mdata"]], header = TRUE)
n_UT <- UT_mdata$nUT_final[1]
print(paste0("Total number of UT k-mers: ", n_UT))

min_Tcount <- UT_mdata$final_lower_cutoff[1]
max_Tcount <- UT_mdata$final_upper_cutoff[1]
print(paste0("Minimum and maximum tumor count: ", min_Tcount, ", ", max_Tcount))

noise_count_tables_file_list <- snakemake@input[["noise_count_tables"]]
print("List of noise count tables")
print(noise_count_tables_file_list)
donor_means_file_list <- snakemake@input[["donor_means"]]

########## 1. Prepare data for modeling ##########
comb_donor_data <- NULL
for (donor in donor_list){
    data_i = read.table(noise_count_tables_file_list[grepl(donor, noise_count_tables_file_list, fixed = TRUE)], header = TRUE, sep = "\t")
    mean_i = read.csv(donor_means_file_list[grepl(donor, donor_means_file_list, fixed = TRUE)], header = TRUE, sep = ",")

    data_i <- data_i |> 
        filter(between(GC, gc_lower, gc_upper)) |> # GC content filtering
        filter(between(UT_count, min_Tcount, max_Tcount)) # count filtering
    if (nrow(data_i) == 0){ # in case of no noise k-mers, create a pseudo set with one noise k-mer
        print("No noise k-mers observed in this donor, adding a pseudocount")
        data_i <- tibble(GC = 50, 
                         UT_count = min_Tcount, 
                         donor_count = 1, 
                         n_kmers = 1)
    }
    data_i_long <- rep(data_i$donor_count, data_i$n_kmers)
    data_i_long <- c(data_i_long, rep(0, n_UT - sum(data_i$n_kmers)))
    
    f_data_i <- tibble(donor_mean = mean_i$mean[1], 
                       count = data_i_long)
    comb_donor_data <- rbind(comb_donor_data, f_data_i)
}


########## 2. Set up and specify the model ##########
# combine input data
data_list = list(n = nrow(comb_donor_data),
                 mean_count = comb_donor_data$donor_mean, 
                 donor_count = comb_donor_data$count)

# set initial values
initf1 <- function() {
      list("mean_noise" = 0.01, "phi_noise" = 0.1)
    }

########## 3. Run the MCMC sampler ##########
set.seed(1)
start_time <- Sys.time()
out <- rstan::stan(file = snakemake@input[["model"]],
        data = data_list,
        chains = 4,
        include = TRUE, 
        pars = c("mean_noise", "phi_noise"), 
        iter = 1100,
        warmup = 100, 
        cores = 4, 
        seed = 1 , 
        init = initf1)
end_time <- Sys.time()
print(paste0("Time used for burn in and posterior sampling: ", end_time - start_time))
flush.console()

########## 4. Post processing ##########
print(out, pars=c("mean_noise", "phi_noise", "lp__"), probs=c(.1,.5,.9))

# model summary
print("Model summary")
fit_summary <- summary(out, probs = c(0.025, 0.5, 0.975))
sum_mod <- fit_summary$summary
print(sum_mod)

# parameter estiamtes with 95% CI's
print("Model summary statistics for the mean and phi parameters")
mean_estimate_mu <- sum_mod["mean_noise", "mean"][[1]]
lower_CI_mu <- sum_mod["mean_noise", "2.5%"][[1]]
upper_CI_mu <- sum_mod["mean_noise", "97.5%"][[1]]
print("mean_noise: Mean, 2.5% and 97.5% quantiles")
print(c(mean_estimate_mu, lower_CI_mu, upper_CI_mu))

mean_estimate_phi <- sum_mod["phi_noise", "mean"][[1]]
lower_CI_phi <- sum_mod["phi_noise", "2.5%"][[1]]
upper_CI_phi <- sum_mod["phi_noise", "97.5%"][[1]]
print("phi_noise: Mean, 2.5% and 97.5% quantiles")
print(c(mean_estimate_phi, lower_CI_phi, upper_CI_phi))

# effective sample size
n_eff_mu <- sum_mod["mean_noise", "n_eff"][[1]]
n_eff_phi <- sum_mod["phi_noise", "n_eff"][[1]]
print("Effective sample Size: mean_noise, phi_noise")
print(c(n_eff_mu, n_eff_phi))

# Rhat values
r_hat_mu <- sum_mod["mean_noise", "Rhat"][[1]]
r_hat_phi <- sum_mod["phi_noise", "Rhat"][[1]]
print("Rhat: mean_noise, phi_noise")
print(c(r_hat_mu, r_hat_phi))

########## 5. write results and convergence diagnostics ##########
res <- tibble(mean_mu = mean_estimate_mu, 
              lower_CI_mu = lower_CI_mu, 
              upper_CI_mu = upper_CI_mu, 
              n_eff_mu = n_eff_mu, 
              r_hat_mu = r_hat_mu,
              mean_phi = mean_estimate_phi, 
              lower_CI_phi = lower_CI_phi, 
              upper_CI_phi = upper_CI_phi, 
              n_eff_phi = n_eff_phi, 
              r_hat_phi = r_hat_phi)

p_trace <- rstan::traceplot(out, inc_warmup = TRUE, nrow = 2, pars = c("mean_noise", "phi_noise"))
p_ac <- stan_ac(out, pars = c("mean_noise", "phi_noise"), nrow = 2)
p_dens <- stan_hist(out, nrow = 2, pars = c("mean_noise", "phi_noise"))

write.csv(res, snakemake@output[["param_estimates"]], row.names = FALSE)
sink()

png(filename = snakemake@output[["traceplot"]])
print(p_trace)
dev.off()

png(filename = snakemake@output[["autocorr"]])
print(p_ac)
dev.off()

png(filename = snakemake@output[["param_density"]])
print(p_dens)
dev.off()

sink(snakemake@output[["summary_txt"]])
print(sum_mod)
sink()
