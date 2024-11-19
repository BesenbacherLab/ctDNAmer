sink(snakemake@log[[1]], append=TRUE)

library(tidyverse)
library(coda)
library(rstan)
rstan_options(auto_write = TRUE)

######################### Input #########################

# input models
post_3cat_mod = snakemake@input[["post_3cat_mod"]]
post_2cat_mod = snakemake@input[["post_2cat_mod"]] 

# input parameters
gc_lower <- as.integer(snakemake@params[["gc_lower"]])
gc_upper <- as.integer(snakemake@params[["gc_upper"]])
wT_mean <- as.numeric(snakemake@params[["wT_mean"]])
wT_lb <- as.numeric(snakemake@params[["wT_lb"]])
TF_prior_beta_b <- as.integer(snakemake@params[["TF_prior_beta_b"]])
gl_comp_mean <- as.numeric(snakemake@params[["gl_comp_mean"]])

# germline component weight and variance scaling parameter
preop_estimates <- read.csv(snakemake@input[["preop_est"]])
w_gl <- preop_estimates$wgl_mean[1]
print(paste0("w_gl: ", w_gl))

gl_phi <- preop_estimates$glphi_mean[1]
print(paste0("gl_phi: ", gl_phi))

# input data - kmers
f <- snakemake@input[["kmer_data"]]
data = read.table(f)
colnames(data) <- c("kmer", "tumor", "gc_content", "cfDNA")
n_UT <- nrow(data)
print(paste0("Total number of UT k-mers: ", n_UT))

# input data - cfDNA mean
f <- snakemake@input[["cfDNA_mean"]]
cfDNA_mean_count <- read.table(f, header = T, sep = ",") 
cfDNA_mean_count <- cfDNA_mean_count |> 
    filter(between(as.numeric(gc_content), gc_lower, gc_upper)) |> 
    select(gc_content, mean, var)
print("Head cfDNA mean")
print(head(cfDNA_mean_count))
data <- left_join(data, cfDNA_mean_count, by = c("gc_content"))

# input data - noise rate
f <- snakemake@input[["noise_rate"]]
noise_est <- read.csv(f)
noise_mu = noise_est$mean_mu[1]
noise_phi = noise_est$mean_phi[1]
print(paste0("Noise rate estimates, mu and phi: ", noise_mu, ", ", noise_phi))

# var = mean + (mean**2)/phi; phi = (mean**2)/(var - mean) --> larger phi, smaller variance, 
# using cfDNA mean and var for phi lower bound calculation
cfDNA_max_mean <- cfDNA_mean_count %>% filter(mean == max(mean))
cfDNA_max_mean_m <- cfDNA_max_mean$mean[1]
cfDNA_max_mean_v <- cfDNA_max_mean$var[1]

t_phi_lb = (cfDNA_max_mean_m**2)/(cfDNA_max_mean_v-cfDNA_max_mean_m)
print(paste0("t phi lower bound: ", t_phi_lb)) 

######################### 3 cat model #########################

# modeling
set.seed(1)
data_list = list(n = n_UT,
                 noise_mu = noise_mu, 
                 noise_phi = noise_phi,
                 cfDNA_mean = data$mean, 
                 c_cfDNA = data$cfDNA,
                 max_count_cfDNA_p1 = max(data$cfDNA) + 1, 
                 wT_mean = wT_mean,
                 wT_lb = wT_lb, 
                 w_gl = w_gl, 
                 gl_phi = gl_phi, 
                 t_phi_lb = t_phi_lb,
                 TF_prior_beta_b = TF_prior_beta_b, 
                 gl_comp_mean = gl_comp_mean)

initf1 <- function() {list("TF" = rbeta(1, 1, TF_prior_beta_b),
                           "t_phi" = rexp(1, 1) + t_phi_lb,
                           "w_t" = rbeta(1, n_UT*wT_mean, n_UT-(n_UT*wT_mean)))} 

set.seed(1)
start_time <- Sys.time()
out <- rstan::stan(file = post_3cat_mod, 
        data = data_list,
        chains = 4,
        include = TRUE, pars = c("TF", "w_t", "t_phi", "log_lik", "gl_assign", "t_assign", "n_assign"), #
        iter = 2500,
        warmup = 500, 
        cores = 4, 
        seed = 1, 
        init = initf1)

end_time <- Sys.time()
print(end_time - start_time)

# post processing
print("Extracting MCMC samples")
list_of_draws <- extract(out, pars = c("TF", "w_t", "t_phi", "log_lik", "gl_assign", "t_assign", "n_assign"))
flush.console()

print("Model summary")
sum_mod <- summary(out, pars = c("TF", "w_t", "t_phi", "log_lik[1]", "lp__"))
sum_mod <- sum_mod$summary
print(sum_mod)

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


######################### 2 cat model #########################

# modeling
set.seed(1)
data_list = list(n = n_UT, 
                 noise_mu = noise_mu, 
                 noise_phi = noise_phi,
                 cfDNA_mean = data$mean, 
                 c_cfDNA = data$cfDNA,
                 max_count_cfDNA_p1 = max(data$cfDNA) + 1, 
                 w_gl = w_gl, 
                 gl_phi = gl_phi, 
                 gl_comp_mean = gl_comp_mean)


set.seed(1)
start_time <- Sys.time()
out2 <- rstan::stan(file = post_2cat_mod, 
        algorithm="Fixed_param",
        data = data_list,
        chains = 4,
        include = TRUE, pars = c("log_lik", "gl_assign", "n_assign"), #
        iter = 2500,
        warmup = 500, 
        cores = 4, 
        seed = 1)

end_time <- Sys.time()
print(end_time - start_time)

# post processing
print("Extracting MCMC samples")
list_of_draws2 <- extract(out2, pars = c("log_lik", "gl_assign", "n_assign"))
flush.console()

print("Model summary")
sum_mod2 <- summary(out2, pars = c("log_lik[1]", "lp__"))
sum_mod2 <- sum_mod2$summary
print(sum_mod2)

# save component assignments
unique_cfDNA_val = sort(unique(data$cfDNA))
cfDNA_assignments_res2 <- NULL
for (i in 1:length(unique_cfDNA_val)){
    c_cfDNA = unique_cfDNA_val[i]
    n_cfDNA = sum(data$cfDNA == c_cfDNA)
    
    gl_assign = list_of_draws2$gl_assign[ ,c_cfDNA+1]
    n_assign = list_of_draws2$n_assign[ ,c_cfDNA+1]
    res_i = tibble(count = c_cfDNA, 
                   n_kmers = n_cfDNA,
                   germline = round(mean(gl_assign)), 
                   noise = round(mean(n_assign)))
    cfDNA_assignments_res2 <- rbind(cfDNA_assignments_res2, res_i)
}
print("Mean component assignments, df head")
print(head(cfDNA_assignments_res2))

######################### Model comparison #########################
log_lik1_sum <- sum(list_of_draws$log_lik)
log_lik2_sum <- sum(list_of_draws2$log_lik)
ll_diff = log_lik1_sum - log_lik2_sum

print("Log lik 1")
print(log_lik1_sum)
print("Log lik 2")
print(log_lik2_sum)
print("Log lik diff")
print(ll_diff)

######################### Write results #########################
res_3cat <- tibble(tf_mean = mean_estimate_tf, tf_lower_CI = lower_CI_tf, tf_upper_CI = upper_CI_tf, tf_n_eff = n_eff_tf, tf_r_hat = r_hat_tf,
                   wt_mean = mean_estimate_wt, wt_lower_CI = lower_CI_wt, wt_upper_CI = upper_CI_wt, wt_n_eff = n_eff_wt, wt_r_hat = r_hat_wt,
                   tphi_mean = mean_estimate_t_phi, tphi_lower_CI = lower_CI_t_phi, tphi_upper_CI = upper_CI_t_phi, tphi_n_eff = n_eff_t_phi, tphi_r_hat = r_hat_t_phi, 
                   log_lik_sum = log_lik1_sum)
write.csv(res_3cat, snakemake@output[["estimates"]], row.names = FALSE)
write.csv(cfDNA_assignments_res, snakemake@output[["cfDNA_assignments_df"]], row.names = FALSE)

res_2cat <- tibble(log_lik_sum = log_lik2_sum)
write.csv(res_2cat, snakemake@output[["estimates_2cat"]], row.names = FALSE)
write.csv(cfDNA_assignments_res2, snakemake@output[["cfDNA_assignments_df_2cat"]], row.names = FALSE)
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

sink(snakemake@output[["summary_txt_2cat"]])
print(sum_mod2)
sink()