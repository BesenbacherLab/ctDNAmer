data{
    int n; // number of UT k-mers
    int n_wt_p; //sample size w_t prior
    real cfDNA_mean[n]; // mean k-mer count in cfDNA sample
    real noise_mu; // noise negbin mu parameter
    real noise_phi; // noise negbin phi parameter
    int c_cfDNA[n]; // k-mer count in cfDNA
    int max_count_cfDNA_p1; // max cfDNA count plus 1, for component assignments saving
    real wT_mean; // tumor component weight prior mean
    real wT_lb; // tumor component weight lower bound
    real w_gl; // germline component weight from 3cat preop model
    real gl_phi; // germline negbin phi parameter from 3cat preop model
    real t_phi_lb_scale; //tumor component phi lower bound scaling factor
    real cfDNA_mean_max; // min cfDNA mean
    int TF_prior_beta_b; //TF beta prior b parameter
    real gl_comp_mean; //germline component mean 
    real t_phi_a; // t_phi gamma prior a parameter
    real t_phi_b; // t_phi gamma prior b parameter
}

parameters{
    real<lower=0, upper=1> TF; // tumor fraction
    real<lower=wT_lb, upper=1> w_t; // tumor component weight
    real<lower=pow((cfDNA_mean_max*TF), 2)/((cfDNA_mean_max*TF*t_phi_lb_scale) - (cfDNA_mean_max*TF))> t_phi;
}

model{
    TF ~ beta(1, TF_prior_beta_b);
    w_t ~ beta(n_wt_p*wT_mean, (n_wt_p - (wT_mean*n_wt_p))) T[wT_lb, ];

    t_phi ~ gamma(t_phi_a, t_phi_b) T[pow((cfDNA_mean_max*TF), 2)/((cfDNA_mean_max*TF*t_phi_lb_scale) - (cfDNA_mean_max*TF)), ];
    
    for (i in 1:n){
        vector[3] contributions;
        contributions[1] = log(w_t) + neg_binomial_2_lpmf(c_cfDNA[i] | TF*cfDNA_mean[i], t_phi);
        contributions[2] = log(w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | gl_comp_mean*cfDNA_mean[i], gl_phi);
        contributions[3] = log(1 - w_t - w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | noise_mu*cfDNA_mean[i], noise_phi);
        
        target += log_sum_exp(contributions);
    }
}

generated quantities {
    vector[n] log_lik;
    vector[n] post_prob_t;
    vector[n] post_prob_gl;
    vector[n] post_prob_n;
    vector[max_count_cfDNA_p1] t_assign;
    vector[max_count_cfDNA_p1] gl_assign;
    vector[max_count_cfDNA_p1] n_assign;
    
    t_assign = rep_vector(0, max_count_cfDNA_p1);
    gl_assign = rep_vector(0, max_count_cfDNA_p1);
    n_assign = rep_vector(0, max_count_cfDNA_p1);
    
    for (i in 1:n) {
        vector[3] contr;
        contr[1] = log(w_t) + neg_binomial_2_lpmf(c_cfDNA[i] | TF*cfDNA_mean[i], t_phi);
        contr[2] = log(w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | gl_comp_mean*cfDNA_mean[i], gl_phi);
        contr[3] = log(1 - w_t - w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | noise_mu*cfDNA_mean[i], noise_phi);
        log_lik[i] = log_sum_exp(contr);

        post_prob_t[i] = exp(contr[1] - log_sum_exp(contr)); 
        post_prob_gl[i] = exp(contr[2] - log_sum_exp(contr));
        post_prob_n[i] = exp(contr[3] - log_sum_exp(contr));
        
        if (post_prob_t[i] > post_prob_gl[i] && post_prob_t[i] > post_prob_n[i])
            t_assign[c_cfDNA[i] + 1] += 1;
        else if (post_prob_gl[i] > post_prob_t[i] && post_prob_gl[i] > post_prob_n[i])
            gl_assign[c_cfDNA[i] + 1] += 1;
        else
            n_assign[c_cfDNA[i] + 1] += 1;
    }
}
