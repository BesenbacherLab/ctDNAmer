data{
    int n; // number of UT k-mers
    real cfDNA_mean[n]; // mean k-mer count in cfDNA sample
    real noise_mu; // noise negbin mu parameter
    real noise_phi; // noise negbin phi parameter
    int c_cfDNA[n]; // k-mer count in cfDNA
    int max_count_cfDNA_p1; // max cfDNA count plus one, for component assignment saving
    real w_gl; // germline component weight from 3cat preop model
    real gl_phi; // germline negbin phi parameter from 3cat preop model
    real gl_comp_mean; //germline component mean
}

model{
    for (i in 1:n){
        vector[2] contributions;
        contributions[1] = log(w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | gl_comp_mean*cfDNA_mean[i], gl_phi);
        contributions[2] = log(1 - w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | noise_mu*cfDNA_mean[i], noise_phi);
        
        target += log_sum_exp(contributions);
    }
}

generated quantities {
    vector[n] log_lik;
    vector[n] post_prob_gl;
    vector[n] post_prob_n;
    vector[max_count_cfDNA_p1] gl_assign;
    vector[max_count_cfDNA_p1] n_assign;
    
    gl_assign = rep_vector(0, max_count_cfDNA_p1);
    n_assign = rep_vector(0, max_count_cfDNA_p1);
    
    for (i in 1:n) {
        vector[2] contr;
        contr[1] = log(w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | gl_comp_mean*cfDNA_mean[i], gl_phi);
        contr[2] = log(1 - w_gl) + neg_binomial_2_lpmf(c_cfDNA[i] | noise_mu*cfDNA_mean[i], noise_phi);
        log_lik[i] = log_sum_exp(contr);

        post_prob_gl[i] = exp(contr[1] - log_sum_exp(contr));
        post_prob_n[i] = exp(contr[2] - log_sum_exp(contr));
        
        if (post_prob_gl[i] > post_prob_n[i])
            gl_assign[c_cfDNA[i] + 1] += 1;
        else
            n_assign[c_cfDNA[i] + 1] += 1;
    }
}
