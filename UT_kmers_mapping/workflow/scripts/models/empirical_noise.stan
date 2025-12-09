data{
    int n; // number of observations -> n_UT*n_donors
    vector[n] mean_count; // mean k-mer count of donor sample
    int donor_count[n]; // donor counts
}
parameters{
    real<lower=0> mean_noise; // noise distribution mean
    real<lower=0> phi_noise; // noise distribution phi
}
model{
    mean_noise ~ beta(1, 30); // mean Prior
    phi_noise ~ exponential(1); // phi Prior

    donor_count ~ neg_binomial_2(mean_noise*mean_count, phi_noise);
}
