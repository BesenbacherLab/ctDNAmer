rule model_tf_preop_cfDNA:
    input:
        kmer_data="results/patients/{pt}/{preop_cfDNA_ID}/UT_cfDNA_annotation.txt",
        cfDNA_mean="results/patients/{pt}/{preop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean.csv",
        noise_rate=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/empirical_noise/estimates.csv"),
        mod=os.path.join(workflow.basedir, "scripts/models/tf_preop.stan"),
    output:
        estimates="results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/preop_tf_estimate.csv",
        cfDNA_assignments_df="results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/preop_cat_assignments.csv",
        traceplot="results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/preop_traceplot.png",
        autocorr="results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/preop_autocorr.png",
        density="results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/preop_param_density.png",
        summary_txt="results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/preop_modeling_summary.txt",
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        "logs/patients/{pt}/{preop_cfDNA_ID}/model_tf_preop_cfDNA.out"
    params:
        wT_mean=config["baseline_model"]["t_w_mean"],
        wT_lb=config["baseline_model"]["t_w_lb"],
        TF_prior_beta_b=config["TF_prior_n"],
        w_gl_prior_beta_b=config["baseline_model"]["gl_w_beta_b"],
        t_phi_lb_scale=config["t_phi_lb_scale"],
        wt_prior_n_scale=config["baseline_model"]["t_w_prior_n_scale"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_tf_preop_stan.R"


rule model_tf_postop_cfDNA:
    input:
        kmer_data="results/patients/{pt}/{postop_cfDNA_ID}/UT_cfDNA_annotation.txt",
        cfDNA_mean="results/patients/{pt}/{postop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean.csv",
        noise_rate=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/empirical_noise/estimates.csv"),
        preop_est= aggregate_preop_estimates,
        mod=os.path.join(workflow.basedir, "scripts/models/tf_postop.stan"),
    output:
        estimates="results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/postop_tf_estimate.csv",
        cfDNA_assignments_df="results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/postop_cat_assignments.csv",
        traceplot="results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/postop_traceplot.png",
        autocorr="results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/postop_autocorr.png",
        density="results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/postop_param_density.png",
        summary_txt="results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/postop_modeling_summary.txt",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        "logs/patients/{pt}/{postop_cfDNA_ID}/model_tf_postop_.out"
    params:
        TF_prior_beta_b=config["TF_prior_n"],
        t_phi_lb_scale=config["t_phi_lb_scale"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_tf_postop_stan.R"

