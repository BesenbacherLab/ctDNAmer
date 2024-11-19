rule model_tf_preop_cfDNA:
    input:
        kmer_data = "results/patients/{pt}/{preop_cfDNA_ID}/UT_cfDNA_annotation.txt",
        cfDNA_mean = "results/patients/{pt}/{preop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv",
        noise_rate = lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/empirical_noise/estimates.csv"),
        pre_3cat_mod = os.path.join(workflow.basedir, "scripts/models/tf_preop_3cat.stan"),
        pre_2cat_mod = os.path.join(workflow.basedir, "scripts/models/tf_preop_2cat.stan"),
    output:
        estimates = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/3cat/preop_tf_estimate.csv",
        cfDNA_assignments_df = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/3cat/preop_cat_assignments.csv",
        traceplot = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/3cat/preop_traceplot.png",
        autocorr = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/3cat/preop_autocorr.png",
        density = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/3cat/preop_param_density.png",
        summary_txt = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/3cat/preop_modeling_summary.txt", 
        estimates_2cat = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/2cat/preop_tf_estimate.csv",
        cfDNA_assignments_df_2cat = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/2cat/preop_cat_assignments.csv",
        traceplot_2cat = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/2cat/preop_traceplot.png",
        autocorr_2cat = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/2cat/preop_autocorr.png",
        density_2cat = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/2cat/preop_param_density.png",
        summary_txt_2cat = "results/patients/{pt}/{preop_cfDNA_ID}/tf_estimation/2cat/preop_modeling_summary.txt", 
    resources: 
        mem_mb = 30000, 
        runtime=lambda wildcards, attempt: attempt * 180,
    threads: 4
    log:
        "logs/patients/{pt}/{preop_cfDNA_ID}/model_tf_preop_cfDNA.out",
    params:
        gc_lower = 20,
        gc_upper = 80,
        wT_mean = 0.01,
        wT_lb = 0.001,
        TF_prior_beta_b = 100,
        w_gl_prior_beta_b = 100,
        gl_comp_mean = 0.5,
    conda: 
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_tf_preop_stan.R"


rule model_tf_postop_cfDNA:
    input:
        kmer_data = "results/patients/{pt}/{postop_cfDNA_ID}/UT_cfDNA_annotation.txt",
        cfDNA_mean = "results/patients/{pt}/{postop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv",
        noise_rate = lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/empirical_noise/estimates.csv"),
        preop_est =  aggregate_preop_estimates,
        post_3cat_mod = os.path.join(workflow.basedir, "scripts/models/tf_postop_3cat.stan"),
        post_2cat_mod = os.path.join(workflow.basedir, "scripts/models/tf_postop_2cat.stan"),
    output: 
        estimates = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/3cat/postop_tf_estimate.csv",
        cfDNA_assignments_df = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/3cat/postop_cat_assignments.csv",
        traceplot = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/3cat/postop_traceplot.png",
        autocorr = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/3cat/postop_autocorr.png",
        density = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/3cat/postop_param_density.png",
        summary_txt = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/3cat/postop_modeling_summary.txt", 
        estimates_2cat = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/2cat/postop_tf_estimate.csv",
        cfDNA_assignments_df_2cat = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/2cat/postop_cat_assignments.csv",
        summary_txt_2cat = "results/patients/{pt}/{postop_cfDNA_ID}/tf_estimation/2cat/postop_modeling_summary.txt", 
    resources: 
        mem_mb = 30000, 
        runtime=lambda wildcards, attempt: attempt * 180,
    threads: 4
    log:
        "logs/patients/{pt}/{postop_cfDNA_ID}/model_tf_postop_cfDNA.out",
    params:
        gc_lower = 20,
        gc_upper = 80,
        wT_mean = 0.01,
        wT_lb = 0.001,
        TF_prior_beta_b = 100,
        gl_comp_mean = 0.5, 
    conda: 
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_tf_postop_stan.R"