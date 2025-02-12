rule model_tf_preop_cfDNA:
    input:
        kmer_data="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/minct{ct_cutoff}/UT_cfDNA_annotation.txt",
        cfDNA_mean="results/patients/{pt}/{preop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv",
        noise_rate=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor_count_filter_test/empirical_noise/minct" + str(wildcards.ct_cutoff) + "/estimates.csv"),
        mod=os.path.join(workflow.basedir, "scripts/models/tf_preop.stan"),
    output:
        estimates="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_tf_estimate.csv",
        cfDNA_assignments_df="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_cat_assignments.csv",
        traceplot="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_traceplot.png",
        autocorr="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_autocorr.png",
        density="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_param_density.png",
        summary_txt="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/preop_modeling_summary.txt",
    resources:
        mem_mb=100000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        "logs/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/model_tf_preop_cfDNA/minct{ct_cutoff}.out",
    params:
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
        wT_mean=config["baseline_model"]["t_w_mean"],
        wT_lb=config["baseline_model"]["t_w_lb"],
        TF_prior_beta_b=config["TF_prior_n"],
        w_gl_prior_beta_b=config["baseline_model"]["gl_w_beta_b"],
        t_phi_lb_scale=config["t_phi_lb_scale"],
        wt_prior_n_scale=config["baseline_model"]["t_w_prior_n_scale"],
    conda:
        "../../envs/R4_1.yaml"
    script:
        "../../scripts/model_tf_preop_stan.R"


rule model_tf_postop_cfDNA:
    input:
        kmer_data="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/minct{ct_cutoff}/UT_cfDNA_annotation.txt",
        cfDNA_mean="results/patients/{pt}/{postop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv",
        noise_rate=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor_count_filter_test/empirical_noise/minct" + str(wildcards.ct_cutoff) + "/estimates.csv"),
        preop_est=aggregate_count_filtering_test_preop_estimates,
        mod=os.path.join(workflow.basedir, "scripts/models/tf_postop.stan"),
    output:
        estimates="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_tf_estimate.csv",
        cfDNA_assignments_df="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_cat_assignments.csv",
        traceplot="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_traceplot.png",
        autocorr="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_autocorr.png",
        density="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_param_density.png",
        summary_txt ="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/minct{ct_cutoff}/postop_modeling_summary.txt",
    resources:
        mem_mb=100000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        "logs/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/model_tf_postop_cfDNA/minct{ct_cutoff}.out",
    params:
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
        TF_prior_beta_b=config["TF_prior_n"],
        t_phi_lb_scale=config["t_phi_lb_scale"],
    conda:
        "../../envs/R4_1.yaml"
    script:
        "../../scripts/model_tf_postop_stan.R"


rule preop_aggregate_tf_estimates_across_count_filtered_ut_sets:
    input:
        estimates=count_filtered_preop_estimates,
    output:
        combined_estimates="results/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/preop_combined_tf_estimates.csv",
    resources:
        mem_mb=2000,
        runtime=lambda wildcards, attempt: attempt * 60,
    log:
        "logs/patients/{pt}/{preop_cfDNA_ID}/unique_tumor_count_filter_test/preop_aggregate_tf_estimates_across_count_filtered_ut_sets.out",
    conda:
        "../../envs/R4_1.yaml"
    script:
        "../../scripts/combine_count_filtered_ut_sets_tf_estimates.R"


rule postop_aggregate_tf_estimates_across_count_filtered_ut_sets:
    input:
        estimates=count_filtered_postop_estimates,
    output:
        combined_estimates="results/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/tf_estimation/postop_combined_tf_estimates.csv",
    resources:
        mem_mb=2000,
        runtime=lambda wildcards, attempt: attempt * 60,
    log:
        "logs/patients/{pt}/{postop_cfDNA_ID}/unique_tumor_count_filter_test/postop_aggregate_tf_estimates_across_count_filtered_ut_sets.out"
    conda:
        "../../envs/R4_1.yaml"
    script:
        "../../scripts/combine_count_filtered_ut_sets_tf_estimates.R"