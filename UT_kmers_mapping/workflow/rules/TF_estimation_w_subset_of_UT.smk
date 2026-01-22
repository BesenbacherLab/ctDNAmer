#####################################################################################
#                                unique tumor module                                #
#####################################################################################

rule sort_ut_kmers:
    input:
        UT_kmers="results/UT_kmers_mapping/patients/{pt}/UT_filtered_1mm.txt",
    output:
        UT_sorted=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/UT_filtered_sorted.txt",
    resources:
        mem_mb=5000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/unique_tumor/sort_ut_kmers.out"
    conda:
        "../envs/py3_12.yaml"
    shell:
        "sort -k1 -o {output.UT_sorted} {input.UT_kmers}"

#####################################################################################
#                               empirical noise module                              #
#####################################################################################

rule filter_intersection_counts:
    input:
        combined_intersection="results/patients/{pt}/empirical_noise/{donor}/combined_intersection.txt",
        UT_sorted=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/UT_filtered_sorted.txt",
    output:
        filtered_intersection=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/{{donor}}/combined_intersection_filtered.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/empirical_noise/combine_intersection_counts_{{donor}}.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "join --check-order -o 1.1 2.2 2.3 {input.UT_sorted} {input.combined_intersection} > {output.filtered_intersection} 2> {log}"

rule create_intersection_tables:
    input:
        kmers=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/{{donor}}/combined_intersection_filtered.txt",
    output:
        count_table=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/{{donor}}/combined_int_gc_content_count_table.txt",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/empirical_noise/create_intersection_table_{{donor}}.out"
    params:
        k=k(),
        max_count=config["count_UB"],
    conda:
        "../envs/py3_12.yaml"
    script:
        "../scripts/gc_content_count_table_noise_kmers.py"

rule model_empirical_noise:
    input:
        UT_mdata="results/patients/{pt}/unique_tumor/UT_n_and_ci.txt",
        noise_count_tables=expand("results/UT_kmers_mapping/patients/{{pt}}/" + f"{UT_subset}/empirical_noise/" + "{donor}/combined_int_gc_content_count_table.txt", donor=donors.index.tolist()),
        donor_means=expand("results/donors/{donor}/mean_count.csv", donor=donors.index.tolist()),
        model=os.path.join(workflow.basedir, "scripts/models/empirical_noise.stan"),
    output:
        param_estimates=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/estimates.csv",
        traceplot=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/traceplot.png",
        autocorr=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/autocorrelation.png",
        summary_txt=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/modeling_summary.txt",
        param_density=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/empirical_noise/param_density.png",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/empirical_noise/model_empirical_noise.out"
    params:
        donor_list=donors.index.tolist(),
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_empirical_noise_stan.R"

#####################################################################################
#                                ut annotation module                               #
#####################################################################################

rule filter_combined_ut_cfDNA:
    input:
        combined="results/patients/{pt}/{cfDNA_ID}/UT_cfDNA_annotation.txt",
        UT_sorted=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/UT_filtered_sorted.txt",
    output:
        combined=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{cfDNA_ID}}/UT_cfDNA_annotation.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/{{cfDNA_ID}}/filter_combined_ut_cfDNA.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "join --check-order -o 1.1 2.2 2.3 2.4 {input.UT_sorted} {input.combined} > {output.combined} 2> {log}" 

#####################################################################################
#                                TF estimation module                               #
#####################################################################################

rule model_tf_preop_cfDNA:
    input:
        kmer_data=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/UT_cfDNA_annotation.txt",
        cfDNA_mean_gc="results/patients/{pt}/{preop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv",
        noise_rate=lambda wildcards: str("results/UT_kmers_mapping/patients/" + str(wildcards.pt) + f"/{UT_subset}/empirical_noise/estimates.csv"),
        mod=os.path.join(workflow.basedir, "scripts/models/tf_preop.stan"),
    output:
        estimates=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/tf_estimation/preop_tf_estimate.csv",
        cfDNA_assignments_df=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/tf_estimation/preop_cat_assignments.csv",
        traceplot=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/tf_estimation/preop_traceplot.png",
        autocorr=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/tf_estimation/preop_autocorr.png",
        density=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/tf_estimation/preop_param_density.png",
        summary_txt=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/tf_estimation/preop_modeling_summary.txt",
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/{{preop_cfDNA_ID}}/model_tf_preop_cfDNA.out"
    params:
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
        TF_prior_beta_b=config["TF_prior_n"],
        t_phi_lb_scale=config["t_phi_lb_scale"],
        wT_mean=config["baseline_model"]["t_w_mean"],
        w_gl_prior_beta_b=config["baseline_model"]["gl_w_beta_b"],
        wt_prior_n_scale=config["baseline_model"]["t_w_prior_n_scale"], 
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_tf_preop_stan.R"

rule model_tf_postop_cfDNA:
    input:
        kmer_data=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/UT_cfDNA_annotation.txt",
        cfDNA_mean_gc="results/patients/{pt}/{postop_cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv",
        noise_rate=lambda wildcards: str("results/UT_kmers_mapping/patients/" + str(wildcards.pt) + f"/{UT_subset}/empirical_noise/estimates.csv"),
        preop_est= aggregate_preop_estimates,
        mod=os.path.join(workflow.basedir, "scripts/models/tf_postop.stan"),
    output:
        estimates=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/tf_estimation/postop_tf_estimate.csv",
        cfDNA_assignments_df=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/tf_estimation/postop_cat_assignments.csv",
        traceplot=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/tf_estimation/postop_traceplot.png",
        autocorr=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/tf_estimation/postop_autocorr.png",
        density=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/tf_estimation/postop_param_density.png",
        summary_txt=f"results/UT_kmers_mapping/patients/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/tf_estimation/postop_modeling_summary.txt",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        f"logs/UT_kmers_mapping/{{pt}}/{UT_subset}/{{postop_cfDNA_ID}}/model_tf_postop_.out"
    params:
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
        TF_prior_beta_b=config["TF_prior_n"],
        t_phi_lb_scale=config["t_phi_lb_scale"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/model_tf_postop_stan.R"