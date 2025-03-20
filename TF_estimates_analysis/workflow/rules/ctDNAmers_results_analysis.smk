

rule plot_TF_estimates_summary: # F2C
    input:
        ctDNAmers_res="results/TF_estimates_analysis/ctDNAmers_TF/units_ctDNAmersTF_ctDNA_det.csv",
    output:
        TF_sum="results/TF_estimates_analysis/plots/ctDNAmers_TF_estimates_summary.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/plot_TF_estimates_summary.out"
    params: 
        cohort = config["cohort"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/plot_TF_estimates_summary.R"


rule plot_lead_time:  #F2E
    input:
        ctDNAmers_res="results/TF_estimates_analysis/ctDNAmers_TF/units_ctDNAmersTF_ctDNA_det.csv",
    output:
        lead_time="results/TF_estimates_analysis/plots/ctDNAmers_lead_time.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/plot_lead_time.out"
    params: 
        cohort = config["cohort"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/plot_lead_time.R"
