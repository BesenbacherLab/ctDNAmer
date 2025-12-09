

rule plot_TF_estimates_summary: # F2C
    input:
        ctDNAmers_res=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF_ctDNA_det.csv",
    output:
        TF_sum=f"{pref}/TF_estimates_analysis/plots/{cohort}/ctDNAmers_TF_estimates_summary.png",
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
        ctDNAmers_res=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF_ctDNA_det.csv",
    output:
        lead_time=f"{pref}/TF_estimates_analysis/plots/{cohort}/ctDNAmers_lead_time.png",
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
