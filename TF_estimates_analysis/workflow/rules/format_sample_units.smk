

rule format_sample_units:
    input:
        samples=config["samples"],
        units=config["samples_cfDNA"],
        clin_data=config["clinical_data"],
        interventions=config["interventions"],
    output:
        formatted_units=f"{pref}/TF_estimates_analysis/metadata/formatted_units.csv",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/format_sample_units.out"
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/format_sample_units.R"


rule aggregate_ctDNAmers_TF_estimates:
    input:
        units=f"{pref}/TF_estimates_analysis/metadata/formatted_units.csv",
    output:
        results=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF.csv",
        ex_cfDNA=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/excluded_units.csv",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/aggregate_ctDNAmers_TF_estimates.out"
    params: 
        cohort = config["cohort"],
        TF_subfd = config["TF_subfd"],
        res_subfd = config["res_subfd"],
        ctDNA_mers_dir = config["ctDNA_mers_dir"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/aggregate_ctDNAmers_estimates.R"


rule find_detection_cutpoint_ctDNAmers_TF:
    input:
        units=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF.csv",
    output:
        results=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF_ctDNA_det.csv",
        results_noTF_zero=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF_ctDNA_det_no_zero_TF.csv",
        detection_stats=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/detection_stats.csv",
        conf_matrix=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/conf_matrix.csv",
        TF_ROC=f"{pref}/TF_estimates_analysis/plots/{cohort}/ctDNAmers_TF_ROC.png",
        TF_barplot=f"{pref}/TF_estimates_analysis/plots/{cohort}/ctDNAmers_TF_barplot.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/find_detection_cutpoint_ctDNAmers_TF.out"
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/ctDNA_detection_labels.R"

if cohort == "frfr":
    rule aggregate_clonalSNVs_AF_estimates:
        input:
            units=f"{pref}/TF_estimates_analysis/metadata/formatted_units.csv",
        output:
            results=f"{pref}/TF_estimates_analysis/clonal_SNVs/units_clonalSNVsmeanAF.csv",
        resources:
            mem_mb=1000,
            runtime=lambda wildcards, attempt: attempt * 180,
        log:
            "logs/TF_estimates_analysis/aggregate_clonalSNVs_AF_estimates.out"
        params: 
            cohort = config["cohort"],
            clonalSNVs_dir = config["clonalSNVs_dir"],
        conda:
            "../envs/R4_1.yaml"
        script:
            "../scripts/aggregate_clonalSNVs_AF_estimates.R"



    rule find_detection_cutpoint_clonalSNVs:
        input:
            units=f"{pref}/TF_estimates_analysis/clonal_SNVs/units_clonalSNVsmeanAF.csv",
        output:
            results=f"{pref}/TF_estimates_analysis/clonal_SNVs/units_clonal_SNVs_ctDNA_det.csv",
            results_noTF_zero=f"{pref}/TF_estimates_analysis/clonal_SNVs/units_clonal_SNVs_ctDNA_det_no_zero_TF.csv",
            detection_stats=f"{pref}/TF_estimates_analysis/clonal_SNVs/detection_stats.csv",
            conf_matrix=f"{pref}/TF_estimates_analysis/clonal_SNVs/conf_matrix.csv",
            TF_ROC=f"{pref}/TF_estimates_analysis/plots/clonal_SNVs_TF_ROC.png",
            TF_barplot=f"{pref}/TF_estimates_analysis/plots/clonal_SNVs_TF_barplot.png",
        resources:
            mem_mb=1000,
            runtime=lambda wildcards, attempt: attempt * 180,
        log:
            "logs/TF_estimates_analysis/find_detection_cutpoint_clonalSNVs.out"
        conda:
            "../envs/R4_1.yaml"
        script:
            "../scripts/ctDNA_detection_labels.R"


