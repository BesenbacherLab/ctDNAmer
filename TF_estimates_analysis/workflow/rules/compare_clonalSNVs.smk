
if cohort == "frfr":
    rule clonalSNV_and_ctDNAmers_comp_ROC_and_correlation:
        input:
            ctDNA_mers_TF=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF_ctDNA_det_no_zero_TF.csv",
            ctDNA_mers_det_stats=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/detection_stats.csv",
            clonal_SNV_AF=f"{pref}/TF_estimates_analysis/clonal_SNVs/units_clonal_SNVs_ctDNA_det_no_zero_TF.csv",
            clonal_SNV_det_stats=f"{pref}/TF_estimates_analysis/clonal_SNVs/detection_stats.csv",
        output:
            ROC_curve=f"{pref}/TF_estimates_analysis/plots/clonalSNVs_ctDNAmers_ROC.png",
            correlation_plot=f"{pref}/TF_estimates_analysis/plots/clonalSNVs_ctDNAmers_correlation.png",
        resources:
            mem_mb=1000,
            runtime=lambda wildcards, attempt: attempt * 180,
        log:
            "logs/TF_estimates_analysis/clonalSNV_and_ctDNAmers_comp_ROC_and_correlation.out"
        conda:
            "../envs/R4_1.yaml"
        script:
            "../scripts/plot_roc_and_correlation.R"


    rule clonalSNV_and_ctDNAmers_comp_pt_timelines:
        input:
            ctDNA_mers_TF=f"{pref}/TF_estimates_analysis/ctDNAmers_TF/{cohort}/units_ctDNAmersTF_ctDNA_det_no_zero_TF.csv",
            clonal_SNV_AF=f"{pref}/TF_estimates_analysis/clonal_SNVs/units_clonal_SNVs_ctDNA_det_no_zero_TF.csv",
            imagings = config["imagings"],
            interventions = config["interventions"]
        output:
            recur_pt_timelines=f"{pref}/TF_estimates_analysis/plots/clonalSNVs_ctDNAmers_recur_pt_timelines.png",
            no_recur_pt_timelines=f"{pref}/TF_estimates_analysis/plots/clonalSNVs_ctDNAmers_no_recur_pt_timelines.png",
        resources:
            mem_mb=1000,
            runtime=lambda wildcards, attempt: attempt * 180,
        log:
            "logs/TF_estimates_analysis/clonalSNV_and_ctDNAmers_comp_pt_timelines.out"
        conda:
            "../envs/R4_1.yaml"
        script:
            "../scripts/plot_pt_timelines.R"