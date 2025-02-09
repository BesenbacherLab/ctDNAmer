rule install_CNAqc_mobster:
    output:
        installation_check=temp("results/CNAqc_installation_check.txt"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/install_CNAqc_mobster.out"
    conda:
        "../envs/R4_3_CNAqc_mobster.yaml"
    script:
        "../scripts/install_CNAqc_mobster.R"


rule run_CNAqc:
    input:
        installation_check="results/CNAqc_installation_check.txt",
        CNA="results/{pt}/sequenza/{pt}_segments.txt",
        purity_and_ploidy="results/{pt}/sequenza/{pt}_confints_CP.txt",
        SNV_annotated="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.cnaqc_format.vcf",
    output:
        SNVs_for_tracking="results/{pt}/SNVs_QC/QC_passed_low_VAF_rm_SNVs.txt",
        segments="results/{pt}/SNVs_QC/CNAqc/segments.png",
        segments_genome_cov="results/{pt}/SNVs_QC/CNAqc/segments_genome_cov.png",
        data_histograms="results/{pt}/SNVs_QC/CNAqc/data_histograms_VAF_DP_NV.png",
        peaks_QC="results/{pt}/SNVs_QC/CNAqc/peaks_QC.png",
        SNVs_QCPASS_11_VAF="results/{pt}/SNVs_QC/CNAqc/QCPASS_11GT_VAF_filters_plot.png",
        mobster_best_fit="results/{pt}/SNVs_QC/mobster/mobster_best_fit.png",
        mobster_best_5_plot="results/{pt}/SNVs_QC/mobster/mobster_best_5_plot.png",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 720,
    params:
        VAF_min=lambda wildcards: str(samples.VAF_min[samples.sample_ID == wildcards.pt].item()),
        VAF_max=lambda wildcards: str(samples.VAF_max[samples.sample_ID == wildcards.pt].item()),
        mobster_cutoff=lambda wildcards: str(samples.mobster_cutoff[samples.sample_ID == wildcards.pt].item()),
    log:
        "logs/patients/{pt}/CNAqc/run_CNAqc.out"
    conda:
        "../envs/R4_3_CNAqc_mobster.yaml"
    script:
        "../scripts/run_CNAqc_and_mobster.R"