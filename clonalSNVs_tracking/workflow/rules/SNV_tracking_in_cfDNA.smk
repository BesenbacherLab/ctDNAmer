rule filter_clonal_SNVs_based_on_primary_tumor_quality:
    input:
        bam_file=get_tumor_bam,
        snvs_to_track="results/{pt}/SNVs_QC/QC_passed_low_VAF_rm_SNVs.txt",
        reference_fasta=reference_fasta,
    output:
        outfile="results/{pt}/SNVs_QC/QC_passed_low_VAF_rm_SNVs_Tfiltered.txt",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 4320,
    log:
        "logs/patients/{pt}/clonal_mutation_tracking/filter_clonal_SNVs_based_on_primary_tumor_quality.out"
    conda:
        "../envs/py3_12.yaml"
    script:
        "../scripts/primary_tumor_quality_filter_clonal_SNVs.py"


rule track_SNVs_in_cfDNA:
    input:
        bam_file=get_cfDNA_bam,
        snvs_to_track="results/{pt}/SNVs_QC/QC_passed_low_VAF_rm_SNVs_Tfiltered.txt",
        reference_fasta=reference_fasta,
    output:
        outfile="results/{pt}/{cfDNA_ID}/SNVs_in_cfDNA.txt",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 2880,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/track_snvs_in_cfDNA.out"
    conda:
        "../envs/py3_12.yaml"
    script:
        "../scripts/track_clonal_SNVs_in_cfDNA.py"


rule tracked_SNVs_analysis_and_plotting:
    input:
        unpack(get_tracked_SNVs_input)
    output:
        mean_VAF="results/{pt}/clonal_SNVs_in_cfDNA/mean_AF.tsv",
        mean_VAF_timeline="results/{pt}/clonal_SNVs_in_cfDNA/timeline_mean_AF.png",
        mean_VAF_dist="results/{pt}/clonal_SNVs_in_cfDNA/timeline_AF_distributions.png",
    resources:
        mem_mb=5000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tracked_SNVs_analysis_and_plotting.out"
    params:
        pt_id="{pt}",
        cfDNA_timepoints = get_cfDNA_timepoints,
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/calculate_and_plot_mean_VAF.R"