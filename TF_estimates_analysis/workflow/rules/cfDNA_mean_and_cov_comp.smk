

rule plot_cfDNA_mean_and_cov_comparison: # SF7
    input:
        samples=config["samples"],
        units=config["samples_cfDNA"],
        wgs_cov=config["wgs_coverage_mean_and_sd"],
    output:
        cfDNA_mean_vs_cov=f"{pref}/TF_estimates_analysis/plots/{cohort}/cfDNA_kmer_mean_vs_seq_cov_mean.png",
        cfDNA_sd_vs_cov=f"{pref}/TF_estimates_analysis/plots/{cohort}/cfDNA_kmer_sd_vs_seq_cov_sd.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/plot_cfDNA_mean_and_cov_comparison.out"
    params: 
        res_dir=config["ctDNA_mers_dir"],
        cohort=config["cohort"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/plot_cfDNA_mean_vs_cov.R"