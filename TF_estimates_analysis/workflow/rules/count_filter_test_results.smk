

rule plot_count_filter_test_results: # SF6
    input:
        cft_samples=config["samples_count_filter_test"],
        cft_units=config["units_count_filter_test"],
    output:
        cft_res="results/TF_estimates_analysis/plots/count_filter_test_TF_estimates.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/plot_count_filter_test_results.out"
    params: 
        res_dir = config["count_filter_test_results_dir"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/plot_count_filter_test_results.R"