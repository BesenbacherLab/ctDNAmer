

rule plot_emp_noise_data: # SF1, SF4B
    input:
        samples=config["samples"],
        donors=config["donors"],
        units=config["samples_cfDNA"],
        clin_data=config["clinical_data"],
    output:
        emp_noise_data=f"{pref}/TF_estimates_analysis/plots/{cohort}/emp_noise_data.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/plot_emp_noise_data.out"
    params: 
        res_dir=config["ctDNA_mers_dir"],
        cohort=config["cohort"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/plot_emp_noise_data.R"


rule plot_emp_noise_estimates: # SF2, SF4A
    input:
        samples=config["samples"],
    output:
        emp_noise_estimates=f"{pref}/TF_estimates_analysis/plots/{cohort}/emp_noise_estimates.png",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/TF_estimates_analysis/plot_emp_noise_estimates.out"
    params: 
        res_dir=config["ctDNA_mers_dir"],
        cohort=config["cohort"],
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/plot_emp_noise_estimates.R"

