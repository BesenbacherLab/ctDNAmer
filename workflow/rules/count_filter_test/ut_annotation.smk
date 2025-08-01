rule intersect_cfDNA_UT:
    input:
        kmers1="results/patients/{pt}/{cfDNA_ID}/kmers.kmc_pre",
        kmers2=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT.kmc_pre"),
        kmers2_suf=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT.kmc_suf"),
    output:
        intersection=temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_pre"),
        intersection_suf=temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_suf"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/intersect_cfDNA_UT_ci2.out",
    params:
        input1_suf_rm="results/patients/{pt}/{cfDNA_ID}/kmers",
        input2_suf_rm=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT"),
        output_suf_rm="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection",
        count_lower1="1",
        count_lower2="2",
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools -t{threads} -v simple {params.input1_suf_rm} -ci{params.count_lower1} -cx1000000000 {params.input2_suf_rm} -ci{params.count_lower2} -cx1000000000 intersect {params.output_suf_rm} -ci1 -cs1000000000 -cx1000000000 -ocleft 2> {log}"


rule dump_cfDNA_UT_intersection:
    input:
        kmers="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_pre",
        kmers_suf="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_suf",
    output:
        dump="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.txt",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/dump_cfDNA_and_UT_UTci2_intersection.out",
    params:
        input_suf_rm="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection",
        count_lower="1",
        count_higher="1000000000",
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.input_suf_rm} -ci{params.count_lower} -cx{params.count_higher} dump {output.dump} -ci{params.count_lower} -cx{params.count_higher} -cs1000000000 2> {log}"


rule annotate_UT_set_with_cfDNA_counts:
    input:
        UT_kmers="results/patients/{pt}/unique_tumor_count_filter_test/UT_kmers_filtered_sorted_minct{ct_cutoff}.txt",
        cfDNA_kmers="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.txt",
    output:
        combined="results/patients/{pt}/{cfDNA_ID}/unique_tumor_count_filter_test/minct{ct_cutoff}/UT_cfDNA_annotation.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/unique_tumor_count_filter_test/annotate_UT_set_with_cfDNA_counts/minct{ct_cutoff}.out"
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "join -a1 -e0 -j1 -o auto {input.UT_kmers} {input.cfDNA_kmers} > {output.combined} 2> {log}"
