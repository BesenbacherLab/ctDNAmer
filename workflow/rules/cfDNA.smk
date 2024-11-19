rule count_cfdna_kmers:
    input:
        unpack(get_cfDNA_input),
    output:
        kmers = "results/patients/{pt}/{cfDNA_ID}/kmers.kmc_pre",
        kmers_suf = "results/patients/{pt}/{cfDNA_ID}/kmers.kmc_suf",
    resources:
        mem_mb=100000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/count_kmers.out",
    params:
        k = k(),
        tmpdir = "kmc_bins/count_kmers/patients/{pt}/{cfDNA_ID}/",
        output_suf_rm = "results/patients/{pt}/{cfDNA_ID}/kmers",
        input_format = get_cfDNA_fileformat,
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        'kmc -k{params.k} -m100 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.input_format} @{input.input_files} {params.output_suf_rm} {params.tmpdir} 2> {log}'


rule intersect_cfDNA_and_iGL_kmers:
    input:
        kmers1 = "results/patients/{pt}/{cfDNA_ID}/kmers.kmc_pre", 
        kmers1_suf = "results/patients/{pt}/{cfDNA_ID}/kmers.kmc_suf",
        kmers2 = lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/germline/kmers.kmc_pre")
    output:
        intersection = temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.kmc_pre"),
        intersection_suf = temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.kmc_suf")
    resources: 
        mem_mb = 20000, 
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/intersect_cfDNA_and_iGL_kmers.out",
    params:
        input1_suf_rm = "results/patients/{pt}/{cfDNA_ID}/kmers", 
        input2_suf_rm = lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/germline/kmers"),
        output_suf_rm = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc",
        count_lower1 = "1",
        count_lower2 = "1",
    conda: 
        "../envs/kmc3_2.yaml"
    shell:
        "kmc_tools -t{threads} -v simple {params.input1_suf_rm} -ci{params.count_lower1} -cx1000000000 {params.input2_suf_rm} -ci{params.count_lower2} -cx1000000000 intersect {params.output_suf_rm} -ci1 -cs1000000000 -cx1000000000 -ocleft 2> {log}" 


rule make_cfDNA_iGL_intersection_count_summary:
    input:
        intersection = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.kmc_pre",
        intersection_suf = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.kmc_suf",
    output:
        count_sum = temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc_sum.txt")
    resources: 
        mem_mb = 20000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/make_cfDNA_iGL_intersection_count_summary.out", 
    params:
        input_suf_rm = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc"
    conda:
        "../envs/kmc3_2.yaml"
    shell: 
        "kmc_tools transform {params.input_suf_rm} -ci1 -cx1000000000 histogram {output.count_sum} -ci1 -cx1000000000 2> {log}"


rule filter_cfDNA_iGL_intersection_count_summary:
    input:
        count_sum = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc_sum.txt"
    output:
        count_sum_filt = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc_sum_0rm.txt"
    resources: 
        mem_mb = 10000, 
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/make_cfDNA_iGL_intersection_count_summary.out", 
    conda:
        "../envs/kmc3_2.yaml"
    shell: 
        "awk -F '\t' '$2!=0' {input.count_sum} > {output.count_sum_filt} 2> {log}"


rule calculate_cfDNA_iGL_intersection_mean:
    input:
        data = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc_sum_0rm.txt"
    output:
        cfDNA_mean = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNc_mean.csv"
    resources: 
        mem_mb = 1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/calculate_cfDNA_iGL_intersection_mean.out",
    params:
        max_count = 1000
    conda:
        "../envs/R4_1.yaml"
    script: 
        "../scripts/calculate_cfDNA_mean_count.R"


rule dump_cfDNA_iGL_intersection:
    input:
        kmc_pre = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.kmc_pre",
        kmc_suf = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.kmc_suf",
    output:
        dump = temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.txt")
    resources: 
        mem_mb = 10000, 
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/dump_cfDNA_iGL_intersection.out",
    params:
        input_suf_rm = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc",
        count_lower = "2",
        count_higher = "100",
    conda: 
        "../envs/kmc3_2.yaml"
    shell: 
        "kmc_tools transform {params.input_suf_rm} -ci{params.count_lower} -cx{params.count_higher} dump {output.dump} -ci{params.count_lower} -cx{params.count_higher} -cs1000000000 2> {log}"


rule create_cfDNA_indGL_intersection_GC_content_count_tables:
    input:
        kmers = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc.txt"
    output:
        count_table = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_gc_content_counts.txt"
    resources: 
        mem_mb = 10000, 
        runtime=lambda wildcards, attempt: attempt * 1440,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/create_cfDNA_indGL_intersection_GC_content_count_tables.out",
    params:
        max_count = 100,
        k = k(),
    conda:
        "../envs/R4_1.yaml"
    script: 
        "../scripts/gc_content_count_table.py"


rule calculate_cfDNA_iGL_intersection_mean_gc_strat:
    input:
        data = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_gc_content_counts.txt"
    output:
        cfDNA_mean_gc_strat = "results/patients/{pt}/{cfDNA_ID}/cfDNA_iGL_int_cfDNAc_mean_gc_strat.csv"
    resources: 
        mem_mb = 10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/calculate_cfDNA_iGL_intersection_mean_gc_strat.out",
    params:
        max_count = 1000,
        max_cutoff = 100,
    conda:
        "../envs/R4_1.yaml"
    script: 
        "../scripts/calculate_cfDNA_mean_count_gc_strat.R"
