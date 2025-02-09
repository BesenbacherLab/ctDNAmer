rule glu_count_kmers:
    input:
        unpack(get_germline_union_input),
    output:
        kmers="results/germline_union/{glu_type}/{glu_sample}/kmers.kmc_pre",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/glu/count_kmers/{glu_type}/{glu_sample}.out"
    params:
        tmpdir="kmc_bins/count_kmers/germline_union/{glu_type}/{glu_sample}/",
        o_suf_rm="results/germline_union/{glu_type}/{glu_sample}/kmers",
        k=k(),
        i_format=get_glu_fileformat,
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        'kmc -k{params.k} -m50 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.i_format} @{input.input_files} {params.o_suf_rm} {params.tmpdir} 2> {log}'


rule glu_source_make_count_summary:
    input:
        kmers="results/germline_union/{glu_type}/{glu_sample}/kmers.kmc_pre",
    output:
        c_sum=temp("results/germline_union/{glu_type}/{glu_sample}/kmers_count_sum.txt"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 120,
    log:
        "logs/glu/make_count_summary/{glu_type}/{glu_sample}.out"
    params:
        i_suf_rm="results/germline_union/{glu_type}/{glu_sample}/kmers",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram {output.c_sum} -ci1 -cx1000000000 2> {log}"


rule glu_source_filter_count_summary:
    input:
        c_sum="results/germline_union/{glu_type}/{glu_sample}/kmers_count_sum.txt",
    output:
        c_sum_filt="results/germline_union/{glu_type}/{glu_sample}/kmers_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 60,
    log:
        "logs/glu/filter_count_summary/{glu_type}/{glu_sample}.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"


rule loop_union_germline:
    input:
        get_glu_input,
    output:
        combined_file=temp("results/germline_union/germline_union_{glu_samples_comb}.kmc_pre"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 720,
    threads: 4
    log:
        "logs/glu/create_union/{glu_samples_comb}.out"
    params:
        i_suf_rm=recurse_glu_samples,
        o_suf_rm="results/germline_union/germline_union_{glu_samples_comb}",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "kmc_tools -t{threads} -v simple {params.i_suf_rm[1]} -ci1 -cx1000000000 {params.i_suf_rm[0]} -ci1 -cx1000000000 union {params.o_suf_rm} -ci1 -cx1000000000 -cs1000000000 -ocmin 2> {log}"


rule glu_make_count_summary:
    input:
        kmers="results/germline_union/germline_union_{glu_samples_comb}.kmc_pre",
    output:
        c_sum=temp("results/germline_union/germline_union_{glu_samples_comb}_count_sum.txt"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/glu/make_glu_count_summary_{glu_samples_comb}.out"
    params:
        i_suf_rm="results/germline_union/germline_union_{glu_samples_comb}",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram {output.c_sum} -ci1 -cx1000000000 2> {log}"


rule glu_filter_count_summary:
    input:
        c_sum="results/germline_union/germline_union_{glu_samples_comb}_count_sum.txt",
    output:
        c_sum_filt="results/germline_union/germline_union_{glu_samples_comb}_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/glu/filter_glu_count_summary_{glu_samples_comb}.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"


rule glu_union_last:
    input:
        get_glu_input_last,
    output:
        combined_file="results/germline_union/final_germline_union_{glu_samples_comb_last}.kmc_pre",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 720,
    threads: 4
    log:
        "logs/glu/create_union/{glu_samples_comb_last}.out"
    params:
        i_suf_rm=get_final_glu_prefix,
        o_suf_rm="results/germline_union/final_germline_union_{glu_samples_comb_last}",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "kmc_tools -t{threads} -v simple {params.i_suf_rm[1]} -ci1 -cx1000000000 {params.i_suf_rm[0]} -ci1 -cx1000000000 union {params.o_suf_rm} -ci1 -cx1000000000 -cs1000000000 -ocmin 2> {log}"


rule final_glu_make_count_summary:
    input:
        kmers="results/germline_union/final_germline_union_{glu_samples_comb_last}.kmc_pre",
    output:
        c_sum=temp("results/germline_union/final_germline_union_{glu_samples_comb_last}_count_sum.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/glu/final_glu_make_count_summary_{glu_samples_comb_last}.out"
    params:
        i_suf_rm="results/germline_union/final_germline_union_{glu_samples_comb_last}",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram {output.c_sum} -ci1 -cx1000000000 2> {log}"


rule final_glu_filter_count_summary:
    input:
        c_sum="results/germline_union/final_germline_union_{glu_samples_comb_last}_count_sum.txt",
    output:
        c_sum_filt="results/germline_union/final_germline_union_{glu_samples_comb_last}_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/glu/final_glu_filter_count_summary_{glu_samples_comb_last}.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"