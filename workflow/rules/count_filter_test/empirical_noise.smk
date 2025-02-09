
rule count_donor_kmers:
    input:
        unpack(get_donor_input),
    output:
        kmers_pre="results/donors/{donor}/kmers.kmc_pre",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/donors/{donor}/count_kmers.out"
    params:
        k=k(),
        tmpdir="kmc_bins/count_kmers/donors/{donor}/",
        o_suf_rm="results/donors/{donor}/kmers",
        input_format=get_donor_fileformat,
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        'kmc -k{params.k} -m50 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.input_format} @{input.input_files} {params.o_suf_rm} {params.tmpdir} 2> {log}'


rule make_count_summary_donors:
    input:
        kmers="results/donors/{donor}/kmers.kmc_pre",
    output:
        count_sum=temp("results/donors/{donor}/kmers_count_sum.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/donors/{donor}/make_count_summary_donors.out"
    params:
        i_suf_rm="results/donors/{donor}/kmers",
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram {output.count_sum} -ci1 -cx1000000000 2> {log}"


rule filter_count_summary_donors:
    input:
        count_sum="results/donors/{donor}/kmers_count_sum.txt",
    output:
        count_sum_filt="results/donors/{donor}/kmers_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/donors/{donor}/filter_count_summary_donors.out"
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.count_sum} > {output.count_sum_filt} 2> {log}"


rule calculate_donor_mean_count:
    input:
        data="results/donors/{donor}/kmers_count_sum_0rm.txt",
    output:
        mean_count="results/donors/{donor}/mean_count.csv",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/donors/{donor}/calculate_donor_mean_count.out"
    params:
        max_count=1000,
    conda:
        "../../envs/R4_1.yaml"
    script:
        "../../scripts/calculate_cfDNA_mean_count.R"


rule intersect_UT_and_donor_kmers:
    input:
        unique_tumor="results/patients/{pt}/unique_tumor/UT.kmc_pre",
        unique_tumor_suf="results/patients/{pt}/unique_tumor/UT.kmc_suf",
        donor="results/donors/{donor}/kmers.kmc_pre",
    output:
        int_UT_c=temp("results/patients/{pt}/empirical_noise/{donor}/UT_counts.kmc_pre"),
        int_UT_c_suf=temp("results/patients/{pt}/empirical_noise/{donor}/UT_counts.kmc_suf"),
        int_donor_c=temp("results/patients/{pt}/empirical_noise/{donor}/donor_counts.kmc_pre"),
        int_donor_c_suf=temp("results/patients/{pt}/empirical_noise/{donor}/donor_counts.kmc_suf"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/empirical_noise/intersect_UT_and_{donor}_donor_kmers.out",
    params:
        i_UT_suf_rm="results/patients/{pt}/unique_tumor/UT",
        i_donor_suf_rm="results/donors/{donor}/kmers",
        o_suf_rm_UT_c="results/patients/{pt}/empirical_noise/{donor}/UT_counts",
        o_suf_rm_donor_c="results/patients/{pt}/empirical_noise/{donor}/donor_counts",
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools -t{threads} -v simple {params.i_UT_suf_rm} -ci1 -cx100000 {params.i_donor_suf_rm} -ci1 -cx1000000000 intersect {params.o_suf_rm_UT_c} -ci1 -cs1000000000 -cx1000000000 -ocleft intersect {params.o_suf_rm_donor_c} -ci1 -cs1000000000 -cx1000000000 -ocright 2> {log}" 


rule dump_intersection_UT_counts:
    input:
        kmers="results/patients/{pt}/empirical_noise/{donor}/UT_counts.kmc_pre",
        kmers_suf="results/patients/{pt}/empirical_noise/{donor}/UT_counts.kmc_suf",
    output:
        dump=temp("results/patients/{pt}/empirical_noise/{donor}/UT_counts.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/empirical_noise/dump_{donor}_intersection_UT_counts_.out"
    params:
        i_suf_rm="results/patients/{pt}/empirical_noise/{donor}/UT_counts",
        c_lower="1",
        c_higher="1000000000",
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci{params.c_lower} -cx{params.c_higher} dump {output.dump} -ci{params.c_lower} -cx{params.c_higher} -cs1000000000 2> {log}"


rule dump_intersection_donor_counts:
    input:
        kmers="results/patients/{pt}/empirical_noise/{donor}/donor_counts.kmc_pre",
        kmers_suf="results/patients/{pt}/empirical_noise/{donor}/donor_counts.kmc_suf",
    output:
        dump=temp("results/patients/{pt}/empirical_noise/{donor}/donor_counts.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/empirical_noise/dump_{donor}_intersection_donor_counts.out"
    params:
        i_suf_rm="results/patients/{pt}/empirical_noise/{donor}/donor_counts",
        c_lower="1",
        c_higher="1000000000",
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci{params.c_lower} -cx{params.c_higher} dump {output.dump} -ci{params.c_lower} -cx{params.c_higher} -cs1000000000 2> {log}"


rule combine_intersection_counts:
    input:
        UT_counts="results/patients/{pt}/empirical_noise/{donor}/UT_counts.txt",
        donor_counts="results/patients/{pt}/empirical_noise/{donor}/donor_counts.txt",
    output:
        combined=temp("results/patients/{pt}/empirical_noise/{donor}/combined_intersection.txt"),
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/empirical_noise/combine_intersection_counts_{donor}.out"
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "join -a1 -e0 -j1 -o auto {input.UT_counts} {input.donor_counts} > {output.combined} 2> {log}"


rule create_intersection_tables:
    input:
        kmers="results/patients/{pt}/empirical_noise/{donor}/combined_intersection.txt",
    output:
        count_table="results/patients/{pt}/empirical_noise/{donor}/combined_int_gc_content_count_table.txt",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/empirical_noise/create_intersection_table_{donor}.out"
    params:
        k=k(),
        max_count=config["count_UB"],
    conda:
        "../../envs/py3_12.yaml"
    script:
        "../../scripts/gc_content_count_table_noise_kmers.py"


rule model_empirical_noise:
    input:
        UT_mdata="results/patients/{pt}/unique_tumor_count_filter_test/UT_mdata_minct{ct_cutoff}.txt",
        noise_count_tables=expand("results/patients/{{pt}}/empirical_noise/{donor}/combined_int_gc_content_count_table.txt", donor=donors.index.tolist()),
        donor_means=expand("results/donors/{donor}/mean_count.csv", donor=donors.index.tolist()),
        model=os.path.join(workflow.basedir, "scripts/models/empirical_noise.stan"),
    output:
        param_estimates="results/patients/{pt}/unique_tumor_count_filter_test/empirical_noise/minct{ct_cutoff}/estimates.csv",
        traceplot="results/patients/{pt}/unique_tumor_count_filter_test/empirical_noise/minct{ct_cutoff}/traceplot.png",
        autocorr="results/patients/{pt}/unique_tumor_count_filter_test/empirical_noise/minct{ct_cutoff}/autocorrelation.png",
        summary_txt="results/patients/{pt}/unique_tumor_count_filter_test/empirical_noise/minct{ct_cutoff}/modeling_summary.txt",
        param_density="results/patients/{pt}/unique_tumor_count_filter_test/empirical_noise/minct{ct_cutoff}/param_density.png",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 360,
    threads: 4
    log:
        "logs/patients/{pt}/unique_tumor_count_filter_test/model_empirical_noise/minct{ct_cutoff}.out"
    params:
        donor_list=donors.index.tolist(),
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
    conda:
        "../../envs/R4_1.yaml"
    script:
        "../../scripts/model_empirical_noise_stan.R"