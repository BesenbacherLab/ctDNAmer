rule count_germline_kmers:
    input:
        unpack(get_germline_input),
    output:
        kmers="results/patients/{pt}/germline/kmers.kmc_pre",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/germline/count_kmers.out"
    params:
        k=k(),
        tmpdir="kmc_bins/count_kmers/patients/{pt}/germline/",
        ou_suf_rm="results/patients/{pt}/germline/kmers",
        input_format=get_germline_fileformat,
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        'kmc -k{params.k} -m50 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.input_format} @{input.input_files} {params.ou_suf_rm} {params.tmpdir} 2> {log}'


if config["tumor_fileformat"].upper() == "BAM":
    rule quality_filter_tumor_bam:
        input:
            unpack(get_tumor_input),
        output:
            filtered_files="results/patients/{pt}/tumor/input_files_filtered.txt",
        resources:
            mem_mb=10000,
            runtime=lambda wildcards, attempt: attempt * 2880,
        log:
            "logs/patients/{pt}/tumor/quality_filter_tumor_input_bam.out"
        conda:
            "../../envs/py3_12.yaml"
        script:
            "../../scripts/quality_filter_tumor_bam.py"
else:
    rule quality_filter_tumor_fastq:
        input:
            unpack(get_tumor_input),
        output:
            filtered_files = "results/patients/{pt}/tumor/input_files_filtered.txt",
        resources:
            mem_mb = 10000,
            runtime=lambda wildcards, attempt: attempt * 360,
        log:
            "logs/patients/{pt}/tumor/quality_filter_tumor_input_fastq.out"
        conda:
            "../../envs/py3_12.yaml"
        script:
            "../../scripts/quality_filter_tumor_fastq.sh"


rule count_tumor_kmers:
    input:
        input_files="results/patients/{pt}/tumor/input_files_filtered.txt",
    output:
        kmers="results/patients/{pt}/tumor/kmers.kmc_pre",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/tumor/count_kmers.out"
    params:
        k=k(),
        tmpdir="kmc_bins/count_kmers/patients/{pt}/tumor/",
        o_suf_rm="results/patients/{pt}/tumor/kmers",
        input_format=get_tumor_fileformat,
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        'kmc -k{params.k} -m50 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.input_format} @{input.input_files} {params.o_suf_rm} {params.tmpdir} 2> {log}'


rule subtract_individual_germline_from_tumor:
    input:
        tumor="results/patients/{pt}/tumor/kmers.kmc_pre",
        germline="results/patients/{pt}/germline/kmers.kmc_pre",
    output:
        diff=temp("results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_pre"),
        diff_suf=temp("results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_suf"),
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 720,
    log:
        "logs/patients/{pt}/unique_tumor/subtract_indGL_from_tumor.out"
    params:
        i_t_suf_rm="results/patients/{pt}/tumor/kmers",
        i_gl_suf_rm="results/patients/{pt}/germline/kmers",
        o_suf_rm="results/patients/{pt}/unique_tumor/tumor_sub_indGL",
        ci_out="1",
        cx_out=config["count_UB"],
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools simple {params.i_t_suf_rm} -ci1 -cx1000000000 {params.i_gl_suf_rm} -ci1 -cx1000000000 kmers_subtract {params.o_suf_rm} -ci{params.ci_out} -cx{params.cx_out} -cs1000000000 2> {log}"


rule subtract_germline_union_from_tumor:
    input:
        tumor="results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_pre",
        tumor_suf="results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_suf",
        germline=get_final_glu(),
    output:
        diff_pre=temp("results/patients/{pt}/unique_tumor/UT.kmc_pre"),
        diff_suf=temp("results/patients/{pt}/unique_tumor/UT.kmc_suf"),
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 720,
    log:
        "logs/patients/{pt}/unique_tumor/subtract_glu_from_tumor.out"
    params:
        i_t_suf_rm="results/patients/{pt}/unique_tumor/tumor_sub_indGL",
        i_gl_suf_rm=get_final_glu().split(".")[0],
        o_suf_rm="results/patients/{pt}/unique_tumor/UT",
        ci_out="1",
        cx_out=config["count_UB"],
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools simple {params.i_t_suf_rm} -ci1 -cx1000000000 {params.i_gl_suf_rm} -ci1 -cx1000000000 kmers_subtract {params.o_suf_rm} -ci{params.ci_out} -cx{params.cx_out} -cs1000000000 2> {log}"


rule dump_ut_kmers:
    input:
        kmers="results/patients/{pt}/unique_tumor/UT.kmc_pre",
        kmers_suf="results/patients/{pt}/unique_tumor/UT.kmc_suf",
    output:
        dump=temp("results/patients/{pt}/unique_tumor/UT.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/unique_tumor/dump_kmers.out"
    params:
        i_suf_rm="results/patients/{pt}/unique_tumor/UT",
        c_lower=config["count_LB"],
        c_higher=config["count_UB"],
    conda:
        "../../envs/kmc3_2.yaml"
    shell:
        "kmc_tools transform {params.i_suf_rm} -ci{params.c_lower} -cx{params.c_higher} dump {output.dump} -ci{params.c_lower} -cx{params.c_higher} -cs1000000000 2> {log}"


checkpoint create_count_filtered_ut_sets:
    input:
        UT_kmers="results/patients/{pt}/unique_tumor/UT.txt",
    output:
        count_filtered_dir=directory("results/patients/{pt}/unique_tumor_count_filter_test/"),
    resources:
        mem_mb=75000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/unique_tumor_count_filter_test/create_count_filtered_ut_sets.out"
    params:
        out_pref="results/patients/{pt}/unique_tumor_count_filter_test/",
        min_Tcount=5,
        max_Tcount=100,
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
        UT_size_upper_limit=500000,
        UT_min_size=500,
    conda:
        "../../envs/py3_12.yaml"
    script:
        "../../scripts/create_count_filtered_ut_sets.py"


rule sort_count_filtered_ut_kmers:
    input:
        UT_kmers="results/patients/{pt}/unique_tumor_count_filter_test/UT_kmers_filtered_minct{ct_cutoff}.txt",
    output:
        sorted="results/patients/{pt}/unique_tumor_count_filter_test/UT_kmers_filtered_sorted_minct{ct_cutoff}.txt",
    resources:
        mem_mb=5000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/unique_tumor_count_filter_test/sort_count_filtered_ut_kmers/minct{ct_cutoff}.out"
    conda:
        "../../envs/py3_12.yaml"
    shell:
        "sort -k1 -o {output.sorted} {input.UT_kmers}"
