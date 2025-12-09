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
        o_suf_rm="results/patients/{pt}/germline/kmers",
        i_format=get_germline_fileformat,
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers
        scratch_kmctmp=/scratch/$SLURM_JOBID/temp_kmc/
    
        mkdir ${{scratch_kmctmp}}

        kmc -k{params.k} -m50 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.i_format} @{input.input_files} ${{scratch_kmers}} ${{scratch_kmctmp}} 2> {log}

        mv ${{scratch_kmers}}.kmc_pre {params.o_suf_rm}.kmc_pre
        mv ${{scratch_kmers}}.kmc_suf {params.o_suf_rm}.kmc_suf
        rm -r ${{scratch_kmctmp}}
        '''


rule make_germline_count_summary:
    input:
        kmers="results/patients/{pt}/germline/kmers.kmc_pre",
    output:
        c_sum=temp("results/patients/{pt}/germline/kmers_count_sum.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/germline/make_germline_count_summary.out"
    params:
        i_suf_rm="results/patients/{pt}/germline/kmers",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers.txt

        kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram ${{scratch_kmers}} -ci1 -cx1000000000 2> {log}
        
        mv ${{scratch_kmers}} {output.c_sum}
        '''



rule filter_germline_count_summary:
    input:
        c_sum="results/patients/{pt}/germline/kmers_count_sum.txt",
    output:
        c_sum_filt="results/patients/{pt}/germline/kmers_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/germline/filter_germline_count_summary.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"


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
            "../envs/py3_12.yaml"
        script:
            "../scripts/quality_filter_tumor_bam.py"
else:
    rule quality_filter_tumor_fastq:
        input:
            unpack(get_tumor_input),
        output:
            "results/patients/{pt}/tumor/input_files_filtered.txt",
        resources:
            mem_mb=10000,
            runtime=lambda wildcards, attempt: attempt * 360,
        log:
            "logs/patients/{pt}/tumor/quality_filter_tumor_input_fastq.out"
        conda:
            "../envs/py3_12.yaml"
        script:
            "../scripts/quality_filter_tumor_fastq.sh"


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
        i_format=get_tumor_fileformat,
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers
        scratch_kmctmp=/scratch/$SLURM_JOBID/temp_kmc/
    
        mkdir ${{scratch_kmctmp}}

        kmc -k{params.k} -m50 -t{threads} -ci1 -cs1000000000 -cx1000000000 -f{params.i_format} @{input.input_files} ${{scratch_kmers}} ${{scratch_kmctmp}} 2> {log}

        mv ${{scratch_kmers}}.kmc_pre {params.o_suf_rm}.kmc_pre
        mv ${{scratch_kmers}}.kmc_suf {params.o_suf_rm}.kmc_suf
        rm -r ${{scratch_kmctmp}}
        '''


rule make_tumor_count_summary:
    input:
        kmers="results/patients/{pt}/tumor/kmers.kmc_pre",
    output:
        c_sum=temp("results/patients/{pt}/tumor/kmers_count_sum.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/make_tumor_count_summary.out"
    params:
        i_suf_rm="results/patients/{pt}/tumor/kmers",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers.txt

        kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram ${{scratch_kmers}} -ci1 -cx1000000000 2> {log}
        
        mv ${{scratch_kmers}} {output.c_sum}
        '''


rule filter_tumor_count_summary:
    input:
        c_sum="results/patients/{pt}/tumor/kmers_count_sum.txt",
    output:
        c_sum_filt="results/patients/{pt}/tumor/kmers_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/filter_tumor_count_summary.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"


rule calculate_tumor_overall_mean:
    input:
        data="results/patients/{pt}/tumor/kmers_count_sum_0rm.txt",
    output:
        tumor_signal_mean="results/patients/{pt}/tumor/signal_mean_count.csv",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/calculate_tumor_overall_mean.out"
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/calculate_tumor_signal_mean_count.R"


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
        cx_out=str(config["count_UB"]),
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers

        kmc_tools simple {params.i_t_suf_rm} -ci1 -cx1000000000 {params.i_gl_suf_rm} -ci1 -cx1000000000 kmers_subtract ${{scratch_kmers}} -ci{params.ci_out} -cx{params.cx_out} -cs1000000000 2> {log}
        
        mv ${{scratch_kmers}}.kmc_pre {params.o_suf_rm}.kmc_pre
        mv ${{scratch_kmers}}.kmc_suf {params.o_suf_rm}.kmc_suf
        '''


rule make_tumor_sub_indGL_count_summary:
    input:
        kmers="results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_pre",
        kmers_suf="results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_suf",
    output:
        c_sum=temp("results/patients/{pt}/unique_tumor/tumor_sub_indGL_count_sum.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/make_tumor_sub_indGL_count_summary.out"
    params:
        i_suf_rm="results/patients/{pt}/unique_tumor/tumor_sub_indGL",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers.txt

        kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram ${{scratch_kmers}} -ci1 -cx1000000000 2> {log}
        
        mv ${{scratch_kmers}} {output.c_sum}
        '''


rule filter_tumor_sub_indGL_count_summary:
    input:
        c_sum="results/patients/{pt}/unique_tumor/tumor_sub_indGL_count_sum.txt",
    output:
        c_sum_filt="results/patients/{pt}/unique_tumor/tumor_sub_indGL_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/filter_tumor_sub_indGL_count_summary.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"


rule subtract_germline_union_from_tumor:
    input:
        tumor="results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_pre",
        tumor_suf="results/patients/{pt}/unique_tumor/tumor_sub_indGL.kmc_suf",
        germline=get_final_glu(),
    output:
        diff_pre="results/patients/{pt}/unique_tumor/UT.kmc_pre",
        diff_suf="results/patients/{pt}/unique_tumor/UT.kmc_suf",
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
        cx_out=str(config["count_UB"]),
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers

        kmc_tools simple {params.i_t_suf_rm} -ci1 -cx1000000000 {params.i_gl_suf_rm} -ci1 -cx1000000000 kmers_subtract ${{scratch_kmers}} -ci{params.ci_out} -cx{params.cx_out} -cs1000000000 2> {log}
        
        mv ${{scratch_kmers}}.kmc_pre {params.o_suf_rm}.kmc_pre
        mv ${{scratch_kmers}}.kmc_suf {params.o_suf_rm}.kmc_suf
        '''


rule make_tumor_sub_indGL_and_GLU_count_summary:
    input:
        kmers="results/patients/{pt}/unique_tumor/UT.kmc_pre",
        kmers_suf="results/patients/{pt}/unique_tumor/UT.kmc_suf",
    output:
        c_sum=temp("results/patients/{pt}/unique_tumor/UT_count_sum.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/make_tumor_sub_indGL_and_GLU_count_summary.out"
    params:
        i_suf_rm="results/patients/{pt}/unique_tumor/UT",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers.txt

        kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 histogram ${{scratch_kmers}} -ci1 -cx1000000000 2> {log}
        
        mv ${{scratch_kmers}} {output.c_sum}
        '''


rule filter_tumor_sub_indGL_and_GLU_count_summary:
    input:
        c_sum="results/patients/{pt}/unique_tumor/UT_count_sum.txt",
    output:
        c_sum_filt="results/patients/{pt}/unique_tumor/UT_count_sum_0rm.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/tumor/filter_tumor_sub_indGL_and_GLU_count_summary.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "awk -F '\t' '$2!=0' {input.c_sum} > {output.c_sum_filt} 2> {log}"


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
        c_lower=str(config["count_LB"]),
        c_higher=str(config["count_UB"]),
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers.txt

        kmc_tools transform {params.i_suf_rm} -ci{params.c_lower} -cx{params.c_higher} dump ${{scratch_kmers}} -ci{params.c_lower} -cx{params.c_higher} -cs1000000000 2> {log}
        
        mv ${{scratch_kmers}} {output.dump}
        '''


rule filter_ut_kmers:
    input:
        UT_kmers="results/patients/{pt}/unique_tumor/UT.txt",
        tumor_mean="results/patients/{pt}/tumor/signal_mean_count.csv",
    output:
        filtered=temp("results/patients/{pt}/unique_tumor/UT_filtered.txt"),
        filtered_mdata="results/patients/{pt}/unique_tumor/UT_n_and_ci.txt",
        UT_data_gc_table="results/patients/{pt}/unique_tumor/UT_gc_content_count_summary.txt",
    resources:
        mem_mb=75000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/unique_tumor/filter_ut_kmers.out"
    params:
        gc_lower=config["gc_lower"],
        gc_upper=config["gc_upper"],
        max_Tcount=config["count_UB"],
        min_UT_size=config["min_UT_set_size"],
    conda:
        "../envs/py3_12.yaml"
    script:
        "../scripts/filter_UT_kmers.py"


rule sort_ut_kmers:
    input:
        UT_kmers="results/patients/{pt}/unique_tumor/UT_filtered.txt",
    output:
        UT_sorted="results/patients/{pt}/unique_tumor/UT_filtered_sorted.txt",
    resources:
        mem_mb=5000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/unique_tumor/sort_ut_kmers.out"
    conda:
        "../envs/py3_12.yaml"
    shell:
        "sort -k1 -o {output.UT_sorted} {input.UT_kmers}"
