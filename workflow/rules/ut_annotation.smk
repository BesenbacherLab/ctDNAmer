rule intersect_cfDNA_UT:
    input:
        kmers1="results/patients/{pt}/{cfDNA_ID}/kmers.kmc_pre",
        kmers2=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT.kmc_pre"),
        kmers2_suf=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT.kmc_suf"),
    output:
        int_sec=temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_pre"),
        int_sec_suf=temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_suf"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/intersect_cfDNA_UT_ci2.out"
    params:
        i1_suf_rm="results/patients/{pt}/{cfDNA_ID}/kmers",
        i2_suf_rm=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT"),
        o_suf_rm="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection",
        c_lower1="1",
        c_lower2="2",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers

        kmc_tools -t{threads} -v simple {params.i1_suf_rm} -ci{params.c_lower1} -cx1000000000 {params.i2_suf_rm} -ci{params.c_lower2} -cx1000000000 intersect ${{scratch_kmers}} -ci1 -cs1000000000 -cx1000000000 -ocleft 2> {log}
        
        mv ${{scratch_kmers}}.kmc_pre {params.o_suf_rm}.kmc_pre
        mv ${{scratch_kmers}}.kmc_suf {params.o_suf_rm}.kmc_suf
        '''


rule dump_cfDNA_UT_intersection:
    input:
        kmers="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_pre",
        kmers_suf="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.kmc_suf",
    output:
        dump=temp("results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.txt"),
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 60 * 12,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/dump_cfDNA_and_UT_UTci2_intersection.out"
    params:
        i_suf_rm="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection",
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        '''
        set -e
        tmp_dir=$(mktemp -d --tmpdir=/scratch/$SLURM_JOBID)
        scratch_kmers=/scratch/$SLURM_JOBID/kmers.txt

        kmc_tools transform {params.i_suf_rm} -ci1 -cx1000000000 dump ${{scratch_kmers}} -ci1 -cx1000000000 -cs1000000000 2> {log}
        
        mv ${{scratch_kmers}} {output.dump}
        '''


rule annotate_UT_set_with_cfDNA_counts:
    input:
        UT_kmers=lambda wildcards: str("results/patients/" + str(wildcards.pt) + "/unique_tumor/UT_filtered_sorted.txt"),
        cfDNA_kmers="results/patients/{pt}/{cfDNA_ID}/cfDNA_UT_UTci2_intersection.txt",
    output:
        combined="results/patients/{pt}/{cfDNA_ID}/UT_cfDNA_annotation.txt",
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/patients/{pt}/{cfDNA_ID}/annotate_UT_set_with_cfDNA_counts.out"
    conda:
        "../envs/kmc3_2.yaml"
    shell:
        "join -a1 -e0 -j1 -o auto {input.UT_kmers} {input.cfDNA_kmers} > {output.combined} 2> {log}"
