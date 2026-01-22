
rule create_fasta_of_unique_tumor_kmers: 
    input: 
        kmers = "results/patients/{pt}/unique_tumor/UT_filtered_sorted.txt",
    output: 
        fasta = "results/UT_kmers_mapping/patients/{pt}/UT_filtered_sorted.fa",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/UT_kmers_mapping/{pt}/create_fasta_of_unique_tumor_kmers.out",
    conda: 
        "../envs/py3_12.yaml"
    script:
        "../scripts/fasta_from_txt.py"

rule map_unique_tumor_kmers_to_reference:
    input:
        fasta = "results/UT_kmers_mapping/patients/{pt}/UT_filtered_sorted.fa",
        reference = "/home/coroperv/MomaReference/BACKUP/hg38/software_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.GRC_exclusions_masked.fna",
    output:
        sam = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.sam",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/UT_kmers_mapping/{pt}/map_unique_tumor_kmers_to_reference.out",
    conda: 
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.reference} {input.fasta} > {output.sam}"

rule create_bam_of_unique_tumor_kmers:
    input:
        sam = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.sam",
    output: 
        bam = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.bam",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/UT_kmers_mapping/{pt}/create_bam_of_unique_tumor_kmers.out",
    conda: 
        "../envs/conda_samtools.yaml"
    shell:
        "samtools view -S -b {input.sam} > {output.bam}"

rule sort_bam_of_unique_tumor_kmers:
    input:
        bam = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.bam",
    output: 
        sorted = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.sorted.bam",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/UT_kmers_mapping/{pt}/sort_bam_of_unique_tumor_kmers.out",
    conda: 
        "../envs/conda_samtools.yaml"
    shell:
        "samtools sort {input.bam} -o {output.sorted}"
    
rule index_bam_of_unique_tumor_kmers:
    input:
        bam = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.sorted.bam",
    output: 
        index = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.sorted.bam.bai",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/UT_kmers_mapping/{pt}/index_bam_of_unique_tumor_kmers.out",
    conda: 
        "../envs/conda_samtools.yaml"
    shell:
        "samtools index {input.bam}"

rule analyse_mapping_char_and_overlap_w_mutect: 
    input:
        unpack(get_mutect_and_delly_input),
        kmers = "results/UT_kmers_mapping/patients/{pt}/UT_filtered.sorted.bam",
    output:
        summary = "results/UT_kmers_mapping/patients/{pt}/UT_filtered_alignment_sum.tsv",
        one_mm_out = "results/UT_kmers_mapping/patients/{pt}/UT_filtered_1mm.txt",
        more_unique_out = "results/UT_kmers_mapping/patients/{pt}/UT_filtered_more_unique.txt",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 60 * 6,
    log:
        "logs/UT_kmers_mapping/{pt}/analyse_mapping_char_and_overlap_w_mutect.out",
    conda: 
        "../envs/py3_12.yaml"
    script:
        "../scripts/analyse_mapping_characteristics.py"

