rule filter_SNVs:
    input:
        unpack(get_mutect_input)
    output:
        filtered_mutect_calls="results/{pt}/SNVs_prep/mutect_calls_filtered.recode.vcf",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/SNVs_prep/filter_SNVs.out",
    params:
        out_pref="results/{pt}/SNVs_prep/mutect_calls_filtered",
    conda: 
        "../envs/vcftools0_1_16.yaml"
    shell:
        "vcftools --gzvcf {input.m_calls} --out {params.out_pref} --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chr10 --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 --chr chr21 --chr chr22 --chr chrX --chr chrY --remove-indels --remove-filtered-all --recode --recode-INFO-all  2> {log}"


rule run_vep:
    input:
        ref_fasta=reference_fasta,
        m_calls="results/{pt}/SNVs_prep/mutect_calls_filtered.recode.vcf",
    output:
        m_calls_annot="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.vcf",
        vep_stats="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.txt_summary.html",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/SNVs_prep/run_vep.out",
    params:
        cache_dir=vep_cache
    threads: 4
    conda:
        "../envs/vep105_0.yaml"
    shell:
        '''vep --cache -i {input.m_calls} --format vcf -o {output.m_calls_annot} --fork {threads} --cache --offline --verbose --dir_cache {params.cache_dir} --merged --stats_file {output.vep_stats} --fasta {input.ref_fasta} --vcf --show_ref_allele --pick_allele --symbol --fields "Allele,Gene,Feature,Feature_type,Consequence,IMPACT,DISTANCE,STRAND,SYMBOL,SYMBOL_SOURCE" 2> {log}'''


rule run_vep_coding_only:
    input:
        ref_fasta=reference_fasta,
        m_calls="results/{pt}/SNVs_prep/mutect_calls_filtered.recode.vcf"
    output: 
        m_calls_annot="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.coding_only.vcf",
        vep_stats="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.coding_only.txt_summary.html"
    resources: 
        mem_mb=50000, 
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/SNVs_prep/run_vep_coding_only.out", 
    params:
        cache_dir=vep_cache
    threads: 4
    conda: 
        "../envs/vep105_0.yaml"
    shell:
        '''vep --cache -i {input.m_calls} --format vcf -o {output.m_calls_annot} --fork {threads} --cache --offline --verbose --dir_cache {params.cache_dir} --merged --stats_file {output.vep_stats} --fasta {input.ref_fasta} --vcf --show_ref_allele --pick_allele --coding_only --symbol --fields "Allele,Gene,Feature,Feature_type,Consequence,IMPACT,DISTANCE,STRAND,SYMBOL,SYMBOL_SOURCE" 2> {log}'''


rule prep_annotated_SNVs_for_cnaqc:
    input:
        m_calls_annot="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.vcf",
        m_calls_annot_coding="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.coding_only.vcf",
    output:
        m_calls_cnaqc="results/{pt}/SNVs_prep/mutect_calls_filtered.cnaqc_format.vcf",
        m_calls_annot_cnaqc="results/{pt}/SNVs_prep/mutect_calls_filtered.annotated.cnaqc_format.vcf",
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 180,
    log:
        "logs/patients/{pt}/SNVs_prep/prep_annotated_SNVs_for_cnaqc.out",
    conda:
        "../envs/R4_1.yaml"
    script:
        "../scripts/prepare_annotated_SNVs_for_cnaqc.R"