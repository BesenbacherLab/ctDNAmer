rule run_sequenza_create_gc_wiggle_track:
    input:
        ref_fa=reference_fasta,
    output:
        gc_wiggle_track=get_sequenza_wiggle_track_from_ref_fasta(),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 1440,
    log:
        "logs/sequenza/create_gc_wiggle_track_from_ref_fasta_window50.out",
    conda:
        "../envs/sequenza_utils3_0_0.yaml"
    shell:
        "sequenza-utils gc_wiggle -w 50 --fasta {input.ref_fa} -o {output.gc_wiggle_track}  2> {log}"


rule sequenza_create_seqz_file:
    input:
        normal_bam=get_normal_bam,
        tumor_bam=get_tumor_bam,
        fasta=reference_fasta,
        gc_track=get_sequenza_wiggle_track_from_ref_fasta(),
    output:
        "results/{pt}/sequenza/out.seqz.qz"
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 8640,
    log:
        "logs/patients/{pt}/sequenza/sequenza_create_seqz_file.out"
    conda:
        "../envs/sequenza_utils3_0_0.yaml"
    shell:
        "sequenza-utils bam2seqz -n {input.normal_bam} -t {input.tumor_bam} --fasta {input.fasta} -gc {input.gc_track} -o {output} --qformat illumina  2> {log}"


rule run_sequenza_bin_seqz_file:
    input:
        "results/{pt}/sequenza/out.seqz.qz",
    output:
        "results/{pt}/sequenza/small_w50.seqz.qz",
    resources:
        mem_mb=10000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/sequenza/sequenza_bin_seqz_file.out"
    conda:
        "../envs/sequenza_utils3_0_0.yaml"
    shell:
        "sequenza-utils seqz_binning --seqz {input} -w 50 -o {output}  2> {log}"


rule install_sequenza:
    output:
        installation_check=temp("results/sequenza_R_installation_check.txt"),
    resources:
        mem_mb=20000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/install_sequenza.out"
    conda:
        "../envs/sequenza_R.yaml"
    script:
        "../scripts/install_sequenza.R"


rule run_sequenza_Rfit:
    input:
        installation_check="results/sequenza_R_installation_check.txt",
        seq_binned="results/{pt}/sequenza/small_w50.seqz.qz",
    output:
        output_segments="results/{pt}/sequenza/{pt}_segments.txt",
        purity_and_ploidy="results/{pt}/sequenza/{pt}_confints_CP.txt",
    resources:
        mem_mb=50000,
        runtime=lambda wildcards, attempt: attempt * 360,
    log:
        "logs/patients/{pt}/sequenza/run_sequenza_Rfit.out"
    params:
        pt_id="{pt}",
        output_dir="results/{pt}/sequenza/",
        female=lambda wildcards: str(samples.female[samples.sample_ID == wildcards.pt].item()),
        cellularity_min=lambda wildcards: str(samples.tp_min[samples.sample_ID == wildcards.pt].item()),
        cellularity_max=lambda wildcards: str(samples.tp_max[samples.sample_ID == wildcards.pt].item()),
        ploidy_min=lambda wildcards: str(samples.pl_min[samples.sample_ID == wildcards.pt].item()),
        ploidy_max=lambda wildcards: str(samples.pl_max[samples.sample_ID == wildcards.pt].item()),
    conda:
        "../envs/sequenza_R.yaml"
    script:
        "../scripts/run_sequenza.R"