rule quantify_run_oarfish:
    input:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.bam",
    output:
        "results/quantify/run_oarfish/{sample}.v{version}.quant",
    params:
        output_path="results/quantify/run_oarfish/{sample}.v{version}",
        filter_group=config["oarfish_filter_group"],
    threads: config["oarfish_threads"]
    resources:
        mem_mb=config["quantify_mem_mb"],
    log:
        "logs/quantify/run_oarfish/{sample}.v{version}.out",
    conda:
        "../envs/standalone/oarfish.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        conda list > {log};
        oarfish --threads {threads} \
                --filter-group {params.filter_group} \
                --model-coverage \
                --alignments {input} \
                --output {params.output_path} &>> {log}
        """


rule quantify_all_oarfish:
    input:
        expand(
            "results/quantify/run_oarfish/{sample}.v{version}.quant",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),


rule quantify_run_oarfish_bootstrapped:
    input:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.bam",
    output:
        quant="results/quantify/run_oarfish_bootstrapped/{sample}.v{version}.quant",
        infreps="results/quantify/run_oarfish_bootstrapped/{sample}.v{version}.infreps.pq",
    params:
        output_path="results/quantify/run_oarfish_bootstrapped/{sample}.v{version}",
        filter_group=config["oarfish_filter_group"],
        num_bootstraps=config["num_bootstraps"],
    threads: config["oarfish_threads"]
    resources:
        mem_mb=config["quantify_mem_mb"],
    log:
        "logs/quantify/run_oarfish_bootstrapped/{sample}.v{version}.out",
    conda:
        "../envs/standalone/oarfish.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        conda list > {log};
        oarfish --threads {threads} \
                --filter-group {params.filter_group} \
                --model-coverage \
                --num-bootstraps {params.num_bootstraps} \
                --alignments {input} \
                --output {params.output_path} &>> {log}
        """


rule quantify_all_oarfish_bootstrapped:
    input:
        expand(
            "results/quantify/run_oarfish_bootstrapped/{sample}.v{version}.infreps.pq",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),


rule create_salmon_index:
    input:
        transcriptome="results/transcriptome/gencode.v{version}.transcriptome.fa",
    output:
        directory("results/salmon_index/v{version}"),
    threads: config["quantify_isoforms_threads"]
    log:
        "logs/salmon_index/v{version}.log",
    conda:
        "../envs/standalone/salmon.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        salmon index -t {input.transcriptome} -i {output} -p {threads} &> {log}
        """


rule quantify_run_salmon_illumina:
    input:
        index="results/salmon_index/v{version}",
        first_reads=lambda wc: illumina_reads(wc, 1),
        second_reads=lambda wc: illumina_reads(wc, 2),
    output:
        "results/quantify/run_salmon_illumina/v{version}/{sample}/quant.sf",
    params:
        output_name="results/quantify/run_salmon_illumina/v{version}/{sample}",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
        num_bootstraps=config["num_bootstraps"],
    log:
        "logs/quantify/run_salmon_illumina/v{version}/{sample}.out",
    threads: config["quantify_isoforms_threads"]
    conda:
        "../envs/standalone/salmon.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        conda list > {log};
        salmon quant -i {input.index} -l A \
            -1 {input.first_reads} -2 {input.second_reads} \
            --validateMappings -o {params.output_name} -p {threads} \
            --seqBias --gcBias --fldMean {params.fragment_length_mean} \
            --fldSD {params.fragment_length_sd} \
            --numBootstraps {params.num_bootstraps} &>> {log}
        """


rule quantify_all_salmon_illumina:
    input:
        expand(
            "results/quantify/run_salmon_illumina/v{version}/{sample}/quant.sf",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),


rule create_kallisto_index:
    input:
        transcriptome="results/transcriptome/gencode.v{version}.transcriptome.fa",
    output:
        "results/kallisto_index/v{version}/index.idx",
    log:
        "logs/kallisto_index/v{version}.log",
    conda:
        "../envs/standalone/kallisto.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        kallisto index -i {output} {input.transcriptome} &> {log}
        """
