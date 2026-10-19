rule metric_run_kvalue:
    input:
        binary="results/tools/miniquant/miniQuant_linux/miniQuant",
        gtf="results/download/transcriptome/gencode.v{version}.{kind}.primary_assembly.annotation.gtf.gz",
    output:
        "results/metric/kvalue/{kind}/v{version}/kvalues.tsv",
    params:
        plain_gtf="results/metric/kvalue/{kind}/v{version}.annotation.gtf",
        output_dir="results/metric/kvalue/{kind}/v{version}",
        threads=config["metric"]["kvalue_threads"],
    log:
        "logs/metric/kvalue/{kind}/v{version}.log",
    wildcard_constraints:
        version=r"\d+",
        kind="comprehensive|basic",
    shell:
        """
        zcat {input.gtf} > {params.plain_gtf}
        {input.binary} kvalue -a {params.plain_gtf} -o {params.output_dir} -t {params.threads} &> {log}
        rm {params.plain_gtf}
        awk -v ver={wildcards.version} 'BEGIN {{OFS="\\t"}} NR==1 {{print $0, "version"}} NR>1 {{print $0, ver}}' {output} > {output}.tmp
        mv {output}.tmp {output}
        """

kvalue = expand(
            "results/metric/kvalue/comprehensive/v{version}/kvalues.tsv",
            version=config["gencode"]["versions"],
        )

kvalue_basic = expand(
            "results/metric/kvalue/basic/v{version}/kvalues.tsv",
            version=config["gencode"]["versions"],
        )


rule metric_samtools_stats_kinnex:
    input:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.sorted.bam",
    output:
        "results/metric/data_quality/kinnex/v{version}/{sample}.stats.txt",
    log:
        "logs/metric/data_quality/kinnex/v{version}_{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["metric"]["samtools_stats_threads"]
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        samtools stats -@ {threads} {input} > {output} 2> {log}
        """


rule metric_samtools_stats_illumina:
    input:
        "results/align/run_illumina_star/v{version}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/metric/data_quality/illumina/v{version}/{sample}.stats.txt",
    log:
        "logs/metric/data_quality/illumina/v{version}_{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["metric"]["samtools_stats_threads"]
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        samtools stats -@ {threads} {input} > {output} 2> {log}
        """


data_quality = expand(
    "results/metric/data_quality/kinnex/v{version}/{sample}.stats.txt",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
) + expand(
    "results/metric/data_quality/illumina/v{version}/{sample}.stats.txt",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)


rule metric_ambiguity_oarfish:
    input:
        quant="results/quantify/run_oarfish/{sample}.v{version}.quant",
        ambig_info="results/quantify/run_oarfish/{sample}.v{version}.ambig_info.tsv",
        gene_map="results/aggregate/oarfish/v{version}/{sample}.transcript_counts.tsv",
    output:
        transcript="results/metric/ambiguity/oarfish/v{version}/{sample}.transcript_ambiguity.tsv",
        gene="results/metric/ambiguity/oarfish/v{version}/{sample}.gene_ambiguity.tsv",
    params:
        cell_line=config["cell_line"],
        transcript_id_col_ix=1,
        unique_col_ix=1,
        ambig_col_ix=2,
    conda:
        "../envs/standalone/pandas.yaml"
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/aggregate_ambiguity.py"


rule metric_ambiguity_salmon:
    input:
        quant="results/quantify/run_salmon_illumina/v{version}/{sample}/quant.sf",
        ambig_info="results/quantify/run_salmon_illumina/v{version}/{sample}/aux_info/ambig_info.tsv",
        gene_map="results/aggregate/salmon/v{version}/{sample}.transcript_counts.tsv",
    output:
        transcript="results/metric/ambiguity/salmon/v{version}/{sample}.transcript_ambiguity.tsv",
        gene="results/metric/ambiguity/salmon/v{version}/{sample}.gene_ambiguity.tsv",
    params:
        cell_line=config["cell_line"],
        transcript_id_col_ix=1,
        unique_col_ix=1,
        ambig_col_ix=2,
    conda:
        "../envs/standalone/pandas.yaml"
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/aggregate_ambiguity.py"


oarfish_gene_ambiguity = expand(
    "results/metric/ambiguity/oarfish/v{version}/{sample}.gene_ambiguity.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

salmon_gene_ambiguity = expand(
    "results/metric/ambiguity/salmon/v{version}/{sample}.gene_ambiguity.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)


rule metric_salmon_infrv:
    input:
        quant_sf="results/quantify/run_salmon_illumina/v{version}/{sample}/quant.sf",
    output:
        "results/metric/infrv/salmon/v{version}/{sample}.infrv.tsv",
    params:
        sample=lambda wc: wc.sample,
    conda:
        "../envs/r_tximport.yaml"
    log:
        "logs/metric/infrv/salmon/v{version}_{sample}.log",
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/compute_salmon_infrv.R"


salmon_infrv = expand(
    "results/metric/infrv/salmon/v{version}/{sample}.infrv.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)


rule metric_oarfish_infrv:
    input:
        quant="results/quantify/run_oarfish_bootstrapped/{sample}.v{version}.quant",
        infreps="results/quantify/run_oarfish_bootstrapped/{sample}.v{version}.infreps.pq",
    output:
        "results/metric/infrv/oarfish/v{version}/{sample}.infrv.tsv",
    params:
        sample=lambda wc: wc.sample,
    conda:
        "../envs/standalone/pandas.yaml"
    log:
        "logs/metric/infrv/oarfish/v{version}_{sample}.log",
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/compute_oarfish_infrv.py"


oarfish_infrv = expand(
    "results/metric/infrv/oarfish/v{version}/{sample}.infrv.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)


rule metric_gene_perplexity:
    input:
        transcript_counts="results/aggregate/{platform}/v{version}/{sample}.transcript_counts.tsv",
    output:
        "results/metric/perplexity/{platform}/v{version}/{sample}.perplexity.tsv",
    params:
        sample=lambda wc: wc.sample,
        platform=lambda wc: wc.platform,
    conda:
        "../envs/standalone/pandas.yaml"
    log:
        "logs/metric/perplexity/{platform}/v{version}_{sample}.log",
    wildcard_constraints:
        version=r"\d+",
        platform="oarfish|salmon",
    script:
        "../scripts/compute_gene_perplexity.py"


oarfish_perplexity = expand(
    "results/metric/perplexity/oarfish/v{version}/{sample}.perplexity.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

salmon_perplexity = expand(
    "results/metric/perplexity/salmon/v{version}/{sample}.perplexity.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)


rule metric_tccn_condition_number:
    input:
        gtf="results/download/transcriptome/gencode.v{version}.comprehensive.primary_assembly.annotation.gtf.gz",
    output:
        "results/metric/tccn/v{version}/condition_number.tsv",
    params:
        vendor_dir="workflow/scripts/vendor_tccn",
        scratch_file=lambda wc: f"results/metric/tccn/v{wc.version}/.scratch_gene.txt",
        read_length=235,  # K-value's default
        min_isoforms=2,
        max_isoforms=15,  # 2^n - 1 combinations per gene
        max_genes=1000,
        per_gene_timeout_sec=15,
    conda:
        "../envs/standalone/pandas.yaml"
    log:
        "logs/metric/tccn/v{version}.log",
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/compute_tccn_cn.py"


tccn_condition_number = expand(
    "results/metric/tccn/v{version}/condition_number.tsv",
    version=config["gencode"]["versions"],
)
