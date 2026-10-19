rule aggregate_counts_oarfish:
    input:
        quant="results/quantify/run_oarfish/{sample}.v{version}.quant",
        gtf="results/download/transcriptome/gencode.v{version}.comprehensive.primary_assembly.annotation.gtf.gz",
    output:
        transcript="results/aggregate/oarfish/v{version}/{sample}.transcript_counts.tsv",
        gene="results/aggregate/oarfish/v{version}/{sample}.gene_counts.tsv",
    params:
        cell_line=lambda wc: wc.sample,
        transcript_id_col_ix=1,
        count_id_col_ix=3,
        length_col_ix=2,
    conda:
        "../envs/standalone/pandas.yaml"
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/aggregate_counts.py"


rule aggregate_all_counts_oarfish:
    input:
        expand(
            "results/aggregate/oarfish/v{version}/{sample}.gene_counts.tsv",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),


rule aggregate_counts_salmon:
    input:
        quant="results/quantify/run_salmon_illumina/v{version}/{sample}/quant.sf",
        gtf="results/download/transcriptome/gencode.v{version}.comprehensive.primary_assembly.annotation.gtf.gz",
    output:
        transcript="results/aggregate/salmon/v{version}/{sample}.transcript_counts.tsv",
        gene="results/aggregate/salmon/v{version}/{sample}.gene_counts.tsv",
    params:
        cell_line=lambda wc: wc.sample,
        transcript_id_col_ix=1,
        count_id_col_ix=5,
        length_col_ix=2,
    conda:
        "../envs/standalone/pandas.yaml"
    wildcard_constraints:
        version=r"\d+",
    script:
        "../scripts/aggregate_counts.py"


rule aggregate_all_counts_salmon:
    input:
        expand(
            "results/aggregate/salmon/v{version}/{sample}.gene_counts.tsv",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),
