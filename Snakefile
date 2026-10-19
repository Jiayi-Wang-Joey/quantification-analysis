configfile: "config/config.yaml"
configfile: "config/tksm.yaml"


include: "workflow/rules/download.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/quantify.smk"
include: "workflow/rules/metric.smk"
include: "workflow/rules/aggregate.smk"
include: "workflow/rules/gene_track.smk"
include: "workflow/rules/simulate.smk"

quant = expand(
    "results/quantify/run_oarfish/{sample}.v{version}.quant",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

oarfish_gene_counts = expand(
    "results/aggregate/oarfish/v{version}/{sample}.gene_counts.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

salmon_gene_counts = expand(
    "results/aggregate/salmon/v{version}/{sample}.gene_counts.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

oarfish_transcript_counts = expand(
    "results/aggregate/oarfish/v{version}/{sample}.transcript_counts.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

salmon_transcript_counts = expand(
    "results/aggregate/salmon/v{version}/{sample}.transcript_counts.tsv",
    sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
    version=config["gencode"]["versions"],
)

eva = {
    "kvalue": kvalue + kvalue_basic,
    "quant": quant,
    "concordance": oarfish_gene_counts + salmon_gene_counts + kvalue,
    "concordance_transcript": oarfish_transcript_counts + salmon_transcript_counts + kvalue + quant,
    "ard_distribution": oarfish_gene_counts + salmon_gene_counts + oarfish_transcript_counts + salmon_transcript_counts + kvalue,
    "data_quality": data_quality,
    "ambiguity": oarfish_gene_ambiguity + salmon_gene_ambiguity + oarfish_gene_counts + salmon_gene_counts + oarfish_transcript_counts + salmon_transcript_counts + kvalue,
    "infrv": oarfish_infrv + salmon_infrv + oarfish_transcript_counts + salmon_transcript_counts + kvalue,
    "perplexity": oarfish_perplexity + salmon_perplexity + oarfish_gene_counts + salmon_gene_counts + kvalue,
    "tccn": tccn_condition_number + kvalue,
    "simulation": sim_dte_plot_files,
    "simreads": sim_reads_plot_files,
    "simmard": sim_mard_plot_files,
}
VAL = eva.keys()


plt = []
for val in VAL:
    x = glob_wildcards("workflow/plots/plt-"+val+"_{x}.R").x
    plt += expand("plts/{val}-{plt}.pdf", val=val, plt=x)

rule all: 
    input: 
        [x for x in eva.values()],
        plt


# figures: plts/{group}-{name}.pdf from workflow/plots/plt-{group}_{name}.R
for val in VAL:
    rule:
        priority: 90
        conda:  "workflow/envs/r_env.yaml"
        input:  expand("workflow/plots/plt-{val}_{{plt}}.R", val=val), x=eva[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/{val}-{{plt}}.pdf", val=val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val=val)
        shell: "R CMD BATCH --no-restore --no-save '--args {params} {output[0]}' {input[0]} {log[0]}"
