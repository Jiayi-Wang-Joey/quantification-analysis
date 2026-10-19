_gene_track_sample = config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")


rule gene_track_plot:
    input:
        gtf="results/download/transcriptome/gencode.v50.comprehensive.primary_assembly.annotation.gtf.gz",
        kinnex_bam=f"results/align/run_kinnex/{_gene_track_sample}.v50.kinnex.sorted.bam",
        kinnex_bai=f"results/align/run_kinnex/{_gene_track_sample}.v50.kinnex.sorted.bam.bai",
        illumina_bam=f"results/align/run_illumina_star/v50/{_gene_track_sample}/{_gene_track_sample}_Aligned.sortedByCoord.out.bam",
        illumina_bai=f"results/align/run_illumina_star/v50/{_gene_track_sample}/{_gene_track_sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "plts/gene_track-outliers.pdf",
    params:
        genes=[
            {"gene_id": "ENSG00000198840", "gene_name": "MT-ND3 (low K-value, high MARD)", "kvalue": 1.0},
            {"gene_id": "ENSG00000183049", "gene_name": "CAMK1D (high K-value, low MARD)", "kvalue": 1049.64},
        ],
    conda:
        "../envs/standalone/pysam.yaml"
    log:
        "logs/gene_track/outliers.log",
    script:
        "../scripts/gene_track_plot.py"
