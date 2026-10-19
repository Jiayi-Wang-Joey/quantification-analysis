rule rename_genome_contigs:
    input:
        "results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.gencode_names.fa",
    log:
        "logs/download/genome/rename_contigs.log",
    shell:
        """
        awk '
        /^>/ {{
            line = $0
            split(line, toks, " ")
            first = substr(toks[1], 2)
            rl = ""
            if (match(line, /rl:[A-Za-z]+/)) {{
                rl = substr(line, RSTART+3, RLENGTH-3)
            }}
            ac = ""
            if (match(line, /AC:[^ ]+/)) {{
                ac = substr(line, RSTART+3, RLENGTH-3)
            }}
            if (rl == "Chromosome" || rl == "Mitochondrion") {{
                print ">" first
            }} else {{
                print ">" ac
            }}
            next
        }}
        {{ print }}
        ' {input} > {output} 2> {log}
        """


rule extract_transcriptome:
    input:
        transcriptome="results/download/transcriptome/gencode.v{version}.comprehensive.primary_assembly.annotation.gtf.gz",
        genome="results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.gencode_names.fa",
    output:
        "results/transcriptome/gencode.v{version}.transcriptome.fa",
    log:
        "logs/transcriptome/v{version}.log",
    conda:
        "../envs/align.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        zcat {input.transcriptome} | gffread - -w {output} -g {input.genome} &> {log}
        """


rule align_kinnex:
    input:
        transcriptome="results/transcriptome/gencode.v{version}.transcriptome.fa",
        reads=lambda wc: pacbio_reads(wc),
    output:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.bam",
    log:
        "logs/align/{sample}.v{version}.log",
    conda:
        "../envs/align.yaml"
    threads: config["align"]["map_threads"]
    resources:
        mem_mb=config["align"]["mem_mb"],
    params:
        n_secondary_alignments=config["align"]["n_secondary_alignments"],
        align_sort_bam_threads=config["align"]["sort_threads"],
        memory=config["align"]["sort_memory_gb"],
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        minimap2 --eqx -N {params.n_secondary_alignments} -ax map-hifi \
            -t {threads} {input.transcriptome} \
            {input.reads} 2>> {log} | \
            samtools sort -n -@ {params.align_sort_bam_threads} \
            -m{params.memory}g -o {output} - &>> {log}
        """


rule align_kinnex_sort_coord:
    input:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.bam",
    output:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.sorted.bam",
    log:
        "logs/align/{sample}.v{version}.sort_coord.log",
    conda:
        "../envs/align.yaml"
    params:
        threads=config["align"]["sort_threads"],
        memory=config["align"]["sort_memory_gb"],
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        samtools sort -@ {params.threads} -m{params.memory}g -o {output} {input} &> {log}
        """


rule align_kinnex_index:
    input:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.sorted.bam",
    output:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.sorted.bam.bai",
    log:
        "logs/align/{sample}.v{version}.kinnex_index.log",
    conda:
        "../envs/standalone/samtools.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        samtools --version > {log};
        samtools index {input} >> {log}
        """


rule align_all_kinnex:
    input:
        expand(
            "results/align/run_kinnex/{sample}.v{version}.kinnex.sorted.bam.bai",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),


rule star_genome_index:
    input:
        genome="results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.gencode_names.fa",
        gtf="results/download/transcriptome/gencode.v{version}.comprehensive.primary_assembly.annotation.gtf.gz",
    output:
        directory("results/star_index/v{version}"),
    params:
        plain_gtf="results/star_index/v{version}.annotation.gtf",
        read_length=config["illumina_read_length"],
    log:
        "logs/star_index/v{version}.log",
    conda:
        "../envs/standalone/star.yaml"
    threads: config["align_map_bam_threads"]
    resources:
        mem_mb=config["align"]["star_index_mem_mb"],
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        STAR --version > {log};
        zcat {input.gtf} > {params.plain_gtf}
        mkdir -p {output}
        STAR --runMode genomeGenerate --runThreadN {threads} \
            --genomeDir {output} --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {params.plain_gtf} \
            --sjdbOverhang {params.read_length} &>> {log}
        rm {params.plain_gtf}
        """


rule align_run_illumina_star:
    input:
        genome="results/star_index/v{version}",
        reads_first=f"results/download/longbench/{config['longbench']['shortread_r1']}",
        reads_second=f"results/download/longbench/{config['longbench']['shortread_r2']}",
    output:
        "results/align/run_illumina_star/v{version}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    params:
        outfile_prefix="results/align/run_illumina_star/v{version}/{sample}/{sample}_",
    threads: config["align_map_bam_threads"]
    resources:
        mem_mb=config["align"]["star_index_mem_mb"],
    log:
        "logs/align/run_illumina_star/v{version}_{sample}.out",
    conda:
        "../envs/standalone/star.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        STAR --version > {log};
        STAR --genomeDir {input.genome} --readFilesIn {input.reads_first} {input.reads_second} \
            --runThreadN {threads} --outFileNamePrefix {params.outfile_prefix} \
            --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif &>> {log}
        """


rule align_index_illumina_star:
    input:
        "results/align/run_illumina_star/v{version}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/align/run_illumina_star/v{version}/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/align/index_illumina_star/v{version}_{sample}.out",
    conda:
        "../envs/standalone/samtools.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        samtools --version > {log};
        samtools index {input} >> {log}
        """


rule align_all_illumina_star:
    input:
        expand(
            "results/align/run_illumina_star/v{version}/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),


rule align_convert_illumina_mapped_bams_to_fastq:
    input:
        "results/align/run_illumina_star/v{version}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        singletons="results/align/convert_mapped_bams_to_fastq/v{version}/{sample}/{sample}.singletons.fastq.gz",
        first_reads="results/align/convert_mapped_bams_to_fastq/v{version}/{sample}/{sample}.reads_1.fastq.gz",
        second_reads="results/align/convert_mapped_bams_to_fastq/v{version}/{sample}/{sample}.reads_2.fastq.gz",
    params:
        memory=config["align_sort_bam_memory_gb"],
        compression_level=config["compression_level"],
    threads: config["stall_io_threads"]
    log:
        "logs/align/convert_mapped_bams_to_fastq/v{version}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    wildcard_constraints:
        version=r"\d+",
    shell:
        """
        samtools --version > {log};
        samtools sort -n -m{params.memory}g -@ {threads} {input} 2>> {log} | \
            samtools fastq -@ {threads} -1 {output.first_reads} \
            -2 {output.second_reads} -s {output.singletons} -c {params.compression_level} \
            - &>> {log}
        """


rule align_all_illumina_fastq:
    input:
        expand(
            "results/align/convert_mapped_bams_to_fastq/v{version}/{sample}/{sample}.reads_1.fastq.gz",
            sample=[config["longbench"]["pacbio_file"].removesuffix("_bulk_PB.fastq.gz")],
            version=config["gencode"]["versions"],
        ),
