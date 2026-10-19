def gencode_url(wildcards):
    suffix = "" if wildcards.kind == "comprehensive" else ".basic"
    return (
        f"{config['gencode']['base_url']}/release_{wildcards.version}/"
        f"gencode.v{wildcards.version}.primary_assembly{suffix}.annotation.gtf.gz"
    )


rule download_transcriptome:
    output:
        "results/download/transcriptome/gencode.v{version}.{kind}.primary_assembly.annotation.gtf.gz",
    params:
        url=gencode_url,
    log:
        "logs/download/transcriptome/v{version}_{kind}.log",
    wildcard_constraints:
        version=r"\d+",
        kind="comprehensive|basic",
    shell:
        """
        wget -O {output} {params.url} &> {log}
        """


rule download_all_transcriptomes:
    input:
        expand(
            "results/download/transcriptome/gencode.v{version}.{kind}.primary_assembly.annotation.gtf.gz",
            version=config["gencode"]["versions"],
            kind=config["gencode"]["kinds"],
        ),


rule download_genome:
    output:
        "results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    params:
        url=config["genome_url"],
        download_path="results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    log:
        "logs/download/genome/out.log",
    conda:
        "../envs/standalone/curl.yaml"
    shell:
        """
        wget -O {params.download_path} {params.url} &> {log};
        gunzip {params.download_path}
        """


rule download_miniquant:
    output:
        "results/tools/miniquant/miniQuant_linux/miniQuant",
    params:
        url="https://github.com/Augroup/miniQuant/releases/download/v1.4.1/miniQuant_linux_latest.tar.gz",
        download_path="results/tools/miniquant/miniQuant_linux_latest.tar.gz",
        extract_dir="results/tools/miniquant",
    log:
        "logs/download/miniquant/out.log",
    conda:
        "../envs/standalone/curl.yaml"
    shell:
        """
        wget -O {params.download_path} {params.url} &> {log};
        tar -zxvf {params.download_path} -C {params.extract_dir} &>> {log};
        chmod +x {output};
        rm {params.download_path}
        """


rule download_longbench_pacbio:
    output:
        f"results/download/longbench/{config['longbench']['pacbio_file']}",
    params:
        url=f"{config['longbench']['base_url']}/{config['longbench']['pacbio_file']}",
    log:
        "logs/download/longbench/pacbio_out.log",
    shell:
        """
        wget -c -O {output} {params.url} &> {log}
        gzip -t {output} &>> {log}
        """


rule download_longbench_shortread:
    output:
        r1=f"results/download/longbench/{config['longbench']['shortread_r1']}",
        r2=f"results/download/longbench/{config['longbench']['shortread_r2']}",
    params:
        r1_url=f"{config['longbench']['base_url']}/{config['longbench']['shortread_r1']}",
        r2_url=f"{config['longbench']['base_url']}/{config['longbench']['shortread_r2']}",
    log:
        "logs/download/longbench/shortread_out.log",
    shell:
        """
        wget -c -O {output.r1} {params.r1_url} &> {log}
        wget -c -O {output.r2} {params.r2_url} &>> {log}
        gzip -t {output.r1} &>> {log}
        gzip -t {output.r2} &>> {log}
        """
