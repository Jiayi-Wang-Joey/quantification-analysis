import re

SIM = config["simulate_dte_parametric"]
SIM_SAMPLES = [f"{c}_rep{i}" for c in ("control", "treatment") for i in range(1, SIM["n_reps"] + 1)]
READS = config["simulate_reads"]
SIM_READS = SIM["outdir"] + "/reads"
TKSM_MODELS = "results/sim_dte/tksm_models"
TKSM = config["tksm"]
TKSM_MODULES = {"Tsb": "transcribe", "Trc": "truncate", "plA": "polyA", "Shf": "shuffle",
                "Flp": "flip", "Tag": "tag", "Seq": "sequence"}


def depth_label(d):
    """5000000 -> '5M'"""
    return f"{d // 1_000_000}M" if d % 1_000_000 == 0 else str(d)


def depth_value(label):
    """'5M' -> 5000000"""
    return int(label[:-1]) * 1_000_000 if label.endswith("M") else int(label)


KINNEX_DEPTHS = [depth_label(d) for d in READS["kinnex_depths"]]
ILLUMINA_DEPTHS = [depth_label(d) for d in READS["illumina_depths"]]


def sim_sample(wc):
    """'sim50_kinnex5M_control_rep1' -> ('50', 'kinnex', '5M', 'control_rep1'); None for real samples"""
    m = re.fullmatch(r"sim(\d+)_(kinnex|illumina)(\d+M?)_(.+)", wc.sample)
    return m.groups() if m else None


def pacbio_reads(wc):
    if sim_sample(wc):
        tv, _, d, s = sim_sample(wc)
        return f"{SIM_READS}/v{tv}/kinnex_{d}/{s}.kinnex.fastq.gz"
    return f"results/download/longbench/{wc.sample}_bulk_PB.fastq.gz"


def illumina_reads(wc, mate):
    if sim_sample(wc):
        tv, _, d, s = sim_sample(wc)
        return f"{SIM_READS}/v{tv}/illumina_{d}/{s}.illumina_R{mate}.fastq.gz"
    return (
        f"results/align/convert_mapped_bams_to_fastq/v{wc.version}/{wc.sample}/"
        f"{wc.sample}.reads_{mate}.fastq.gz"
    )


def read_seed(wc):
    return SIM["seed"] * 100000 + int(wc.tv) * 100 + SIM_SAMPLES.index(wc.sample) + 1


rule simulate_dte_parametric:
    input:
        gtf=expand(
            "results/download/transcriptome/gencode.v{version}.comprehensive.primary_assembly.annotation.gtf.gz",
            version=SIM["versions"],
        ),
    output:
        design=f"{SIM['outdir']}/de_design.tsv",
        control=expand(SIM["outdir"] + "/v{version}/control_abundance_v{version}.tsv", version=SIM["versions"]),
        treatment=expand(SIM["outdir"] + "/v{version}/treatment_abundance_v{version}.tsv", version=SIM["versions"]),
        replicates=expand(
            SIM["outdir"] + "/v{version}/{sample}_v{version}.tsv", version=SIM["versions"], sample=SIM_SAMPLES
        ),
        labels=expand(SIM["outdir"] + "/v{version}/de_labels_v{version}.tsv", version=SIM["versions"]),
    params:
        p=SIM,
    log:
        "logs/simulate/simulate_dte_parametric.log",
    resources:
        mem_mb=24000,
    conda:
        "../envs/standalone/r_simulate.yaml"
    script:
        "../scripts/simulate_dte_parametric.R"


sim_dte_plot_files = [
    f"{SIM['outdir']}/v{max(SIM['versions'])}/{name}_v{max(SIM['versions'])}.tsv"
    for name in ["control_abundance", "treatment_abundance", "de_labels"] + SIM_SAMPLES
]


rule simulate_tksm_models:
    input:
        reads=f"results/download/longbench/{config['longbench']['pacbio_file']}",
        transcriptome="results/transcriptome/gencode.v50.transcriptome.fa",
    output:
        kde=f"{TKSM_MODELS}/kde_truncation.json",
        error_model=f"{TKSM_MODELS}/error_model",
        qscore_model=f"{TKSM_MODELS}/qscore_model",
        identity=f"{TKSM_MODELS}/identity.txt",
    params:
        n_lines=TKSM["models"]["training_reads"] * 4,
        outdir=TKSM_MODELS,
        grid_end=TKSM["models"]["truncation"]["grid_end"],
        bandwidth=TKSM["models"]["truncation"]["bandwidth"],
        end_ratio=TKSM["models"]["truncation"]["end_ratio"],
        prepare=(
            "awk -F'\\t' -v OFS='\\t' '$5 == \"+\" {$9 = $7; print}'"
            if TKSM["models"]["truncation"]["five_prime_only"] else "cat"
        ),
    threads: TKSM["threads"]
    resources:
        mem_mb=32000,
    log:
        "logs/simulate/tksm_models.log",
    conda:
        "../envs/standalone/tksm.yaml"
    shell:
        """
        (set +o pipefail; gzip -dc {input.reads} | head -n {params.n_lines}) > {params.outdir}/H69_sub.fastq
        minimap2 -c -x map-hifi --secondary=no -t {threads} {input.transcriptome} {params.outdir}/H69_sub.fastq \\
            2> {log} | awk -F'\\t' '/tp:A:P/' > {params.outdir}/H69_sub.paf
        {params.prepare} {params.outdir}/H69_sub.paf > {params.outdir}/H69_sub.truncation_training.paf
        tksm model-truncation -i {params.outdir}/H69_sub.truncation_training.paf -o {output.kde} -t {threads} \\
            --grid-end {params.grid_end} -b {params.bandwidth} --end-ratio {params.end_ratio} &>> {log}
        badread error_model --reference {input.transcriptome} --reads {params.outdir}/H69_sub.fastq \\
            --alignment {params.outdir}/H69_sub.paf > {output.error_model} 2>> {log}
        badread qscore_model --reference {input.transcriptome} --reads {params.outdir}/H69_sub.fastq \\
            --alignment {params.outdir}/H69_sub.paf > {output.qscore_model} 2>> {log}
        python workflow/scripts/vendor_tksm/badread_identity.py {params.outdir}/H69_sub.paf 2>> {log} | tail -1 | \\
            sed -E 's/mean: ([0-9.]+), max: ([0-9.]+), stddev: ([0-9.]+)/\\1,\\2,\\3/' > {output.identity}
        grep -Eq '^[0-9.]+,[0-9.]+,[0-9.]+$' {output.identity}
        """


rule simulate_plain_gtf:
    input:
        "results/download/transcriptome/gencode.v{tv}.comprehensive.primary_assembly.annotation.gtf.gz",
    output:
        temp("results/sim_dte/gtf/gencode.v{tv}.gtf"),
    shell:
        "zcat {input} > {output}"


rule simulate_kinnex_molecules:
    input:
        abundance=SIM["outdir"] + "/v{tv}/{sample}_v{tv}.tsv",
    output:
        SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.kinnex_molecules.tsv",
    params:
        depth=lambda wc: depth_value(wc.depth),
        seed=read_seed,
    log:
        "logs/simulate/kinnex_molecules/v{tv}_{depth}_{sample}.log",
    wildcard_constraints:
        tv=r"\d+",
        depth=r"\d+M?",
        sample="|".join(SIM_SAMPLES),
    conda:
        "../envs/standalone/pandas.yaml"
    script:
        "../scripts/draw_kinnex_molecules.py"


def tksm_commands(wc, input, output, threads):
    steps = [next(iter(m.items())) for m in TKSM["pipeline"]]
    if steps[0][0] != "Tsb" or steps[-1][0] != "Seq":
        raise ValueError("config/tksm.yaml pipeline must start with Tsb and end with Seq")
    tmp = f"{SIM_READS}/v{wc.tv}/kinnex_{wc.depth}/tmp_{wc.sample}"
    log = f"logs/simulate/kinnex_reads/v{wc.tv}_{wc.depth}_{wc.sample}.log"
    seed = read_seed(wc)
    cmds, prev = [], None
    for i, (name, cfg) in enumerate(steps):
        cfg = cfg or {}
        extra = cfg.get("params", "")
        out = f"{tmp}/{i}_{TKSM_MODULES[name]}.mdf"
        if name == "Tsb":
            cmds.append(f"tksm transcribe -g {input.gtf} -a {input.molecules} "
                        f"--molecule-count {depth_value(wc.depth)} {extra} -s {seed} -o {out}")
        elif name == "Seq":
            identity = cfg.get("identity") or f"$(cat {input.identity})"
            cmds.append(f"tksm sequence -i {prev} -r {input.genome} -o {output.reads} "
                        f"--badread-identity {identity} --badread-error-model {input.error_model} "
                        f"--badread-qscore-model {input.qscore_model} --threads {threads} {extra}")
            out = prev
        else:
            model = f"--kde-model {input.kde}" if name == "Trc" else ""
            cmds.append(f"tksm {TKSM_MODULES[name]} -i {prev} -o {out} {model} {extra} -s {seed}")
        prev = out
    cmds = [f"{c} &>> {log}" for c in cmds]
    cmds.append(f"awk -F'\\t' '/^\\+/ {{ id = substr($1, 2); match($3, /tid=[^;]+/); "
                f"print id \"\\t\" substr($3, RSTART + 4, RLENGTH - 4) }}' {prev} | gzip > {output.molecules}")
    return "\n".join(cmds)


rule simulate_kinnex_reads:
    input:
        molecules=SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.kinnex_molecules.tsv",
        gtf="results/sim_dte/gtf/gencode.v{tv}.gtf",
        genome="results/download/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.gencode_names.fa",
        kde=f"{TKSM_MODELS}/kde_truncation.json",
        error_model=f"{TKSM_MODELS}/error_model",
        qscore_model=f"{TKSM_MODELS}/qscore_model",
        identity=f"{TKSM_MODELS}/identity.txt",
    output:
        reads=SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.kinnex.fastq.gz",
        molecules=SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.kinnex_molecule_tid.tsv.gz",
    params:
        tmp=SIM_READS + "/v{tv}/kinnex_{depth}/tmp_{sample}",
        commands=tksm_commands,
    threads: TKSM["threads"]
    resources:
        mem_mb=32000,
    log:
        "logs/simulate/kinnex_reads/v{tv}_{depth}_{sample}.log",
    wildcard_constraints:
        tv=r"\d+",
        depth=r"\d+M?",
        sample="|".join(SIM_SAMPLES),
    conda:
        "../envs/standalone/tksm.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        : > {log}
        {params.commands}
        rm -r {params.tmp}
        """


rule simulate_illumina_reads:
    input:
        abundance=SIM["outdir"] + "/v{tv}/{sample}_v{tv}.tsv",
        transcriptome="results/transcriptome/gencode.v{tv}.transcriptome.fa",
    output:
        r1=SIM_READS + "/v{tv}/illumina_{depth}/{sample}.illumina_R1.fastq.gz",
        r2=SIM_READS + "/v{tv}/illumina_{depth}/{sample}.illumina_R2.fastq.gz",
    params:
        depth=lambda wc: depth_value(wc.depth),
        read_length=READS["illumina_read_length"],
        error_rate=READS["illumina_error_rate"],
        chunk_pairs=READS["illumina_chunk_pairs"],
        seed=read_seed,
        tmpdir=SIM_READS + "/v{tv}/illumina_{depth}/tmp_{sample}",
    resources:
        mem_mb=16000,
    log:
        "logs/simulate/illumina_reads/v{tv}_{depth}_{sample}.log",
    wildcard_constraints:
        tv=r"\d+",
        depth=r"\d+M?",
        sample="|".join(SIM_SAMPLES),
    conda:
        "../envs/standalone/polyester.yaml"
    script:
        "../scripts/simulate_illumina_polyester.R"


rule simulate_true_counts_kinnex:
    input:
        reads=SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.kinnex.fastq.gz",
        molecules=SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.kinnex_molecule_tid.tsv.gz",
    output:
        SIM_READS + "/v{tv}/kinnex_{depth}/{sample}.true_counts.tsv",
    wildcard_constraints:
        tv=r"\d+",
        depth=r"\d+M?",
        sample="|".join(SIM_SAMPLES),
    shell:
        """
        # inputs in order: 1 molecule -> transcript table, 2 Kinnex reads
        awk -F'\\t' -v OFS='\\t' '
            FNR == 1 {{ f++ }}
            f == 1 {{ tid[$1] = $2; next }}
            FNR % 4 == 1 {{ match($0, /molecule_id=[^ ]+/); n[tid[substr($0, RSTART + 12, RLENGTH - 12)]]++ }}
            END {{ print "transcript_id", "true_reads"; for (t in n) print t, n[t] }}' \\
            <(zcat {input.molecules}) <(zcat {input.reads}) > {output}
        """


rule simulate_true_counts_illumina:
    input:
        SIM_READS + "/v{tv}/illumina_{depth}/{sample}.illumina_R1.fastq.gz",
    output:
        SIM_READS + "/v{tv}/illumina_{depth}/{sample}.true_counts.tsv",
    wildcard_constraints:
        tv=r"\d+",
        depth=r"\d+M?",
        sample="|".join(SIM_SAMPLES),
    shell:
        """
        zcat {input} | awk -v OFS='\\t' '
            NR % 4 == 1 {{ split($0, a, "/"); split(a[2], b, ";"); n[b[1]]++ }}
            END {{ print "transcript_id", "true_reads"; for (t in n) print t, n[t] }}' > {output}
        """


READ_SAMPLES = READS.get("samples", SIM_SAMPLES)
sim_kinnex_samples = [f"sim{v}_kinnex{d}_{s}" for v in READS["versions"] for d in KINNEX_DEPTHS for s in READ_SAMPLES]
sim_illumina_samples = [f"sim{v}_illumina{d}_{s}" for v in READS["versions"] for d in ILLUMINA_DEPTHS for s in READ_SAMPLES]


def matched(pattern, samples):
    """fill {sample} and {version}: each simulated sample with its own truth version"""
    return [pattern.format(sample=s, version=re.match(r"sim(\d+)_", s).group(1)) for s in samples]
sim_read_samples = sim_kinnex_samples + sim_illumina_samples


def sim_true_counts(sample):
    tv, platform, d, s = re.fullmatch(r"sim(\d+)_(kinnex|illumina)(\d+M?)_(.+)", sample).groups()
    return f"{SIM_READS}/v{tv}/{platform}_{d}/{s}.true_counts.tsv"


rule simulate_fit_parameters:
    # rough fit of the truth's distributions to real H69
    input:
        gtf="results/download/transcriptome/gencode.v50.comprehensive.primary_assembly.annotation.gtf.gz",
        quant=["results/quantify/run_oarfish/H69.v50.quant", "results/quantify/run_salmon_illumina/v50/H69/quant.sf"],
    output:
        table="results/sim_dte/calibration/fit_parameters_v50.tsv",
        pdf="plts/simulation-calibration.pdf",
    params:
        min_gene_reads=10,
        iso_gene_reads=100,
        min_share=SIM["min_share"],
        current={k: SIM[k] for k in ("sdlog", "lambda", "alpha")},
    resources:
        mem_mb=16000,
    log:
        "logs/simulate/fit_parameters.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/fit_sim_parameters.R"


rule simulate_tx2gene:
    input:
        "results/download/transcriptome/gencode.v{v}.comprehensive.primary_assembly.annotation.gtf.gz",
    output:
        "results/sim_dte/tx2gene/v{v}.tsv",
    wildcard_constraints:
        v=r"\d+",
    shell:
        """
        (echo -e "transcript_id\\tgene_id"
         zcat {input} | awk -F'\\t' '$3 == "transcript" {{
            match($9, /transcript_id "[^"]+"/); t = substr($9, RSTART + 15, RLENGTH - 16)
            match($9, /gene_id "[^"]+"/); g = substr($9, RSTART + 9, RLENGTH - 10)
            print t "\\t" g }}') > {output}
        """


def sim_estimate(wc):
    tv, platform, d, s = sim_sample(wc)
    if platform == "kinnex":
        return f"results/quantify/run_oarfish/{wc.sample}.v{wc.ref}.quant"
    return f"results/quantify/run_salmon_illumina/v{wc.ref}/{wc.sample}/quant.sf"


rule metric_sim_mard:
    input:
        true_counts=lambda wc: sim_true_counts(wc.sample),
        abundance=lambda wc: f"{SIM['outdir']}/v{sim_sample(wc)[0]}/{sim_sample(wc)[3]}_v{sim_sample(wc)[0]}.tsv",
        estimate=sim_estimate,
        tx2gene_truth=lambda wc: f"results/sim_dte/tx2gene/v{sim_sample(wc)[0]}.tsv",
        tx2gene_ref="results/sim_dte/tx2gene/v{ref}.tsv",
        kvalue="results/metric/kvalue/comprehensive/v{ref}/kvalues.tsv",
    output:
        transcripts="results/metric/sim_mard/v{ref}/{sample}.transcripts.tsv",
        summary="results/metric/sim_mard/v{ref}/{sample}.summary.tsv",
    params:
        platform=lambda wc: sim_sample(wc)[1],
        depth=lambda wc: sim_sample(wc)[2],
        ref=lambda wc: wc.ref,
        sample=lambda wc: sim_sample(wc)[3],
        truth=lambda wc: sim_sample(wc)[0],
        min_length=READS["mard_min_length"],
        min_tpm=READS["mard_min_tpm"],
    wildcard_constraints:
        ref=r"\d+",
        sample=r"sim\d+_(kinnex|illumina)\d+M?_.+",
    conda:
        "../envs/standalone/pandas.yaml"
    script:
        "../scripts/compute_sim_mard.py"


sim_mard_summaries = matched("results/metric/sim_mard/v{version}/{sample}.summary.tsv", sim_read_samples)

# inputs of plts/simmard-*.pdf
sim_mard_plot_files = sim_mard_summaries + [f.replace(".summary.tsv", ".transcripts.tsv") for f in sim_mard_summaries]


rule metric_kinnex_read_qc:
    # per-read QC of the first n reads (name-sorted BAM)
    input:
        "results/align/run_kinnex/{sample}.v{version}.kinnex.bam",
    output:
        "results/metric/read_qc/kinnex/v{version}/{sample}.tsv",
    params:
        n=200000,
    log:
        "logs/metric/read_qc/kinnex/v{version}_{sample}.log",
    wildcard_constraints:
        version=r"\d+",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (echo -e "read_length\\ttranscript_length\\tref_start\\tref_end\\tstrand\\tpolya_length\\tnm\\taligned_bases\\tinserted\\tdeleted\\tmean_quality\\tn_alignments"
         set +o pipefail
         samtools view -h {input} 2> {log} | awk -F'\\t' -v OFS='\\t' -v n={params.n} '
            BEGIN {{ for (i = 33; i < 127; i++) ord[sprintf("%c", i)] = i }}
            function emit() {{ if (row != "" && k < n) {{ print row, cnt; row = ""; if (++k >= n) exit }} }}
            /^@SQ/ {{ for (i = 2; i <= NF; i++) {{ if ($i ~ /^SN:/) s = substr($i, 4); if ($i ~ /^LN:/) l = substr($i, 4) }}; len[s] = l; next }}
            /^@/ {{ next }}
            $1 != cur {{ emit(); cur = $1; cnt = 0; row = "" }}
            {{ if (!and($2, 4)) cnt++ }}
            !and($2, 2308) {{
               span = 0; al = 0; ins = 0; del = 0; c = $6
               while (match(c, /^[0-9]+[MIDNSHP=X]/)) {{ op = substr(c, RLENGTH, 1); v = substr(c, 1, RLENGTH - 1) + 0
                   if (op ~ /[MDN=X]/) span += v; if (op ~ /[M=X]/) al += v; if (op == "I") ins += v; if (op == "D") del += v
                   c = substr(c, RLENGTH + 1) }}
               seq = $10; p = 0
               if (and($2, 16)) {{ for (i = 1; i <= length(seq) && substr(seq, i, 1) == "T"; i++) p++ }}
               else {{ for (i = length(seq); i > 0 && substr(seq, i, 1) == "A"; i--) p++ }}
               nm = "NA"; for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6)
               q = "NA"; if ($11 != "*") {{ t = 0; for (i = 1; i <= length($11); i++) t += ord[substr($11, i, 1)] - 33; q = t / length($11) }}
               row = length(seq) OFS len[$3] OFS $4 - 1 OFS $4 - 1 + span OFS (and($2, 16) ? "-" : "+") OFS p OFS nm OFS al OFS ins OFS del OFS q
            }}
            END {{ emit() }}') > {output}
        """


def illumina_qc_reads(wc, mate):
    if sim_sample(wc):
        return illumina_reads(wc, mate)
    return f"results/download/longbench/{config['longbench'][f'shortread_r{mate}']}"


rule metric_illumina_read_qc:
    # first n read pairs aligned to the transcriptome: samtools stats and read-1 positions
    input:
        r1=lambda wc: illumina_qc_reads(wc, 1),
        r2=lambda wc: illumina_qc_reads(wc, 2),
        transcriptome="results/transcriptome/gencode.v{version}.transcriptome.fa",
    output:
        stats="results/metric/read_qc/illumina/v{version}/{sample}.stats.txt",
        positions="results/metric/read_qc/illumina/v{version}/{sample}.positions.tsv",
    params:
        n_lines=500000 * 4,
        tmp="results/metric/read_qc/illumina/v{version}/tmp_{sample}",
    threads: 8
    resources:
        mem_mb=16000,
    log:
        "logs/metric/read_qc/illumina/v{version}_{sample}.log",
    wildcard_constraints:
        version=r"\d+",
    conda:
        "../envs/align.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        [ -s {input.transcriptome}.fai ] || samtools faidx {input.transcriptome} 2> {log}
        (set +o pipefail; zcat {input.r1} | head -n {params.n_lines}) > {params.tmp}/R1.fq
        (set +o pipefail; zcat {input.r2} | head -n {params.n_lines}) > {params.tmp}/R2.fq
        minimap2 -ax sr --secondary=no -t {threads} {input.transcriptome} {params.tmp}/R1.fq {params.tmp}/R2.fq 2>> {log} | \\
            samtools sort -@ {threads} -o {params.tmp}/aln.bam - 2>> {log}
        samtools stats -r {input.transcriptome} {params.tmp}/aln.bam > {output.stats} 2>> {log}
        (echo -e "transcript_length\\tref_start\\tref_end\\tinsert_size"
         samtools view -h -F 2308 -f 64 {params.tmp}/aln.bam | awk -F'\\t' -v OFS='\\t' '
            /^@SQ/ {{ for (i = 2; i <= NF; i++) {{ if ($i ~ /^SN:/) s = substr($i, 4); if ($i ~ /^LN:/) l = substr($i, 4) }}; len[s] = l; next }}
            /^@/ {{ next }}
            {{ span = 0; c = $6
               while (match(c, /^[0-9]+[MIDNSHP=X]/)) {{ op = substr(c, RLENGTH, 1); v = substr(c, 1, RLENGTH - 1) + 0
                   if (op ~ /[MDN=X]/) span += v; c = substr(c, RLENGTH + 1) }}
               print len[$3], $4 - 1, $4 - 1 + span, ($9 < 0 ? -$9 : $9) }}') > {output.positions}
        rm -r {params.tmp}
        """


# inputs of plts/simreads-*.pdf
_K, _I = KINNEX_DEPTHS[-1], ILLUMINA_DEPTHS[-1]
_KS = [s for s in sim_kinnex_samples if f"_kinnex{_K}_" in s]
_IS = [s for s in sim_illumina_samples if f"_illumina{_I}_" in s]
sim_reads_plot_files = (
    [f"results/metric/read_qc/kinnex/v{v}/H69.tsv" for v in READS["versions"]]
    + [f"results/metric/data_quality/kinnex/v{v}/H69.stats.txt" for v in READS["versions"]]
    + [f"results/metric/read_qc/illumina/v{v}/H69.stats.txt" for v in READS["versions"]]
    + [f"results/quantify/run_salmon_illumina/v{v}/H69/aux_info/meta_info.json" for v in READS["versions"]]
    + matched("results/metric/read_qc/kinnex/v{version}/{sample}.tsv", _KS)
    + matched("results/metric/data_quality/kinnex/v{version}/{sample}.stats.txt", _KS)
    + matched("results/quantify/run_oarfish/{sample}.v{version}.quant", _KS)
    + matched("results/metric/read_qc/illumina/v{version}/{sample}.stats.txt", [s for s in _IS if s.endswith("control_rep1")])
    + matched("results/metric/read_qc/illumina/v{version}/{sample}.positions.tsv", [s for s in _IS if s.endswith("control_rep1")])
    + [f"results/metric/read_qc/illumina/v{v}/H69.positions.tsv" for v in READS["versions"]]
    + matched("results/quantify/run_salmon_illumina/v{version}/{sample}/aux_info/meta_info.json", _IS)
    + [sim_true_counts(s) for s in _KS + _IS if s.endswith("control_rep1")]
    + [f"{SIM['outdir']}/v{v}/control_rep1_v{v}.tsv" for v in READS["versions"]]
)
