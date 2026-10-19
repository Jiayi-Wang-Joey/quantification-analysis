import gzip
import re
from typing import Dict, List, Optional, Tuple

import pysam
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

MAX_READS_PER_TRACK = 8
MAX_READS_SCAN = 3000
EXON_DISPLAY_TOTAL = 800.0
INTRON_DISPLAY_WIDTH = 25.0
MIN_EXON_DISPLAY = 4.0


def parse_gtf_exons(gtf_path: str, gene_id_prefix: str):
    exons_by_tx: Dict[str, List[Tuple[int, int]]] = {}
    tx_names: Dict[str, str] = {}
    chrom = None
    strand = None
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "exon":
                continue
            attr = fields[8]
            gene_match = re.search(r'gene_id "([^"]+)"', attr)
            if not gene_match or not gene_match.group(1).startswith(gene_id_prefix):
                continue
            tx_match = re.search(r'transcript_id "([^"]+)"', attr)
            name_match = re.search(r'transcript_name "([^"]+)"', attr)
            transcript_id = tx_match.group(1)
            start, end = int(fields[3]), int(fields[4])
            chrom = fields[0]
            strand = fields[6]
            exons_by_tx.setdefault(transcript_id, []).append((start, end))
            if name_match:
                tx_names[transcript_id] = name_match.group(1)
    for tx in exons_by_tx:
        exons_by_tx[tx] = sorted(exons_by_tx[tx])
    return chrom, strand, exons_by_tx, tx_names


def introns_from_exons(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]


def merge_exon_union(exons_by_tx: Dict[str, List[Tuple[int, int]]]) -> List[Tuple[int, int]]:
    all_exons = sorted(e for exons in exons_by_tx.values() for e in exons)
    merged = []
    for start, end in all_exons:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def build_display_transform(exon_union: List[Tuple[int, int]]):
    total_exon_len = sum(e - s + 1 for s, e in exon_union)
    exon_scale = EXON_DISPLAY_TOTAL / max(1, total_exon_len)
    segments = []
    d_cursor = 0.0
    prev_end = None
    for start, end in exon_union:
        if prev_end is not None and start > prev_end + 1:
            segments.append((prev_end + 1, start - 1, d_cursor, d_cursor + INTRON_DISPLAY_WIDTH, "intron"))
            d_cursor += INTRON_DISPLAY_WIDTH
        exon_len = end - start + 1
        disp_len = max(MIN_EXON_DISPLAY, exon_len * exon_scale)
        segments.append((start, end, d_cursor, d_cursor + disp_len, "exon"))
        d_cursor += disp_len
        prev_end = end

    def g_to_d(pos: int) -> float:
        for g_start, g_end, d_start, d_end, _kind in segments:
            if g_start <= pos <= g_end:
                frac = (pos - g_start) / max(1, (g_end - g_start))
                return d_start + frac * (d_end - d_start)
        if pos < segments[0][0]:
            return segments[0][2]
        return segments[-1][3]

    return g_to_d, segments, d_cursor


def cigar_blocks(read) -> List[Tuple[int, int]]:
    blocks = []
    cur = read.reference_start + 1
    block_start = cur
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            cur += length
        elif op == 3:
            blocks.append((block_start, cur - 1))
            cur += length
            block_start = cur
        elif op == 2:
            cur += length
    blocks.append((block_start, cur - 1))
    return blocks


def classify_short_read(blocks: List[Tuple[int, int]], exons_by_tx) -> set:
    read_introns = [(blocks[i][1] + 1, blocks[i + 1][0] - 1) for i in range(len(blocks) - 1)]
    compatible = set()
    for tx, exons in exons_by_tx.items():
        tx_introns = set(introns_from_exons(exons))
        if not all(ri in tx_introns for ri in read_introns):
            continue
        if all(any(es <= bs and be <= ee for es, ee in exons) for bs, be in blocks):
            compatible.add(tx)
    return compatible


def transcript_interval_to_genomic_blocks(exons, strand, t_start, t_end):
    ordered = sorted(exons) if strand == "+" else sorted(exons, reverse=True)
    blocks = []
    cum = 0
    for e_start, e_end in ordered:
        elen = e_end - e_start + 1
        e_t_start, e_t_end = cum + 1, cum + elen
        ov_start, ov_end = max(t_start, e_t_start), min(t_end, e_t_end)
        if ov_start <= ov_end:
            off_start, off_end = ov_start - cum - 1, ov_end - cum - 1
            if strand == "+":
                g1, g2 = e_start + off_start, e_start + off_end
            else:
                g1, g2 = e_end - off_end, e_end - off_start
            blocks.append((min(g1, g2), max(g1, g2)))
        cum += elen
    return sorted(blocks)


def get_short_reads(bam_path, chrom, region_start, region_end, exons_by_tx):
    reads = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    scanned = 0
    for read in bam.fetch(chrom, region_start - 1, region_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        scanned += 1
        if scanned > MAX_READS_SCAN:
            break
        blocks = cigar_blocks(read)
        compatible = classify_short_read(blocks, exons_by_tx)
        if not compatible:
            continue
        reads.append({"blocks": blocks, "ambiguous": len(compatible) > 1})
    bam.close()
    return reads


def get_long_reads(bam_path, transcript_ids, exons_by_tx, strand):
    bam = pysam.AlignmentFile(bam_path, "rb")
    read_to_tx: Dict[str, set] = {}
    read_primary: Dict[str, Tuple[str, int, int]] = {}
    scanned = 0
    for tx in transcript_ids:
        if tx not in bam.references:
            continue
        for aln in bam.fetch(tx):
            if aln.is_unmapped:
                continue
            scanned += 1
            if scanned > MAX_READS_SCAN:
                break
            read_to_tx.setdefault(aln.query_name, set()).add(tx)
            if not aln.is_secondary and not aln.is_supplementary:
                read_primary[aln.query_name] = (tx, aln.reference_start + 1, aln.reference_end)
    bam.close()
    reads = []
    for qname, tx_set in read_to_tx.items():
        if qname not in read_primary:
            continue
        tx, t_start, t_end = read_primary[qname]
        exons = exons_by_tx.get(tx)
        if not exons:
            continue
        blocks = transcript_interval_to_genomic_blocks(exons, strand, t_start, t_end)
        if not blocks:
            continue
        reads.append({"blocks": blocks, "ambiguous": len(tx_set) > 1})
    return reads


def sample_reads(reads, n):
    ambiguous = [r for r in reads if r["ambiguous"]]
    unambiguous = [r for r in reads if not r["ambiguous"]]
    half = n // 2
    chosen = unambiguous[:half] + ambiguous[: n - half]
    if len(chosen) < n:
        pool = unambiguous[half:] + ambiguous[n - half :]
        chosen += pool[: n - len(chosen)]
    return chosen


def draw_blocks(ax, blocks, g_to_d, y, height, color):
    d_blocks = [(g_to_d(s), g_to_d(e)) for s, e in blocks]
    for i in range(len(d_blocks) - 1):
        x1 = d_blocks[i][1]
        x2 = d_blocks[i + 1][0]
        ax.plot([x1, x2], [y + height / 2, y + height / 2], color=color, linewidth=1)
    for d_start, d_end in d_blocks:
        ax.add_patch(
            Rectangle((d_start, y), max(d_end - d_start, 0.5), height, facecolor=color, edgecolor="none")
        )


def plot_gene_page(pdf, gene_id, gene_name, kvalue, chrom, strand, exons_by_tx, tx_names, short_reads, long_reads):
    exon_union = merge_exon_union(exons_by_tx)
    g_to_d, _segments, total_width = build_display_transform(exon_union)

    tx_order = sorted(exons_by_tx.keys(), key=lambda t: exons_by_tx[t][0][0])
    n_isoforms = len(tx_order)
    n_short = len(short_reads)
    n_long = len(long_reads)

    row_h = 1.0
    fig_h = 1.5 + n_isoforms * row_h + (n_short + n_long) * row_h + 1.5
    fig, ax = plt.subplots(figsize=(11, max(4, fig_h * 0.35)))

    y = 0.0
    for tx in reversed(tx_order):
        exons = exons_by_tx[tx]
        draw_blocks(ax, exons, g_to_d, y, row_h * 0.6, "navy")
        label = tx_names.get(tx, tx)
        ax.text(-10, y + row_h * 0.3, label, ha="right", va="center", fontsize=6)
        y += row_h
    isoform_top = y

    y += 0.5
    ax.axhline(y, color="black", linestyle="--", linewidth=1)
    y += 0.5
    ax.text(-10, y + (n_short * row_h) / 2, "Short\nreads", ha="right", va="center", fontsize=9)
    for r in short_reads:
        color = "grey" if r["ambiguous"] else "#f5c518"
        draw_blocks(ax, r["blocks"], g_to_d, y, row_h * 0.6, color)
        y += row_h
    short_bottom = y

    y += 0.5
    ax.axhline(y, color="black", linestyle="--", linewidth=1)
    y += 0.5
    ax.text(-10, y + (n_long * row_h) / 2, "Long\nreads", ha="right", va="center", fontsize=9)
    for r in long_reads:
        color = "grey" if r["ambiguous"] else "#f5c518"
        draw_blocks(ax, r["blocks"], g_to_d, y, row_h * 0.6, color)
        y += row_h
    long_bottom = y

    ax.set_xlim(-120, total_width + 200)
    ax.set_ylim(-0.5, long_bottom + 0.5)
    ax.set_title(f"{gene_name} (K-value = {kvalue:.2f})  |  {chrom}, {strand} strand", fontsize=12, style="italic")
    ax.axis("off")

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor="grey", label="Ambiguous reads"),
        Rectangle((0, 0), 1, 1, facecolor="#f5c518", label="Unambiguous reads"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=8)

    pdf.savefig(fig)
    plt.close(fig)


def main():
    genes = snakemake.params.genes
    gtf_path = snakemake.input.gtf
    kinnex_bam = snakemake.input.kinnex_bam
    illumina_bam = snakemake.input.illumina_bam
    output_path = snakemake.output[0]

    with PdfPages(output_path) as pdf:
        for gene in genes:
            gene_id, gene_name, kvalue = gene["gene_id"], gene["gene_name"], gene["kvalue"]
            chrom, strand, exons_by_tx, tx_names = parse_gtf_exons(gtf_path, gene_id)
            if chrom is None:
                continue
            exon_union = merge_exon_union(exons_by_tx)
            region_start, region_end = exon_union[0][0], exon_union[-1][1]

            short_reads_all = get_short_reads(illumina_bam, chrom, region_start, region_end, exons_by_tx)
            long_reads_all = get_long_reads(kinnex_bam, list(exons_by_tx.keys()), exons_by_tx, strand)

            short_reads = sample_reads(short_reads_all, MAX_READS_PER_TRACK)
            long_reads = sample_reads(long_reads_all, MAX_READS_PER_TRACK)

            plot_gene_page(
                pdf, gene_id, gene_name, kvalue, chrom, strand, exons_by_tx, tx_names, short_reads, long_reads
            )


main()
