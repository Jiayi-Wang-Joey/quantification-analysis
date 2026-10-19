import gzip
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, snakemake.params.vendor_dir)
from Process_gene_data import process_gene_data
from ComputeC import computeC

GTF = snakemake.input.gtf
OUTPUT = snakemake.output[0]
SCRATCH = snakemake.params.scratch_file
READ_LENGTH = snakemake.params.read_length
MIN_ISOFORMS = snakemake.params.min_isoforms
MAX_ISOFORMS = snakemake.params.max_isoforms
GENCODE_VERSION = snakemake.wildcards.version

# GTF lines per gene, in TCCN's input format
gene_lines = defaultdict(list)
n_transcripts = defaultdict(int)

with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        feature = fields[2]
        if feature not in ("gene", "transcript", "exon"):
            continue
        gene_id = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
        start, end = fields[3], fields[4]
        if feature == "gene":
            gene_lines[gene_id].insert(0, f'gene {start} {end} gene_id "{gene_id}";\n')
        elif feature == "transcript":
            tid = re.search(r'transcript_id "([^"]+)"', fields[8]).group(1)
            gene_lines[gene_id].append(f'transcript {start} {end} gene_id "{gene_id}"; transcript_id "{tid}";\n')
            n_transcripts[gene_id] += 1
        elif feature == "exon":
            tid = re.search(r'transcript_id "([^"]+)"', fields[8]).group(1)
            gene_lines[gene_id].append(f'exon {start} {end} gene_id "{gene_id}"; transcript_id "{tid}";\n')

# 2^n - 1 combinations per gene
candidate_genes = [g for g, n in n_transcripts.items() if MIN_ISOFORMS <= n <= MAX_ISOFORMS]
print(f"total genes: {len(gene_lines)}, candidates with {MIN_ISOFORMS}-{MAX_ISOFORMS} isoforms: {len(candidate_genes)}", flush=True)

# random subset of genes: TCCN is slow per gene
MAX_GENES = getattr(snakemake.params, "max_genes", None)
FORCE_INCLUDE = ["ENSG00000168005"]  # SPINDOC
if MAX_GENES is not None and len(candidate_genes) > MAX_GENES:
    import random
    random.seed(0)
    forced = [g for g in candidate_genes if any(g.startswith(p) for p in FORCE_INCLUDE)]
    rest = [g for g in candidate_genes if g not in forced]
    candidate_genes = forced + random.sample(rest, MAX_GENES - len(forced))
    print(f"sampled down to {len(candidate_genes)} genes ({len(forced)} forced-included)", flush=True)

import signal
import time


class GeneTimeout(Exception):
    pass


def _timeout_handler(signum, frame):
    raise GeneTimeout()


PER_GENE_TIMEOUT_SEC = getattr(snakemake.params, "per_gene_timeout_sec", 15)
signal.signal(signal.SIGALRM, _timeout_handler)

results = []
n_timeout = 0
for i, gene_id in enumerate(candidate_genes):
    n_iso = n_transcripts[gene_id]
    t0 = time.time()
    with open(SCRATCH, "w") as out:
        out.writelines(gene_lines[gene_id])
    signal.alarm(PER_GENE_TIMEOUT_SEC)
    try:
        arrays, gene_dict = process_gene_data(SCRATCH)
        C = computeC(arrays, READ_LENGTH, gene_dict)
        if isinstance(C, str):  # tool returns an error string for some inputs
            continue
        eigenvalues, _ = np.linalg.eig(C.T @ C)
        eigenvalues = eigenvalues.real
        eigenvalues = eigenvalues[eigenvalues > 1e-12]
        if len(eigenvalues) == 0:
            continue
        cn = eigenvalues.max() / eigenvalues.min()
        results.append((gene_id, n_iso, cn))
    except GeneTimeout:
        n_timeout += 1
        print(f"  [{i+1}] {gene_id} (n_iso={n_iso}) TIMED OUT after {PER_GENE_TIMEOUT_SEC}s, skipping", flush=True)
        continue
    except Exception as e:
        print(f"  [{i+1}] {gene_id} (n_iso={n_iso}) failed: {e}", flush=True)
        continue
    finally:
        signal.alarm(0)
    elapsed = time.time() - t0
    if elapsed > 1.0:
        print(f"  [{i+1}] {gene_id} (n_iso={n_iso}) took {elapsed:.1f}s", flush=True)
    if (i + 1) % 20 == 0:
        print(f"{i+1}/{len(candidate_genes)} processed, {len(results)} succeeded, {n_timeout} timed out", flush=True)

with open(OUTPUT, "w") as f:
    f.write("gene_id\tn_isoforms\tCN\tread_length\tgencode_version\n")
    for gene_id, n, cn in results:
        f.write(f"{gene_id}\t{n}\t{cn}\t{READ_LENGTH}\t{GENCODE_VERSION}\n")

print(f"Done: {len(results)}/{len(candidate_genes)} genes succeeded. Written to {OUTPUT}", flush=True)
