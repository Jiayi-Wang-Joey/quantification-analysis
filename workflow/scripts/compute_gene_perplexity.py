import numpy as np
import pandas as pd
from typeguard import typechecked

# perplexity = 2^H over isoform shares of the gene's TPM (Schertzer et al. 2025)


@typechecked
def compute_gene_perplexity(
    transcript_counts_path: str,
    output_path: str,
    sample: str,
    platform: str,
    gencode_version: str,
) -> int:
    df = pd.read_csv(transcript_counts_path, sep="\t")
    tpm_col = f"{sample}_tpm"
    df = df[["gene_id", "transcript_id", tpm_col]].rename(columns={tpm_col: "tpm"})
    df = df[df["tpm"] > 0].copy()

    gene_total = df.groupby("gene_id")["tpm"].transform("sum")
    p = df["tpm"] / gene_total
    df["plogp"] = p * np.log2(p)

    entropy = df.groupby("gene_id")["plogp"].sum() * -1
    gene_potential = df.groupby("gene_id").size()

    out = pd.DataFrame({
        "gene_id": entropy.index,
        "gene_potential": gene_potential.reindex(entropy.index).values,
        "gene_perplexity": (2.0 ** entropy).values,
    })
    out["platform"] = platform
    out["sample"] = sample
    out["gencode_version"] = gencode_version
    out.to_csv(path_or_buf=output_path, sep="\t", header=True, index=False)
    return 0


compute_gene_perplexity(
    transcript_counts_path=snakemake.input.transcript_counts,
    output_path=snakemake.output[0],
    sample=snakemake.params.sample,
    platform=snakemake.params.platform,
    gencode_version=snakemake.wildcards.version,
)
