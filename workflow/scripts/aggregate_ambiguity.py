from typing import List

import pandas as pd
from typeguard import typechecked


@typechecked
def aggregate_ambiguity(
    quant_path: str,
    ambig_info_path: str,
    transcript_id_col_ix: int,
    unique_col_ix: int,
    ambig_col_ix: int,
    gene_map_path: str,
    output_path_transcript: str,
    output_path_gene: str,
    cell_line: str,
    gencode_version: str,
) -> int:
    quant_df = pd.read_csv(quant_path, sep="\t")
    ambig_df = pd.read_csv(ambig_info_path, sep="\t")
    assert len(quant_df) == len(ambig_df), "quant and ambig_info row counts don't match"

    tx_df = pd.DataFrame(
        {
            "transcript_id": quant_df.iloc[:, (transcript_id_col_ix - 1)].values,
            "unique_reads": ambig_df.iloc[:, (unique_col_ix - 1)].astype(float).values,
            "ambig_reads": ambig_df.iloc[:, (ambig_col_ix - 1)].astype(float).values,
        }
    )

    gene_map = pd.read_csv(gene_map_path, sep="\t")[["gene_id", "gene_name", "transcript_id"]].drop_duplicates()
    tx_df = tx_df.merge(gene_map, how="inner", on="transcript_id")
    tx_df["cell_line"] = cell_line
    tx_df["gencode_version"] = gencode_version
    tx_df["ambiguous_fraction"] = tx_df["ambig_reads"] / (tx_df["ambig_reads"] + tx_df["unique_reads"])
    tx_df = tx_df[
        [
            "gene_id", "gene_name", "transcript_id", "cell_line", "gencode_version",
            "unique_reads", "ambig_reads", "ambiguous_fraction",
        ]
    ]
    tx_df.to_csv(path_or_buf=output_path_transcript, sep="\t", header=True, index=False)

    gene_df = tx_df.groupby(["gene_id", "gene_name", "cell_line", "gencode_version"], as_index=False)[
        ["unique_reads", "ambig_reads"]
    ].sum()
    gene_df["ambiguous_fraction"] = gene_df["ambig_reads"] / (gene_df["ambig_reads"] + gene_df["unique_reads"])
    gene_df.to_csv(path_or_buf=output_path_gene, sep="\t", header=True, index=False)

    return 0


aggregate_ambiguity(
    quant_path=snakemake.input.quant,
    ambig_info_path=snakemake.input.ambig_info,
    transcript_id_col_ix=snakemake.params.transcript_id_col_ix,
    unique_col_ix=snakemake.params.unique_col_ix,
    ambig_col_ix=snakemake.params.ambig_col_ix,
    gene_map_path=snakemake.input.gene_map,
    output_path_transcript=snakemake.output.transcript,
    output_path_gene=snakemake.output.gene,
    cell_line=snakemake.params.cell_line,
    gencode_version=snakemake.wildcards.version,
)
