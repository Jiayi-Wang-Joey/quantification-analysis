from typing import List

import pandas as pd
from typeguard import typechecked


@typechecked
def aggregate_counts(
    sample_paths: List[str],
    sample_names: List[str],
    gtf_annotation_path: str,
    output_path_transcript: str,
    output_path_gene: str,
    transcript_id_col_ix: int,
    count_id_col_ix: int,
    length_col_ix: int,
    cell_line: str,
    gencode_version: str,
) -> int:
    df_dict = {}
    tpm_dict = {}
    for sample_path, sample_name in zip(sample_paths, sample_names):
        df = pd.read_csv(sample_path, sep="\t")
        if sample_name == sample_names[0]:
            df_dict["transcript_id"] = df.iloc[
                :, (transcript_id_col_ix - 1)
            ].values.tolist()
        counts = df.iloc[:, (count_id_col_ix - 1)].astype(float)
        lengths = df.iloc[:, (length_col_ix - 1)].astype(float)
        # TPM from raw length for both platforms
        rpk = counts / (lengths / 1000.0)
        tpm = rpk / rpk.sum() * 1e6
        df_dict[sample_name] = counts.values.tolist()
        tpm_dict[f"{sample_name}_tpm"] = tpm.values.tolist()
    df = pd.DataFrame({**df_dict, **tpm_dict})
    tpm_names = [f"{name}_tpm" for name in sample_names]

    gtf_df = pd.read_csv(gtf_annotation_path, sep="\t", header=None, comment="#")
    gtf_df.columns = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    gtf_df = gtf_df[gtf_df["feature"] == "transcript"]
    gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )
    gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
    gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')
    gtf_df["gene_type"] = gtf_df["attribute"].str.extract(r'gene_type "([^"]+)"')
    gtf_df["transcript_type"] = gtf_df["attribute"].str.extract(
        r'transcript_type "([^"]+)"'
    )
    gene_transcript_df = gtf_df[
        ["gene_id", "gene_name", "gene_type", "transcript_id", "transcript_type"]
    ].copy(deep=True)
    del gtf_df

    df = df.merge(
        gene_transcript_df,
        how="inner",
        left_on="transcript_id",
        right_on="transcript_id",
    ).drop_duplicates()

    df["cell_line"] = cell_line
    df["gencode_version"] = gencode_version

    id_cols = [
        "gene_id",
        "gene_name",
        "gene_type",
        "transcript_id",
        "transcript_type",
        "cell_line",
        "gencode_version",
    ]
    df = df[id_cols + sample_names + tpm_names]
    df = df.sort_values(["gene_id", "transcript_id"])
    df.to_csv(path_or_buf=output_path_transcript, sep="\t", header=True, index=False)

    agg_dict = {col: "sum" for col in sample_names + tpm_names}
    for col in ["gene_name", "gene_type", "cell_line", "gencode_version"]:
        agg_dict[col] = "first"
    gene_df = (
        df.drop(columns=["transcript_id", "transcript_type"])
        .groupby("gene_id")
        .agg(agg_dict)
        .reset_index()
    )
    gene_df.to_csv(path_or_buf=output_path_gene, sep="\t", header=True, index=False)
    return 0


aggregate_counts(
    sample_paths=[snakemake.input.quant],
    sample_names=[snakemake.params.cell_line],
    gtf_annotation_path=snakemake.input.gtf,
    output_path_transcript=snakemake.output.transcript,
    output_path_gene=snakemake.output.gene,
    transcript_id_col_ix=snakemake.params.transcript_id_col_ix,
    count_id_col_ix=snakemake.params.count_id_col_ix,
    length_col_ix=snakemake.params.length_col_ix,
    cell_line=snakemake.params.cell_line,
    gencode_version=snakemake.wildcards.version,
)
