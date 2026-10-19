import pandas as pd
from typeguard import typechecked

# fishpond::computeInfRV defaults: InfRV = max(var - mu, 0) / (mu + pc) + shift
PC = 5
SHIFT = 0.01


@typechecked
def compute_oarfish_infrv(
    quant_path: str,
    infreps_path: str,
    output_path: str,
    sample: str,
    gencode_version: str,
) -> int:
    quant_df = pd.read_csv(quant_path, sep="\t")
    reps = pd.read_parquet(infreps_path)
    assert len(quant_df) == len(reps), "quant and infreps row counts don't match"

    mu = reps.mean(axis=1)
    inf_var = reps.var(axis=1, ddof=1)
    inf_rv = (inf_var - mu).clip(lower=0) / (mu + PC) + SHIFT

    out = pd.DataFrame(
        {
            "transcript_id": quant_df["tname"].values,
            "mean_count": mu.values,
            "inf_var": inf_var.values,
            "InfRV": inf_rv.values,
            "platform": "oarfish",
            "sample": sample,
            "gencode_version": gencode_version,
        }
    )
    out.to_csv(path_or_buf=output_path, sep="\t", header=True, index=False)
    return 0


compute_oarfish_infrv(
    quant_path=snakemake.input.quant,
    infreps_path=snakemake.input.infreps,
    output_path=snakemake.output[0],
    sample=snakemake.params.sample,
    gencode_version=snakemake.wildcards.version,
)
