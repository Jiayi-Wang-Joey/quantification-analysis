# multinomial molecule counts from the sample's TPM, written as TPM so that
# `tksm transcribe --molecule-count depth` reproduces them; barcode NA (empty crashes tksm polyA)

import numpy as np
import pandas as pd

ab = pd.read_csv(snakemake.input.abundance, sep="\t", usecols=["target_id", "tpm"])
depth = int(snakemake.params.depth)
rng = np.random.default_rng(int(snakemake.params.seed))

counts = rng.multinomial(depth, ab["tpm"].to_numpy() / ab["tpm"].sum())
out = ab.loc[counts > 0, ["target_id"]].copy()
out["tpm"] = counts[counts > 0] / depth * 1e6
out["cb"] = "NA"
out.to_csv(snakemake.output[0], sep="\t", header=False, index=False)

with open(snakemake.log[0], "w") as log:
    print(f"{depth} molecules over {len(out)} transcripts", file=log)
