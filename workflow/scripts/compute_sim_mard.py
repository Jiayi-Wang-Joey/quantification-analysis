# MARD of one simulated sample against its true TPM (miniQuant): TPM < min_tpm set to 0,
# ARD = |est - true| / (est + true); *_eval columns: transcripts >= min_length, TPM renormalised

import numpy as np
import pandas as pd

unver = lambda ids: ids.str.replace(r"\.\d+(?=(_PAR_Y)?$)", "", regex=True)
p = snakemake.params
MIN_TPM = p.min_tpm

ab = pd.read_csv(snakemake.input.abundance, sep="\t", usecols=["target_id", "length", "tpm"])
ab.columns = ["transcript_id", "length", "tpm_true"]

if p.platform == "kinnex":
    est = pd.read_csv(snakemake.input.estimate, sep="\t", usecols=["tname", "len", "num_reads"])
    est.columns = ["transcript_id", "length", "reads"]
    est["tpm_est"] = est["reads"] / est["reads"].sum() * 1e6
else:
    est = pd.read_csv(snakemake.input.estimate, sep="\t", usecols=["Name", "Length", "TPM"])
    est.columns = ["transcript_id", "length", "tpm_est"]

reads = pd.read_csv(snakemake.input.true_counts, sep="\t")[["transcript_id", "true_reads"]]
t2g_truth = pd.read_csv(snakemake.input.tx2gene_truth, sep="\t")
t2g_ref = pd.read_csv(snakemake.input.tx2gene_ref, sep="\t")
for d in (ab, est, reads, t2g_truth, t2g_ref):
    d["key"] = unver(d["transcript_id"])
gene_truth = dict(zip(t2g_truth["key"], unver(t2g_truth["gene_id"])))
gene_ref = dict(zip(t2g_ref["key"], unver(t2g_ref["gene_id"])))
in_ref = set(t2g_ref["key"])

d = ab[["key", "length", "tpm_true"]].merge(est[["key", "tpm_est"]], on="key", how="outer")
d["length"] = d["length"].fillna(d["key"].map(dict(zip(est["key"], est["length"]))))
d[["tpm_true", "tpm_est"]] = d[["tpm_true", "tpm_est"]].fillna(0)
d["true_reads"] = d["key"].map(dict(zip(reads["key"], reads["true_reads"]))).fillna(0)


def score(d, true_col, est_col, suffix):
    """truncate TPM < MIN_TPM, then status-independent ARD"""
    t = d[true_col].where(d[true_col] >= MIN_TPM, 0)
    e = d[est_col].where(d[est_col] >= MIN_TPM, 0)
    d[f"tpm_true{suffix}"], d[f"tpm_est{suffix}"] = t, e
    d[f"ard{suffix}"] = np.where(t + e > 0, (e - t).abs() / (t + e), np.nan)


score(d, "tpm_true", "tpm_est", "")
d = d[(d["tpm_true"] > 0) | (d["tpm_est"] > 0)].copy()
d["status"] = np.where(d["tpm_true"] == 0, "not_expressed",
                       np.where(d["key"].isin(in_ref), "in_reference", "missing_from_reference"))
d["gene"] = d["key"].map(gene_ref).fillna(d["key"].map(gene_truth))
k = pd.read_csv(snakemake.input.kvalue, sep="\t", usecols=["Gene_id", "K-value"])
d["kvalue"] = d["gene"].map(dict(zip(unver(k["Gene_id"]), k["K-value"])))

d["in_eval"] = d["length"] >= p.min_length
e = d["in_eval"]
d["tpm_true_eval_raw"] = np.where(e, d["tpm_true"] / d.loc[e, "tpm_true"].sum() * 1e6, np.nan)
d["tpm_est_eval_raw"] = np.where(e, d["tpm_est"] / d.loc[e, "tpm_est"].sum() * 1e6, np.nan)
score(d, "tpm_true_eval_raw", "tpm_est_eval_raw", "_eval")
d.loc[~e, ["tpm_true_eval", "tpm_est_eval", "ard_eval"]] = np.nan
d = d.drop(columns=["tpm_true_eval_raw", "tpm_est_eval_raw"])

meta = {"sample": p.sample, "platform": p.platform, "depth": p.depth, "truth": f"v{p.truth}", "reference": f"v{p.ref}"}
d.assign(**meta).to_csv(snakemake.output.transcripts, sep="\t", index=False)


def summarise(d, suffix):
    a, t, e = d[f"ard{suffix}"], d[f"tpm_true{suffix}"], d[f"tpm_est{suffix}"]
    expr = t > 0
    g = d[expr].assign(a=a[expr]).groupby("gene")["a"].mean()
    return {
        f"median_gene_mard{suffix}": g.median(),
        f"mard_expressed{suffix}": a[expr].mean(),
        f"mard_expressed_in_reference{suffix}": a[expr & (d["status"] == "in_reference")].mean(),
        f"mard_tpm_ge_1{suffix}": a[t >= 1].mean(),
        f"pct_estimate_not_expressed{suffix}": 100 * e[d["status"] == "not_expressed"].sum() / e.sum(),
        f"pct_truth_missing_from_reference{suffix}": 100 * t[d["status"] == "missing_from_reference"].sum() / t.sum(),
    }


summary = {**meta, "min_tpm": MIN_TPM, "min_length": p.min_length, "n_transcripts": len(d),
           **summarise(d, ""), **summarise(d[d["in_eval"]], "_eval")}
pd.DataFrame([summary]).to_csv(snakemake.output.summary, sep="\t", index=False)
