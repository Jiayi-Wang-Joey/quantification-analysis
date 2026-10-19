# Vendored from TCCN

`Process_gene_data.py` and `ComputeC.py` are copied unmodified from
https://github.com/HUSTzoulab/TCCN (no LICENSE file present in that repo;
vendored here for internal pipeline use with attribution), accompanying:

Zhang, Guo, Liu, Zou. "Quantifying transcript complexity via the condition
number of gene-specific random matrix." Briefings in Bioinformatics 27(2),
2026. https://academic.oup.com/bib/article/27/2/bbag134/8550916

Condition number recipe (from that repo's `Example.ipynb`), used as-is in
`../compute_tccn_cn.py`:

```python
C = computeC(arrays, read_length, gene_dict)
eigenvalues, _ = np.linalg.eig(C.T @ C)
CN = max(eigenvalues) / min(eigenvalues)
```

Note this is `(sigma_max/sigma_min)^2` (eigenvalues of `C^T C`), not the
direct singular-value ratio K-value uses -- take `sqrt(CN)` before comparing
to K-value on the same scale.

`computeC`'s combination matrix has `2^n - 1` rows for a gene with `n`
isoforms (every non-empty subset of isoforms, via
`generate_combinations_matrix`) -- infeasible much past ~20-25 isoforms.
`compute_tccn_cn.py` caps `n_isoforms` accordingly.
