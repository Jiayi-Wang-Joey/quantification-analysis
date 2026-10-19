suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


paths <- args[[1]]

kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]
oarfish_paths <- paths[grepl("/aggregate/oarfish/", paths)]
salmon_paths <- paths[grepl("/aggregate/salmon/", paths)]

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, gencode_version := as.integer(gencode_version)]
k <- k[, .(gene_id, gencode_version, `K-value`)]

oarfish <- rbindlist(lapply(oarfish_paths, fread, header = TRUE))
setnames(oarfish, c("H69", "H69_tpm"), c("kinnex", "kinnex_tpm"))
oarfish <- oarfish[, .(gene_id, gencode_version, kinnex, kinnex_tpm)]

salmon <- rbindlist(lapply(salmon_paths, fread, header = TRUE))
setnames(salmon, c("H69", "H69_tpm"), c("illumina", "illumina_tpm"))
salmon <- salmon[, .(gene_id, gencode_version, illumina, illumina_tpm)]

merged <- merge(oarfish, salmon, by = c("gene_id", "gencode_version"))
merged <- merge(merged, k, by = c("gene_id", "gencode_version"))

merged <- merged[`K-value` > 90]
# TPM > 1 in both platforms (miniQuant)
merged <- merged[illumina_tpm > 1 & kinnex_tpm > 1]
merged[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]


merged[, ard := abs(illumina - kinnex) / (illumina + kinnex)]

cor_dt <- merged[, .(
    rho = cor(`K-value`, ard, method = "spearman"),
    x = min(`K-value`),
    y = max(ard)
), by = gencode_version]
cor_dt[, label := sprintf("rho = %.2f", rho)]

p <- ggplot(merged, aes(x = `K-value`, y = ard)) +
    geom_point_rast(alpha = 0.4, color = "steelblue") +
    geom_smooth(method = "loess", color = "black", se = TRUE) +
    geom_text(
        data = cor_dt,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
    ) +
    facet_wrap(~gencode_version, ncol = 4) +
    scale_x_log10() +
    labs(
        x = "K-value (log10 scale)",
        y = "Absolute relative difference: |Illumina - Kinnex| / (Illumina + Kinnex)",
        title = "Illumina vs Kinnex disagreement vs K-value (K-value > 90)"
    ) +
    theme_classic()

ggsave(args[[2]], plot = p, width = 9, height = 6)
