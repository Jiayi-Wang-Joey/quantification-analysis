suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
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
oarfish <- oarfish[, .(transcript_id, gencode_version, kinnex, kinnex_tpm)]

salmon <- rbindlist(lapply(salmon_paths, fread, header = TRUE))
setnames(salmon, c("H69", "H69_tpm"), c("illumina", "illumina_tpm"))
salmon <- salmon[, .(gene_id, transcript_id, gencode_version, illumina, illumina_tpm)]

merged <- merge(oarfish, salmon, by = c("transcript_id", "gencode_version"))
merged <- merge(merged, k, by = c("gene_id", "gencode_version"))

# TPM > 1 in both platforms (miniQuant)
merged <- merged[illumina_tpm > 1 & kinnex_tpm > 1]

# per-transcript ARD between platforms; K-value of the gene
merged[, mard := abs(illumina - kinnex) / ((illumina + kinnex) / 2)]

KVALUE_BREAKS <- c(0, 1, 3, 5, 10, 25, 90, Inf)
KVALUE_LABELS <- c("<=1", "1-3", "3-5", "5-10", "10-25", "25-90", ">90")
merged[, kvalue_bin := cut(
    `K-value`, breaks = KVALUE_BREAKS, labels = KVALUE_LABELS, right = TRUE, include.lowest = TRUE
)]

merged[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]

p <- ggplot(merged, aes(x = kvalue_bin, y = mard, color = kvalue_bin)) +
    geom_violin(fill = "white", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
    scale_color_brewer(palette = "Paired") +
    facet_wrap(~gencode_version, ncol = 4) +
    labs(
        x = "K-value bin (gene-level, broadcast per transcript)",
        y = "Transcript-level MARD: |Illumina - Kinnex| / mean(Illumina, Kinnex)",
        color = "K-value bin"
    ) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(args[[2]], plot = p, width = 12, height = 7)
