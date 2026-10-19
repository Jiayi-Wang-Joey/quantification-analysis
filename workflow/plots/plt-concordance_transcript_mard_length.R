suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("transcript_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("transcript_counts\\.tsv$", paths)]
oarfish_quant_paths <- paths[grepl("/quantify/run_oarfish/", paths) & grepl("\\.quant$", paths)]

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, gencode_version := as.integer(gencode_version)]
k <- k[, .(gene_id, gencode_version, `K-value`)]

# transcript length from the oarfish .quant files
length_dt <- rbindlist(lapply(oarfish_quant_paths, function(p) {
    dt <- fread(p, header = TRUE)
    setnames(dt, c("tname", "len"), c("transcript_id", "length"))
    dt[, gencode_version := as.integer(sub(".*\\.v([0-9]+)\\.quant$", "\\1", p))]
    dt[, .(transcript_id, gencode_version, length)]
}))

oarfish <- rbindlist(lapply(oarfish_count_paths, fread, header = TRUE))
setnames(oarfish, c("H69", "H69_tpm"), c("kinnex", "kinnex_tpm"))
oarfish <- oarfish[, .(gene_id, transcript_id, gencode_version, kinnex, kinnex_tpm)]

salmon <- rbindlist(lapply(salmon_count_paths, fread, header = TRUE))
setnames(salmon, c("H69", "H69_tpm"), c("illumina", "illumina_tpm"))
salmon <- salmon[, .(transcript_id, gencode_version, illumina, illumina_tpm)]

merged <- merge(oarfish, salmon, by = c("transcript_id", "gencode_version"))
merged <- merge(merged, k, by = c("gene_id", "gencode_version"))
merged <- merge(merged, length_dt, by = c("transcript_id", "gencode_version"))

# TPM > 1 in both platforms (miniQuant)
merged <- merged[illumina_tpm > 1 & kinnex_tpm > 1]

# transcript-level ARD
merged[, mard := abs(illumina - kinnex) / ((illumina + kinnex) / 2)]

# miniQuant K-value bins, K > 25 split at 90
KVALUE_LABELS <- c("=1", "1-2", "2-3", "3-5", "5-8", "8-13", "13-25", "25-90", ">=90")
merged[, kvalue_bin := fcase(
    `K-value` == 1, "=1",
    `K-value` > 1 & `K-value` <= 2, "1-2",
    `K-value` > 2 & `K-value` <= 3, "2-3",
    `K-value` > 3 & `K-value` <= 5, "3-5",
    `K-value` > 5 & `K-value` <= 8, "5-8",
    `K-value` > 8 & `K-value` <= 13, "8-13",
    `K-value` > 13 & `K-value` <= 25, "13-25",
    `K-value` > 25 & `K-value` < 90, "25-90",
    `K-value` >= 90, ">=90"
)]
merged[, kvalue_bin := factor(kvalue_bin, levels = KVALUE_LABELS)]

# Kinnex short-transcript cut-off
merged[, length_bin := cut(
    length, breaks = c(0, 1250, Inf), labels = c("short (<1.25kb)", "long (>=1.25kb)")
)]

merged[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]

p <- ggplot(merged, aes(x = kvalue_bin, y = mard, color = kvalue_bin)) +
    geom_violin(fill = "white", linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.3) +
    scale_y_log10() +
    scale_color_brewer(palette = "Paired") +
    facet_grid(length_bin ~ gencode_version) +
    labs(
        x = "K-value bin",
        y = "Transcript ARD: |Illumina - Kinnex| / mean(Illumina, Kinnex)",
        color = "K-value bin"
    ) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
    )

ggsave(args[[2]], plot = p, width = 18, height = 7)
