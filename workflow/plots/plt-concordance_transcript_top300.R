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
oarfish <- oarfish[, .(transcript_id, gencode_version, kinnex, kinnex_tpm)]

salmon <- rbindlist(lapply(salmon_paths, fread, header = TRUE))
setnames(salmon, c("H69", "H69_tpm"), c("illumina", "illumina_tpm"))
salmon <- salmon[, .(gene_id, transcript_id, gencode_version, illumina, illumina_tpm)]

# top 300 genes by K-value per version, all their transcripts
n_top <- 300
top_genes <- k[order(-`K-value`), head(.SD, n_top), by = gencode_version]

merged <- merge(oarfish, salmon, by = c("transcript_id", "gencode_version"))
merged <- merge(merged, top_genes, by = c("gene_id", "gencode_version"))

# TPM > 1 in both platforms (miniQuant)
merged <- merged[illumina_tpm > 1 & kinnex_tpm > 1]
merged[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]

cor_dt <- merged[, .(
    rho = cor(illumina, kinnex, method = "spearman"),
    x = min(illumina),
    y = max(kinnex)
), by = gencode_version]
cor_dt[, label := sprintf("rho = %.2f", rho)]

p <- ggplot(merged, aes(x = illumina, y = kinnex, color = `K-value`)) +
    geom_point_rast(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(
        data = cor_dt,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
    ) +
    facet_wrap(~gencode_version, ncol = 4) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_viridis_c(trans = "log10") +
    coord_equal() +
    labs(
        x = "Illumina counts (salmon)",
        y = "Kinnex counts (oarfish)",
        title = sprintf("Illumina vs Kinnex transcript counts (top %d K-value genes)", n_top)
    ) +
    theme_classic()

ggsave(args[[2]], plot = p, width = 9, height = 6)
