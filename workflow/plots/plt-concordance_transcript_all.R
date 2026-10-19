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

# CPM on the full library, before filtering
oarfish <- rbindlist(lapply(oarfish_paths, fread, header = TRUE))
setnames(oarfish, c("H69", "H69_tpm"), c("kinnex", "kinnex_tpm"))
oarfish[, kinnex_cpm := kinnex / sum(kinnex) * 1e6, by = gencode_version]
oarfish <- oarfish[, .(transcript_id, gencode_version, kinnex, kinnex_cpm, kinnex_tpm)]

salmon <- rbindlist(lapply(salmon_paths, fread, header = TRUE))
setnames(salmon, c("H69", "H69_tpm"), c("illumina", "illumina_tpm"))
salmon[, illumina_cpm := illumina / sum(illumina) * 1e6, by = gencode_version]
salmon <- salmon[, .(gene_id, transcript_id, gencode_version, illumina, illumina_cpm, illumina_tpm)]

merged <- merge(oarfish, salmon, by = c("transcript_id", "gencode_version"))
merged <- merge(merged, k, by = c("gene_id", "gencode_version"))

# TPM > 1 in both platforms (miniQuant)
merged <- merged[illumina_tpm > 1 & kinnex_tpm > 1]
merged[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]

raw_dt <- merged[, .(transcript_id, gencode_version, illumina, kinnex)]
cpm_dt <- merged[, .(transcript_id, gencode_version, illumina = illumina_cpm, kinnex = kinnex_cpm)]

make_scatter_plot <- function(dt, subtitle) {
    cor_dt <- dt[, .(
        rho = cor(illumina, kinnex, method = "spearman"),
        x = min(illumina),
        y = max(kinnex)
    ), by = gencode_version]
    cor_dt[, label := sprintf("rho = %.2f", rho)]

    ggplot(dt, aes(x = illumina, y = kinnex)) +
        geom_point_rast(color = "black", size = 0.3, alpha = 0.3) +
        geom_density_2d(color = "lightblue") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        geom_text(
            data = cor_dt,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
        ) +
        facet_wrap(~gencode_version, ncol = 4) +
        scale_x_log10() +
        scale_y_log10() +
        coord_equal() +
        labs(
            x = "Illumina (salmon)",
            y = "Kinnex (oarfish)",
            subtitle = subtitle
        ) +
        theme_classic()
}

pdf(args[[2]], width = 12, height = 6)
print(make_scatter_plot(raw_dt, "Raw counts, all transcripts"))
print(make_scatter_plot(cpm_dt, "CPM, all transcripts"))
invisible(dev.off())
