suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


paths <- args[[1]]

oarfish_infrv_paths <- paths[grepl("/metric/infrv/oarfish/", paths)]
salmon_infrv_paths <- paths[grepl("/metric/infrv/salmon/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("transcript_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("transcript_counts\\.tsv$", paths)]

load_infrv <- function(infrv_paths, count_paths, colname) {
    infrv <- rbindlist(lapply(infrv_paths, fread, header = TRUE))
    infrv <- infrv[, .(transcript_id, gencode_version, InfRV)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "tpm")
    tpm <- tpm[, .(transcript_id, gencode_version, tpm)]

    merged <- merge(infrv, tpm, by = c("transcript_id", "gencode_version"))
    # TPM > 1 (miniQuant)
    merged <- merged[!is.na(InfRV) & tpm > 1, .(transcript_id, gencode_version, InfRV)]
    setnames(merged, "InfRV", colname)
    merged
}

oarfish_dt <- load_infrv(oarfish_infrv_paths, oarfish_count_paths, "InfRV_oarfish")
salmon_dt <- load_infrv(salmon_infrv_paths, salmon_count_paths, "InfRV_salmon")

# transcripts with TPM > 1 in both platforms
plot_dt <- merge(oarfish_dt, salmon_dt, by = c("transcript_id", "gencode_version"))
plot_dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]

cor_dt <- plot_dt[, .(
    rho = cor(InfRV_salmon, InfRV_oarfish, method = "spearman"),
    x = min(InfRV_salmon),
    y = max(InfRV_oarfish)
), by = gencode_version]
cor_dt[, label := sprintf("rho = %.2f", rho)]

p <- ggplot(plot_dt, aes(x = InfRV_salmon, y = InfRV_oarfish)) +
    geom_point_rast(alpha = 0.3, size = 0.5, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_text(
        data = cor_dt,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
    ) +
    facet_wrap(~gencode_version, ncol = 4) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
        x = "InfRV, Illumina (salmon), log10 scale",
        y = "InfRV, Kinnex (oarfish), log10 scale"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(args[[2]], plot = p, width = 12, height = 7)
