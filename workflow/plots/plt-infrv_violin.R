suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

oarfish_infrv_paths <- paths[grepl("/metric/infrv/oarfish/", paths)]
salmon_infrv_paths <- paths[grepl("/metric/infrv/salmon/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("transcript_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("transcript_counts\\.tsv$", paths)]

# drop the top 0.5% so extreme InfRV doesn't flatten the violins
INFRV_DISPLAY_CAP <- 200

load_infrv <- function(infrv_paths, count_paths, platform_label) {
    infrv <- rbindlist(lapply(infrv_paths, fread, header = TRUE))
    infrv <- infrv[, .(transcript_id, gencode_version, InfRV)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "tpm")
    tpm <- tpm[, .(transcript_id, gencode_version, tpm)]

    merged <- merge(infrv, tpm, by = c("transcript_id", "gencode_version"))
    # TPM > 1 (miniQuant)
    merged <- merged[!is.na(InfRV) & tpm > 1 & InfRV <= INFRV_DISPLAY_CAP]
    merged[, platform := platform_label]
    merged[, .(transcript_id, gencode_version, InfRV, platform)]
}

oarfish_dt <- load_infrv(oarfish_infrv_paths, oarfish_count_paths, "Kinnex (oarfish)")
salmon_dt <- load_infrv(salmon_infrv_paths, salmon_count_paths, "Illumina (salmon)")

plot_dt <- rbindlist(list(oarfish_dt, salmon_dt))
plot_dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]
plot_dt[, platform := factor(platform, levels = c("Kinnex (oarfish)", "Illumina (salmon)"))]

p <- ggplot(plot_dt, aes(x = gencode_version, y = InfRV, color = gencode_version)) +
    geom_violin(fill = "white", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
    scale_y_log10() +
    scale_color_brewer(palette = "Paired") +
    facet_wrap(~platform, ncol = 1) +
    labs(
        x = "GENCODE version",
        y = "InfRV (fishpond::computeInfRV formula, log10 scale)",
        color = "GENCODE version"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(args[[2]], plot = p, width = 9, height = 7)
