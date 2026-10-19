suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

oarfish_ambig_paths <- paths[grepl("/metric/ambiguity/oarfish/", paths) & grepl("gene_ambiguity\\.tsv$", paths)]
salmon_ambig_paths <- paths[grepl("/metric/ambiguity/salmon/", paths) & grepl("gene_ambiguity\\.tsv$", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("gene_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("gene_counts\\.tsv$", paths)]

load_ambig <- function(ambig_paths, count_paths, platform_label) {
    ambig <- rbindlist(lapply(ambig_paths, fread, header = TRUE))
    ambig <- ambig[, .(gene_id, gencode_version, ambiguous_fraction)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "tpm")
    tpm <- tpm[, .(gene_id, gencode_version, tpm)]

    merged <- merge(ambig, tpm, by = c("gene_id", "gencode_version"))
    # TPM > 1 (miniQuant)
    merged <- merged[!is.na(ambiguous_fraction) & tpm > 1]
    merged[, platform := platform_label]
    merged[, .(gene_id, gencode_version, ambiguous_fraction, platform)]
}

oarfish_dt <- load_ambig(oarfish_ambig_paths, oarfish_count_paths, "Kinnex (oarfish)")
salmon_dt <- load_ambig(salmon_ambig_paths, salmon_count_paths, "Illumina (salmon)")

plot_dt <- rbindlist(list(oarfish_dt, salmon_dt))
plot_dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]
plot_dt[, platform := factor(platform, levels = c("Kinnex (oarfish)", "Illumina (salmon)"))]

p <- ggplot(plot_dt, aes(x = gencode_version, y = ambiguous_fraction, fill = gencode_version)) +
    geom_violin(draw_quantiles = 0.5, linewidth = 0.4) +
    scale_fill_brewer(palette = "Paired") +
    facet_wrap(~platform, ncol = 1) +
    labs(
        x = "GENCODE version",
        y = "Ambiguous fraction: ambig_reads / (ambig_reads + unique_reads)",
        fill = "GENCODE version",
        subtitle = "Independent of K-value: does ambiguity rise with version, and is long-read less ambiguous? (line = median)"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(args[[2]], plot = p, width = 9, height = 7)
