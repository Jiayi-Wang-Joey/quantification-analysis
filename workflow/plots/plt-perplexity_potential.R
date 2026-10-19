suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

oarfish_perplexity_paths <- paths[grepl("/metric/perplexity/oarfish/", paths)]
salmon_perplexity_paths <- paths[grepl("/metric/perplexity/salmon/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("gene_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("gene_counts\\.tsv$", paths)]

# gene TPM > 1 (no isoform filter)
load_perplexity <- function(perplexity_paths, count_paths, platform_label) {
    dt <- rbindlist(lapply(perplexity_paths, fread, header = TRUE))
    dt <- dt[, .(gene_id, gencode_version, gene_potential, gene_perplexity)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "gene_tpm")
    tpm <- tpm[, .(gene_id, gencode_version, gene_tpm)]

    merged <- merge(dt, tpm, by = c("gene_id", "gencode_version"))
    merged <- merged[gene_tpm > 1]
    merged[, platform := platform_label]
    merged[, .(gene_id, gencode_version, gene_potential, gene_perplexity, platform)]
}

dt <- rbindlist(list(
    load_perplexity(oarfish_perplexity_paths, oarfish_count_paths, "Kinnex (oarfish)"),
    load_perplexity(salmon_perplexity_paths, salmon_count_paths, "Illumina (salmon)")
))

plot_dt <- melt(
    dt, id.vars = c("gene_id", "gencode_version", "platform"),
    measure.vars = c("gene_potential", "gene_perplexity"),
    variable.name = "metric", value.name = "value"
)
plot_dt[, metric := factor(
    metric,
    levels = c("gene_potential", "gene_perplexity"),
    labels = c("Gene potential (isoforms w/ TPM>0)", "Gene perplexity (effective isoforms)")
)]
plot_dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]
plot_dt[, platform := factor(platform, levels = c("Kinnex (oarfish)", "Illumina (salmon)"))]

p <- ggplot(plot_dt, aes(x = gencode_version, y = value, color = metric)) +
    geom_boxplot(outlier.shape = NA, fill = "white", linewidth = 0.4) +
    scale_y_log10() +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~platform, ncol = 1) +
    labs(
        x = "GENCODE version",
        y = "Value per gene (log10 scale)",
        color = "Metric"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(args[[2]], plot = p, width = 10, height = 7)
