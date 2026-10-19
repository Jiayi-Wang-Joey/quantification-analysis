suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


paths <- args[[1]]

oarfish_perplexity_paths <- paths[grepl("/metric/perplexity/oarfish/", paths)]
salmon_perplexity_paths <- paths[grepl("/metric/perplexity/salmon/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("gene_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("gene_counts\\.tsv$", paths)]

# gene TPM > 1
load_perplexity <- function(perplexity_paths, count_paths, colname) {
    dt <- rbindlist(lapply(perplexity_paths, fread, header = TRUE))
    dt <- dt[, .(gene_id, gencode_version, gene_perplexity)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "gene_tpm")
    tpm <- tpm[, .(gene_id, gencode_version, gene_tpm)]

    merged <- merge(dt, tpm, by = c("gene_id", "gencode_version"))
    merged <- merged[gene_tpm > 1, .(gene_id, gencode_version, gene_perplexity)]
    setnames(merged, "gene_perplexity", colname)
    merged
}

oarfish_dt <- load_perplexity(oarfish_perplexity_paths, oarfish_count_paths, "perplexity_oarfish")
salmon_dt <- load_perplexity(salmon_perplexity_paths, salmon_count_paths, "perplexity_salmon")

# genes with TPM > 1 in both platforms
plot_dt <- merge(oarfish_dt, salmon_dt, by = c("gene_id", "gencode_version"))
plot_dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]

cor_dt <- plot_dt[, .(
    rho = cor(perplexity_salmon, perplexity_oarfish, method = "spearman"),
    x = min(perplexity_salmon),
    y = max(perplexity_oarfish)
), by = gencode_version]
cor_dt[, label := sprintf("rho = %.2f", rho)]

p <- ggplot(plot_dt, aes(x = perplexity_salmon, y = perplexity_oarfish)) +
    geom_point_rast(alpha = 0.15, size = 0.4, color = "steelblue") +
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
        x = "Gene perplexity, Illumina (salmon), log10 scale",
        y = "Gene perplexity, Kinnex (oarfish), log10 scale"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(args[[2]], plot = p, width = 12, height = 7)
