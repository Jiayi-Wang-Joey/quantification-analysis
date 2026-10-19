suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


paths <- args[[1]]

kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]
oarfish_perplexity_paths <- paths[grepl("/metric/perplexity/oarfish/", paths)]
salmon_perplexity_paths <- paths[grepl("/metric/perplexity/salmon/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("gene_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("gene_counts\\.tsv$", paths)]

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, gencode_version := as.integer(gencode_version)]
k <- k[, .(gene_id, gencode_version, `K-value`)]

# gene TPM > 1 (no isoform filter)
load_perplexity <- function(perplexity_paths, count_paths, platform_label) {
    dt <- rbindlist(lapply(perplexity_paths, fread, header = TRUE))
    dt <- dt[, .(gene_id, gencode_version, gene_potential, gene_perplexity)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "gene_tpm")
    tpm <- tpm[, .(gene_id, gencode_version, gene_tpm)]

    merged <- merge(dt, tpm, by = c("gene_id", "gencode_version"))
    merged <- merge(merged, k, by = c("gene_id", "gencode_version"))
    merged <- merged[gene_tpm > 1]
    merged[, platform := platform_label]
    merged[, .(gene_id, gencode_version, gene_potential, gene_perplexity, `K-value`, platform)]
}

dt <- rbindlist(list(
    load_perplexity(oarfish_perplexity_paths, oarfish_count_paths, "Kinnex (oarfish)"),
    load_perplexity(salmon_perplexity_paths, salmon_count_paths, "Illumina (salmon)")
))
dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]
dt[, platform := factor(platform, levels = c("Kinnex (oarfish)", "Illumina (salmon)"))]

make_scatter_plot <- function(sub_dt, subtitle) {
    ggplot(sub_dt, aes(x = gene_potential, y = gene_perplexity, color = `K-value`)) +
        geom_point_rast(alpha = 0.4, size = 0.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
        scale_x_log10() +
        scale_y_log10() +
        scale_color_viridis_c(trans = "log10") +
        facet_wrap(~gencode_version, ncol = 4) +
        labs(
            x = "Gene potential: isoforms with TPM>0 (log10 scale)",
            y = "Gene perplexity: effective isoforms (log10 scale)",
            subtitle = subtitle
        ) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA))
}

pdf(args[[2]], width = 12, height = 7)
print(make_scatter_plot(dt[platform == "Illumina (salmon)"], "Illumina (salmon)"))
print(make_scatter_plot(dt[platform == "Kinnex (oarfish)"], "Kinnex (oarfish)"))
invisible(dev.off())
