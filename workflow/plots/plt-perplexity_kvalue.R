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

# gene TPM > 1
make_plot_dt <- function(perplexity_paths, count_paths) {
    dt <- rbindlist(lapply(perplexity_paths, fread, header = TRUE))
    dt <- dt[, .(gene_id, gencode_version, gene_perplexity)]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "gene_tpm")
    tpm <- tpm[, .(gene_id, gencode_version, gene_tpm)]

    merged <- merge(dt, tpm, by = c("gene_id", "gencode_version"))
    merged <- merge(merged, k, by = c("gene_id", "gencode_version"))
    merged <- merged[gene_tpm > 1]
    merged[, gencode_version := factor(
        paste0("v", gencode_version),
        levels = paste0("v", sort(unique(gencode_version)))
    )]
    merged
}

oarfish_dt <- make_plot_dt(oarfish_perplexity_paths, oarfish_count_paths)
salmon_dt <- make_plot_dt(salmon_perplexity_paths, salmon_count_paths)

make_kvalue_perplexity_plot <- function(dt, subtitle) {
    cor_dt <- dt[, .(
        rho = cor(`K-value`, gene_perplexity, method = "spearman"),
        x = min(`K-value`),
        y = max(gene_perplexity)
    ), by = gencode_version]
    cor_dt[, label := sprintf("rho = %.2f", rho)]

    ggplot(dt, aes(x = `K-value`, y = gene_perplexity)) +
        geom_point_rast(alpha = 0.15, size = 0.4, color = "steelblue") +
        geom_text(
            data = cor_dt,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
        ) +
        facet_wrap(~gencode_version, ncol = 4) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
            x = "K-value (log10 scale)",
            y = "Gene perplexity (log10 scale)",
            subtitle = subtitle
        ) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA))
}

pdf(args[[2]], width = 12, height = 6)
print(make_kvalue_perplexity_plot(oarfish_dt, "Kinnex (oarfish)"))
print(make_kvalue_perplexity_plot(salmon_dt, "Illumina (salmon)"))
invisible(dev.off())
