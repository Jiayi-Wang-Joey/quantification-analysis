suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


paths <- args[[1]]

kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]
oarfish_ambig_paths <- paths[grepl("/metric/ambiguity/oarfish/", paths) & grepl("gene_ambiguity\\.tsv$", paths)]
salmon_ambig_paths <- paths[grepl("/metric/ambiguity/salmon/", paths) & grepl("gene_ambiguity\\.tsv$", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("gene_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("gene_counts\\.tsv$", paths)]

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, gencode_version := as.integer(gencode_version)]
k <- k[, .(gene_id, gencode_version, `K-value`)]

oarfish <- rbindlist(lapply(oarfish_ambig_paths, fread, header = TRUE))
oarfish <- oarfish[, .(gene_id, gencode_version, ambiguous_fraction)]

salmon <- rbindlist(lapply(salmon_ambig_paths, fread, header = TRUE))
salmon <- salmon[, .(gene_id, gencode_version, ambiguous_fraction)]

oarfish_tpm <- rbindlist(lapply(oarfish_count_paths, fread, header = TRUE))
setnames(oarfish_tpm, "H69_tpm", "tpm")
oarfish_tpm <- oarfish_tpm[, .(gene_id, gencode_version, tpm)]

salmon_tpm <- rbindlist(lapply(salmon_count_paths, fread, header = TRUE))
setnames(salmon_tpm, "H69_tpm", "tpm")
salmon_tpm <- salmon_tpm[, .(gene_id, gencode_version, tpm)]

make_plot_dt <- function(dt, tpm_dt) {
    merged <- merge(dt, k, by = c("gene_id", "gencode_version"))
    merged <- merge(merged, tpm_dt, by = c("gene_id", "gencode_version"))
    merged <- merged[!is.na(ambiguous_fraction)]
    # TPM > 1 (miniQuant)
    merged <- merged[tpm > 1]
    merged[, gencode_version := factor(
        paste0("v", gencode_version),
        levels = paste0("v", sort(unique(gencode_version)))
    )]
    merged
}

oarfish_dt <- make_plot_dt(oarfish, oarfish_tpm)
salmon_dt <- make_plot_dt(salmon, salmon_tpm)

make_kvalue_ambiguity_plot <- function(dt, subtitle) {
    cor_dt <- dt[, .(
        rho = cor(`K-value`, ambiguous_fraction, method = "spearman"),
        x = min(`K-value`),
        y = max(ambiguous_fraction)
    ), by = gencode_version]
    cor_dt[, label := sprintf("rho = %.2f", rho)]

    ggplot(dt, aes(x = `K-value`, y = ambiguous_fraction, color = gencode_version)) +
        geom_point_rast(alpha = 0.4, size = 0.5) +
        geom_text(
            data = cor_dt,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
        ) +
        facet_wrap(~gencode_version, ncol = 4) +
        scale_x_log10() +
        scale_color_brewer(palette = "Paired") +
        labs(
            x = "K-value (log10 scale)",
            y = "Ambiguous fraction: ambig_reads / (ambig_reads + unique_reads)",
            subtitle = subtitle,
            color = "GENCODE version"
        ) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA))
}

pdf(args[[2]], width = 12, height = 6)
print(make_kvalue_ambiguity_plot(oarfish_dt, "Kinnex (oarfish)"))
print(make_kvalue_ambiguity_plot(salmon_dt, "Illumina (salmon)"))
invisible(dev.off())
