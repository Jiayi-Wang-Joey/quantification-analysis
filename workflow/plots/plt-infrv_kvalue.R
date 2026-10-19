suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


paths <- args[[1]]

kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]
oarfish_infrv_paths <- paths[grepl("/metric/infrv/oarfish/", paths)]
salmon_infrv_paths <- paths[grepl("/metric/infrv/salmon/", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("transcript_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("transcript_counts\\.tsv$", paths)]

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, gencode_version := as.integer(gencode_version)]
k <- k[, .(gene_id, gencode_version, `K-value`)]

oarfish_infrv <- rbindlist(lapply(oarfish_infrv_paths, fread, header = TRUE))
oarfish_infrv <- oarfish_infrv[, .(transcript_id, gencode_version, InfRV)]

salmon_infrv <- rbindlist(lapply(salmon_infrv_paths, fread, header = TRUE))
salmon_infrv <- salmon_infrv[, .(transcript_id, gencode_version, InfRV)]

oarfish_tx <- rbindlist(lapply(oarfish_count_paths, fread, header = TRUE))
setnames(oarfish_tx, "H69_tpm", "tpm")
oarfish_tx <- oarfish_tx[, .(gene_id, transcript_id, gencode_version, tpm)]

salmon_tx <- rbindlist(lapply(salmon_count_paths, fread, header = TRUE))
setnames(salmon_tx, "H69_tpm", "tpm")
salmon_tx <- salmon_tx[, .(gene_id, transcript_id, gencode_version, tpm)]

# gene K-value joined to each transcript
make_plot_dt <- function(infrv_dt, tx_dt) {
    merged <- merge(infrv_dt, tx_dt, by = c("transcript_id", "gencode_version"))
    merged <- merge(merged, k, by = c("gene_id", "gencode_version"))
    # TPM > 1 (miniQuant)
    merged <- merged[!is.na(InfRV) & tpm > 1]
    merged[, gencode_version := factor(
        paste0("v", gencode_version),
        levels = paste0("v", sort(unique(gencode_version)))
    )]
    merged
}

oarfish_dt <- make_plot_dt(oarfish_infrv, oarfish_tx)
salmon_dt <- make_plot_dt(salmon_infrv, salmon_tx)

make_kvalue_infrv_plot <- function(dt, subtitle) {
    cor_dt <- dt[, .(
        rho = cor(`K-value`, InfRV, method = "spearman"),
        x = min(`K-value`),
        y = max(InfRV)
    ), by = gencode_version]
    cor_dt[, label := sprintf("rho = %.2f", rho)]

    ggplot(dt, aes(x = `K-value`, y = InfRV, color = gencode_version)) +
        geom_point_rast(alpha = 0.4, size = 0.5) +
        geom_text(
            data = cor_dt,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
        ) +
        facet_wrap(~gencode_version, ncol = 4) +
        scale_x_log10() +
        scale_y_log10() +
        scale_color_brewer(palette = "Paired") +
        labs(
            x = "K-value (log10 scale)",
            y = "InfRV (fishpond::computeInfRV formula, log10 scale)",
            subtitle = subtitle,
            color = "GENCODE version"
        ) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA))
}

pdf(args[[2]], width = 12, height = 6)
print(make_kvalue_infrv_plot(oarfish_dt, "Kinnex (oarfish)"))
print(make_kvalue_infrv_plot(salmon_dt, "Illumina (salmon)"))
invisible(dev.off())
