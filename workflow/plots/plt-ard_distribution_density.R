suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]
oarfish_tx_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("transcript_counts\\.tsv$", paths)]
salmon_tx_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("transcript_counts\\.tsv$", paths)]

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, gencode_version := as.integer(gencode_version)]
k <- k[, .(gene_id, gencode_version, `K-value`)]

compute_ard <- function(oarfish, salmon, join_cols) {
    merged <- merge(oarfish, salmon, by = join_cols)
    # TPM > 1 in both platforms (miniQuant)
    merged <- merged[illumina_tpm > 1 & kinnex_tpm > 1]
    merged[, ard := abs(illumina - kinnex) / (illumina + kinnex)]
    merged
}


oarfish_tx <- rbindlist(lapply(oarfish_tx_paths, fread, header = TRUE))
setnames(oarfish_tx, c("H69", "H69_tpm"), c("kinnex", "kinnex_tpm"))
oarfish_tx <- oarfish_tx[, .(transcript_id, gencode_version, kinnex, kinnex_tpm)]
salmon_tx <- rbindlist(lapply(salmon_tx_paths, fread, header = TRUE))
setnames(salmon_tx, c("H69", "H69_tpm"), c("illumina", "illumina_tpm"))
salmon_tx <- salmon_tx[, .(gene_id, transcript_id, gencode_version, illumina, illumina_tpm)]
tx_ard <- compute_ard(oarfish_tx, salmon_tx, c("transcript_id", "gencode_version"))
tx_ard <- merge(tx_ard, k, by = c("gene_id", "gencode_version"))

gene_ard <- tx_ard[, .(ard = mean(ard), `K-value` = `K-value`[1]), by = .(gene_id, gencode_version)]


n_top <- 300
top_genes <- k[order(-`K-value`), head(.SD, n_top), by = gencode_version][, .(gene_id, gencode_version)]

make_plot_dt <- function(dt, resolution_label) {
    all_dt <- copy(dt)
    all_dt[, gene_set := "All genes"]

    high_dt <- dt[`K-value` > 90]
    high_dt[, gene_set := "K-value > 90"]

    top_dt <- merge(dt, top_genes, by = c("gene_id", "gencode_version"))
    top_dt[, gene_set := "Top 300"]

    out <- rbindlist(list(all_dt, high_dt, top_dt), use.names = TRUE)
    out[, resolution := resolution_label]
    out[, .(ard, gencode_version, gene_set, resolution)]
}

plot_dt <- rbindlist(list(
    make_plot_dt(gene_ard, "Gene level"),
    make_plot_dt(tx_ard, "Transcript level")
))
plot_dt[, gencode_version := factor(
    paste0("v", gencode_version),
    levels = paste0("v", sort(unique(gencode_version)))
)]
plot_dt[, resolution := factor(resolution, levels = c("Gene level", "Transcript level"))]
plot_dt[, gene_set := factor(gene_set, levels = c("All genes", "K-value > 90", "Top 300"))]

p <- ggplot(plot_dt, aes(x = ard, color = gencode_version)) +
    geom_density() +
    scale_color_brewer(palette = "Paired") +
    facet_grid(resolution ~ gene_set) +
    labs(
        x = "ARD / MARD: |Illumina - Kinnex| / (Illumina + Kinnex)",
        y = "Density",
        color = "GENCODE version"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(args[[2]], plot = p, width = 8, height = 4.5)
