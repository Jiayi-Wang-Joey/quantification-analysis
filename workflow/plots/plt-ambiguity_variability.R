suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

oarfish_ambig_paths <- paths[grepl("/metric/ambiguity/oarfish/", paths) & grepl("gene_ambiguity\\.tsv$", paths)]
salmon_ambig_paths <- paths[grepl("/metric/ambiguity/salmon/", paths) & grepl("gene_ambiguity\\.tsv$", paths)]
oarfish_count_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("gene_counts\\.tsv$", paths)]
salmon_count_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("gene_counts\\.tsv$", paths)]
oarfish_tx_paths <- paths[grepl("/aggregate/oarfish/", paths) & grepl("transcript_counts\\.tsv$", paths)]
salmon_tx_paths <- paths[grepl("/aggregate/salmon/", paths) & grepl("transcript_counts\\.tsv$", paths)]
kvalue_paths <- paths[grepl("/metric/kvalue/", paths)]

MIN_DEPTH <- 20
MIN_VERSIONS <- 5
N_EXAMPLES <- 3

k <- rbindlist(lapply(kvalue_paths, fread, header = TRUE))
k <- k[!is.na(`K-value`)]
setnames(k, c("Gene_id", "version"), c("gene_id", "gencode_version"))
k[, base_gene_id := sub("\\.[0-9]+$", "", gene_id)]
k[, gencode_version := as.integer(gencode_version)]

load_ambig <- function(gene_paths, count_paths) {
    dt <- rbindlist(lapply(gene_paths, fread, header = TRUE))
    dt[, base_gene_id := sub("\\.[0-9]+$", "", gene_id)]
    dt[, total_reads := unique_reads + ambig_reads]

    tpm <- rbindlist(lapply(count_paths, fread, header = TRUE))
    setnames(tpm, "H69_tpm", "tpm")
    tpm <- tpm[, .(gene_id, gencode_version, tpm)]

    merge(dt, tpm, by = c("gene_id", "gencode_version"))
}

oarfish_ambig <- load_ambig(oarfish_ambig_paths, oarfish_count_paths)
salmon_ambig <- load_ambig(salmon_ambig_paths, salmon_count_paths)

compute_variability <- function(dt) {
    # TPM > 1 (miniQuant)
    d <- dt[total_reads >= MIN_DEPTH & tpm > 1]
    out <- d[, .(
        sd_ambig = sd(ambiguous_fraction),
        n_versions = .N,
        gene_name = gene_name[1]
    ), by = base_gene_id]
    out[n_versions >= MIN_VERSIONS]
}

oarfish_var <- compute_variability(oarfish_ambig)
salmon_var <- compute_variability(salmon_ambig)

combined <- merge(
    oarfish_var[, .(base_gene_id, gene_name, sd_oarfish = sd_ambig)],
    salmon_var[, .(base_gene_id, sd_salmon = sd_ambig)],
    by = "base_gene_id", all = TRUE
)
combined[, sd_max := pmax(sd_oarfish, sd_salmon, na.rm = TRUE)]
combined <- combined[order(-sd_max)]

top_genes <- head(combined$base_gene_id, N_EXAMPLES)
top_gene_names <- head(combined$gene_name, N_EXAMPLES)

load_tx <- function(tx_paths, platform_label, target_gene_id) {
    dt <- rbindlist(lapply(tx_paths, fread, header = TRUE))
    dt[, base_gene_id := sub("\\.[0-9]+$", "", gene_id)]
    dt <- dt[base_gene_id == target_gene_id]
    setnames(dt, "H69_tpm", "tpm")
    dt[, platform := platform_label]
    dt[, .(transcript_id, gencode_version, tpm, platform)]
}

make_gene_plot <- function(target_gene_id, target_gene_name) {
    oarfish_tx <- load_tx(oarfish_tx_paths, "Kinnex (oarfish)", target_gene_id)
    salmon_tx <- load_tx(salmon_tx_paths, "Illumina (salmon)", target_gene_id)

    plot_dt <- rbindlist(list(oarfish_tx, salmon_tx))
    plot_dt[, gencode_version := factor(
        paste0("v", gencode_version),
        levels = paste0("v", sort(unique(gencode_version)))
    )]

    gene_k <- k[base_gene_id == target_gene_id]
    kvalue_str <- if (nrow(gene_k) > 0) {
        latest <- gene_k[which.max(gencode_version)]
        sprintf(", K-value = %.1f (v%d)", latest$`K-value`, latest$gencode_version)
    } else {
        ""
    }

    ggplot(plot_dt, aes(x = transcript_id, y = tpm, fill = platform)) +
        geom_col(position = "dodge") +
        facet_wrap(~gencode_version, scales = "free_x", ncol = 4) +
        labs(
            x = "Transcript",
            y = "TPM",
            fill = "Platform",
            subtitle = sprintf("%s (%s)%s", target_gene_name, target_gene_id, kvalue_str)
        ) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            panel.border = element_rect(color = "black", fill = NA)
        )
}

pdf(args[[2]], width = 14, height = 8)
for (i in seq_along(top_genes)) {
    print(make_gene_plot(top_genes[i], top_gene_names[i]))
}
invisible(dev.off())
