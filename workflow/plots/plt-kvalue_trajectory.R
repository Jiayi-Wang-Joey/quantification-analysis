suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

load_kvalue <- function(paths, kind) {
    sub_paths <- paths[grepl(sprintf("/metric/kvalue/%s/", kind), paths)]
    k <- rbindlist(lapply(sub_paths, fread, header = TRUE))
    k[, version := factor(version, levels = sort(unique(version)))]
    k[!is.na(`K-value`)]
}

make_trajectory_plot <- function(k, kind_label) {
    k[, gene_base := sub("\\..*$", "", Gene_id)]

    n_versions <- length(unique(k$version))
    gene_counts <- k[, .N, by = gene_base]
    common_genes <- gene_counts[N == n_versions, gene_base]
    k_common <- k[gene_base %in% common_genes]

    n_highlight <- 300
    swing <- k_common[, .(range = max(`K-value`) - min(`K-value`)), by = gene_base]
    top_genes <- swing[order(-range)][1:n_highlight, gene_base]
    k_highlight <- k_common[gene_base %in% top_genes]

    ggplot(k_common, aes(x = version, y = `K-value`)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.4) +
        geom_line(
            data = k_highlight,
            aes(group = gene_base),
            alpha = 0.15, color = "steelblue"
        ) +
        scale_y_log10() +
        labs(
            x = "GENCODE version",
            y = "K-value (log10 scale)",
            title = sprintf("K-value trajectories (%s)", kind_label)
        ) +
        theme_bw()
}

pdf(args[[2]], width = 10, height = 6)
print(make_trajectory_plot(load_kvalue(paths, "comprehensive"), "Comprehensive"))
print(make_trajectory_plot(load_kvalue(paths, "basic"), "Basic"))
invisible(dev.off())
