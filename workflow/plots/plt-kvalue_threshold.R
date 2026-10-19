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

make_threshold_plot <- function(k, threshold, kind_label) {
    summary_dt <- k[, .(
        total = .N,
        n_high = sum(`K-value` > threshold),
        pct_high = 100 * sum(`K-value` > threshold) / .N
    ), by = version]

    ggplot(summary_dt, aes(x = version, y = pct_high)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = sprintf("%.1f%%", pct_high)), vjust = -0.5, size = 3) +
        labs(
            x = "GENCODE version",
            y = sprintf("%% of genes with K-value > %d", threshold),
            title = sprintf("Genes with K-value > %d (%s)", threshold, kind_label)
        ) +
        theme_bw() +
        ylim(0, max(summary_dt$pct_high) * 1.15)
}

pdf(args[[2]], width = 8, height = 5)
for (kind in list(list(dir = "comprehensive", label = "Comprehensive"), list(dir = "basic", label = "Basic"))) {
    k <- load_kvalue(paths, kind$dir)
    print(make_threshold_plot(k, 25, kind$label))
    print(make_threshold_plot(k, 90, kind$label))
}
invisible(dev.off())
