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

make_distribution_plot <- function(k, kind_label) {
    ggplot(k, aes(x = version, y = `K-value`, col = version)) +
        geom_violin() +
        scale_color_brewer(palette = "Paired") +
        scale_y_log10() +
        labs(
            x = "GENCODE version",
            y = "K-value (log10 scale)",
            title = sprintf("K-value (%s)", kind_label)
        ) +
        theme_bw()
}

pdf(args[[2]], width = 8, height = 5)
print(make_distribution_plot(load_kvalue(paths, "comprehensive"), "Comprehensive"))
print(make_distribution_plot(load_kvalue(paths, "basic"), "Basic"))
invisible(dev.off())
