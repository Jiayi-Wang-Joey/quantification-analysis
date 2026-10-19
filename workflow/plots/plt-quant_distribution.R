suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]
res <- lapply(paths, function(p) {
    dt <- fread(p, header = TRUE)
    dt[, version := sub(".*\\.v([0-9]+)\\.quant$", "\\1", p)]
    dt
})
q <- rbindlist(res)
q[, version := factor(version, levels = sort(unique(as.integer(version))))]
q <- q[num_reads > 0]

p <- ggplot(q, aes(x = version, y = num_reads, col = version)) +
    geom_violin() +
    scale_color_brewer(palette = "Paired") +
    scale_y_log10() +
    labs(
        x = "GENCODE version",
        y = "Number of reads (log10 scale)",
        title = "Oarfish read counts"
    ) +
    theme_bw()

ggsave(args[[2]], plot = p, width = 8, height = 5)
