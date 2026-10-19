suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


paths <- args[[1]]

parse_stats <- function(path) {
    lines <- readLines(path)
    sn <- lines[startsWith(lines, "SN")]
    get_val <- function(key) {
        line <- sn[grepl(paste0("^SN\t", key, ":"), sn)]
        as.numeric(sub(paste0("^SN\t", key, ":\t"), "", line[1]))
    }
    version <- sub(".*/v([0-9]+)/.*", "\\1", path)
    platform <- if (grepl("/kinnex/", path)) "Kinnex" else "Illumina"
    total <- get_val("raw total sequences")
    mapped <- get_val("reads mapped")
    data.table(
        version = version,
        platform = platform,
        total_reads = total / 1e6,
        mapped_reads = mapped / 1e6,
        mapping_rate = 100 * mapped / total,
        avg_length = get_val("average length"),
        error_rate = 100 * get_val("error rate")
    )
}

dt <- rbindlist(lapply(paths, parse_stats))
dt[, version := factor(
    paste0("v", version),
    levels = paste0("v", sort(unique(as.integer(version))))
)]

long_dt <- melt(
    dt,
    id.vars = c("version", "platform"),
    measure.vars = c("total_reads", "mapped_reads", "mapping_rate", "avg_length", "error_rate"),
    variable.name = "metric",
    value.name = "value"
)

metric_labels <- c(
    total_reads = "Total reads (M)",
    mapped_reads = "Mapped reads (M)",
    mapping_rate = "Mapping rate (%)",
    avg_length = "Average Read length",
    error_rate = "Error rate (%)"
)
long_dt[, metric := factor(metric_labels[as.character(metric)], levels = metric_labels)]

p <- ggplot(long_dt, aes(x = platform, y = value, fill = platform)) +
    geom_col() +
    facet_grid(version ~ metric, scales = "free_y") +
    labs(x = NULL, y = "Value", fill = "Technology", title = "Data quality summary") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(args[[2]], plot = p, width = 12, height = 10)
