suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

dt <- fread("gencode_version_commit_mentions.tsv", header = TRUE)
dt[, date := as.IDate(as.POSIXct(date, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"))]
dt[, month := as.IDate(format(date, "%Y-%m-01"))]

release_dates <- data.table(
    version = c(25, 30, 35, 39, 44, 45, 46, 47, 48, 49, 50),
    release_date = as.IDate(c(
        "2016-07-01", "2019-04-01", "2020-08-01", "2021-12-01", "2023-07-01",
        "2024-01-01", "2024-05-01", "2024-10-01", "2025-05-01", "2025-09-01", "2026-06-01"
    ))
)

monthly <- dt[, .N, by = .(version, month)]
monthly[, version := factor(version, levels = release_dates$version)]

# "Paired" palette without the pale yellow
paired_colors <- RColorBrewer::brewer.pal(12, "Paired")[-11]

p <- ggplot(monthly, aes(x = month, y = N, color = version)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    geom_vline(
        data = release_dates,
        aes(xintercept = release_date, color = factor(version, levels = release_dates$version)),
        linetype = "dashed", linewidth = 0.5, alpha = 0.6, show.legend = FALSE
    ) +
    scale_color_manual(values = paired_colors) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
        x = "Month",
        y = "GitHub commits mentioning this GENCODE version",
        color = "GENCODE version",
        title = "Real-world adoption of GENCODE versions over time",
        subtitle = "GitHub commit-search mentions as a usage proxy; dashed lines = each version's actual release date"
    ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA))

ggsave("gencode_usage_over_time.pdf", plot = p, width = 12, height = 6)
