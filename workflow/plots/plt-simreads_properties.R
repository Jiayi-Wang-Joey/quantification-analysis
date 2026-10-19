suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


# simulated reads vs real H69 per GENCODE version: Kinnex (pages 2-10), Illumina (11-16), both (17-18)

paths <- args[[1]]
GREY <- "#8f8e88"
VCOL <- c("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                           panel.grid.major = element_line(linewidth = 0.2),
                                           strip.background = element_blank(),
                                           legend.position = "top",
                                           plot.subtitle = element_text(size = 9)))
version_of <- function(f) paste0("v", sub(".*/v([0-9]+)/.*", "\\1", f))
is_real <- function(f) grepl("/H69[./]", f)
data_label <- function(f) ifelse(is_real(f), "real H69", paste("simulated", version_of(f)))
depth_num <- function(x) as.numeric(sub("M$", "", x)) * ifelse(grepl("M$", x), 1e6, 1)
versions <- sort(unique(version_of(grep("/v[0-9]+/", paths, value = TRUE))), decreasing = TRUE)
LEVELS <- c("real H69", paste("simulated", versions))
COLOURS <- setNames(c(GREY, VCOL[seq_along(versions)]), LEVELS)
lvl <- function(x) factor(x, levels = LEVELS)
colour_scale <- scale_colour_manual(values = COLOURS, name = NULL)
fill_scale <- scale_fill_manual(values = COLOURS, name = NULL)
dodge <- position_dodge(width = 0.7)

read_sn <- function(f) {
    sn <- fread(cmd = paste("grep '^SN'", f, "| cut -f2,3"), header = FALSE, sep = "\t", col.names = c("key", "value"))
    setNames(as.numeric(sn$value), sub(":$", "", sn$key))
}
read_section <- function(f, tag) {
    d <- fread(cmd = paste0("grep '^", tag, "\t' ", f, " | cut -f2-"), header = FALSE, sep = "\t", fill = TRUE)
    d[, `:=`(data = data_label(f), version = version_of(f))]
}

kstats <- grep("data_quality/kinnex/", paths, value = TRUE)
istats <- grep("read_qc/illumina/.*stats\\.txt$", paths, value = TRUE)


# Kinnex data
qc <- rbindlist(lapply(grep("read_qc/kinnex/", paths, value = TRUE), function(f) {
    fread(f)[, `:=`(data = data_label(f), version = version_of(f))]
}))
# position 0 = transcript 5' end on either strand
qc[, `:=`(coverage = (ref_end - ref_start) / transcript_length,
          five_missing = ref_start / transcript_length,
          three_missing = (transcript_length - ref_end) / transcript_length,
          identity = 1 - nm / (aligned_bases + inserted + deleted),
          mismatch_rate = (nm - inserted - deleted) / (aligned_bases + inserted + deleted),
          insertion_rate = inserted / (aligned_bases + inserted + deleted),
          deletion_rate = deleted / (aligned_bases + inserted + deleted))]
qc[, tx_bin := cut(transcript_length, c(0, 500, 1000, 2000, 3000, 5000, Inf),
                   labels = c("< 0.5 kb", "0.5-1 kb", "1-2 kb", "2-3 kb", "3-5 kb", ">= 5 kb"))]

ks <- rbindlist(lapply(kstats, function(f) {
    v <- read_sn(f)
    data.table(data = data_label(f), version = version_of(f),
               `% reads mapped` = 100 * v[["reads mapped"]] / v[["raw total sequences"]],
               `error rate (%)` = 100 * v[["error rate"]],
               `secondary alignments per mapped read` = v[["non-primary alignments"]] / v[["reads mapped"]],
               `split alignments (% of mapped)` = 100 * v[["supplementary alignments"]] / v[["reads mapped"]])
}))

rl <- rbindlist(lapply(kstats, function(f) setnames(read_section(f, "RL")[, 1:2], c("length", "n"))[
    , `:=`(data = data_label(f), version = version_of(f))]))
rl <- rl[, .(n = sum(n)), by = .(data, version, length)]
rl_sum <- rl[, .(mean = sum(length * n) / sum(n), short = 100 * sum(n[length < 300]) / sum(n)), by = .(data, version)]


# Illumina data
is_ <- rbindlist(lapply(istats, function(f) {
    v <- read_sn(f)
    data.table(data = data_label(f), version = version_of(f),
               `% reads mapped` = 100 * v[["reads mapped"]] / v[["raw total sequences"]],
               `error rate (%)` = 100 * v[["error rate"]], `read length` = v[["average length"]],
               `mean insert size` = v[["insert size average"]], `mean base quality` = v[["average quality"]])
}))
ins <- rbindlist(lapply(istats, function(f) setnames(read_section(f, "IS")[, 1:2], c("insert", "pairs"))[, `:=`(data = data_label(f), version = version_of(f))]))
gcf <- rbindlist(lapply(istats, function(f) setnames(read_section(f, "GCF")[, 1:2], c("gc", "n"))[
    , `:=`(data = data_label(f), version = version_of(f))]))
per_cycle <- function(f, tag, fun) {
    d <- fread(cmd = paste0("grep '^", tag, "\t' ", f, " | cut -f2-"), header = FALSE, sep = "\t", fill = TRUE)
    if (!nrow(d)) return(NULL)
    m <- as.matrix(d[, -1]); m[is.na(m)] <- 0
    data.table(cycle = d$V1, value = fun(m), data = data_label(f), version = version_of(f))
}
ffq <- rbindlist(lapply(istats, per_cycle, "FFQ", function(m) as.vector(m %*% (seq_len(ncol(m)) - 1)) / rowSums(m)))
mpc <- rbindlist(lapply(istats, per_cycle, "MPC", function(m) rowSums(m[, -1, drop = FALSE])))  # without N
reads_n <- is_[, .(data, version)][, reads := sapply(istats, function(f) read_sn(f)[["reads mapped"]])]
pos <- rbindlist(lapply(grep("read_qc/illumina/.*positions\\.tsv$", paths, value = TRUE), function(f) {
    fread(f)[, `:=`(data = data_label(f), version = version_of(f))]
}))


# 1. overview table
kq <- qc[, .(`polyA tail (mean A's)` = mean(polya_length), `% reads with polyA >= 10` = 100 * mean(polya_length >= 10),
             `mean transcript coverage` = mean(coverage), `% reads covering >= 90%` = 100 * mean(coverage >= 0.9),
             `5' missing (mean share)` = mean(five_missing), `3' missing (mean share)` = mean(three_missing),
             `read identity (%)` = 100 * mean(identity, na.rm = TRUE),
             `mean base quality` = mean(mean_quality, na.rm = TRUE),
             `% reads aligned once` = 100 * mean(n_alignments == 1)), by = .(data, version)]
over_k <- Reduce(function(a, b) merge(a, b, by = c("data", "version")),
                 list(ks[, lapply(.SD, mean), by = .(data, version)],
                      rl_sum[, .(data, version, `mean read length` = mean, `% reads < 300 bp` = short)], kq))
over_i <- is_[, .(data, version, `% reads mapped`, `error rate (%)`, `read length`, `mean insert size`, `mean base quality`)]
ov <- rbind(melt(over_k, id.vars = c("data", "version"), variable.name = "metric")[, platform := "Kinnex"],
            melt(over_i, id.vars = c("data", "version"), variable.name = "metric")[, platform := "Illumina"])
ov[, metric := paste0(platform, ": ", metric)]
ov[, metric := factor(metric, levels = rev(unique(metric)))]
ov[, label := formatC(value, digits = 3, format = "fg")]
ov[, stripe := as.integer(metric) %% 2 == 0]
p1 <- ggplot(ov, aes(x = lvl(data), y = metric)) +
    geom_tile(aes(fill = stripe), colour = NA, show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "grey94", `FALSE` = "white")) +
    geom_text(aes(label = label), size = 3) +
    facet_wrap(~version, nrow = 1, scales = "free_x") +
    labs(x = NULL, y = NULL, title = "Simulated reads vs real H69: overview",
         subtitle = paste("Kinnex: samtools stats of all reads (simulated: mean of 6 samples) + first 200k reads;",
                          "Illumina: first 500k pairs (control_rep1), minimap2 short-read mode")) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8))


# 2-10. Kinnex
rlb <- copy(rl)[, bin := 10^(floor(log10(pmax(length, 10)) * 40) / 40)]
rlb <- rlb[, .(n = sum(n)), by = .(data, version, bin)][, frac := n / sum(n), by = .(data, version)]
p2 <- ggplot(rlb, aes(x = bin, y = frac, colour = lvl(data))) +
    geom_step(linewidth = 0.6) + colour_scale +
    scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000)) +
    facet_wrap(~version, nrow = 1) +
    labs(x = "read length (bp, log scale)", y = "fraction of reads", title = "Kinnex read length",
         subtitle = paste(rl_sum[, sprintf("%s %s: mean %.0f bp, %.1f%% < 300 bp", data, version, mean, short)], collapse = "\n"))

p3 <- ggplot(qc, aes(x = transcript_length, y = read_length)) +
    rasterise(geom_bin2d(bins = 70), dpi = 200) +
    geom_abline(linewidth = 0.3, colour = "grey30", linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    scale_fill_viridis_c(trans = "log10", name = "reads") +
    facet_wrap(~ paste(data, "| aligned to", version), nrow = 2) +
    labs(x = "length of the transcript the read aligns to (bp)", y = "read length (bp)",
         title = "Kinnex read length vs transcript length (first 200k reads)",
         subtitle = "dashed: read as long as its transcript; below it: truncated reads; just above it: full-length reads with polyA tail")

p4 <- ggplot(qc[polya_length <= 60], aes(x = polya_length, colour = lvl(data))) +
    geom_freqpoly(aes(y = after_stat(density)), binwidth = 1, linewidth = 0.6) + colour_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "polyA tail length (A's at the read's 3' end)", y = "density", title = "Kinnex polyA tail",
         subtitle = paste(kq[, sprintf("%s %s: mean %.1f, %.0f%% >= 10", data, version,
                                       `polyA tail (mean A's)`, `% reads with polyA >= 10`)], collapse = "\n"))

tr <- melt(qc[, .(data, version, `transcript coverage` = coverage,
                  `5' end missing (share of transcript)` = five_missing,
                  `3' end missing (share of transcript)` = three_missing)],
           id.vars = c("data", "version"), variable.name = "what")
p5 <- ggplot(tr, aes(x = pmin(pmax(value, 0), 1), colour = lvl(data))) +
    geom_freqpoly(aes(y = after_stat(density)), binwidth = 0.02, boundary = 0, linewidth = 0.6) + colour_scale +
    facet_grid(version ~ what) +
    labs(x = NULL, y = "density", title = "Kinnex truncation (primary alignments, first 200k reads)",
         subtitle = "real reads' 3' 'loss' reflects alternative 3' ends vs the annotation (they end in a polyA tail), not truncation")

tb <- qc[, .(`mean transcript coverage` = mean(coverage), `% reads covering >= 90%` = 100 * mean(coverage >= 0.9),
             `5' missing (mean share)` = mean(five_missing), `% of reads` = .N), by = .(data, version, tx_bin)]
tb[, `% of reads` := 100 * `% of reads` / sum(`% of reads`), by = .(data, version)]
tbl <- melt(tb, id.vars = c("data", "version", "tx_bin"), variable.name = "metric")
tbl[, metric := factor(metric, levels = c("% of reads", "mean transcript coverage", "% reads covering >= 90%", "5' missing (mean share)"))]
p6 <- ggplot(tbl, aes(x = tx_bin, y = value, colour = lvl(data), group = lvl(data))) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) + colour_scale +
    facet_grid(metric ~ version, scales = "free_y", labeller = labeller(metric = label_wrap_gen(14))) +
    labs(x = "length of the transcript the read aligns to", y = NULL, title = "Kinnex truncation by transcript length",
         subtitle = "top: where the reads come from; below: how complete they are") +
    theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 8))

er <- melt(qc[, .(data, version, `read identity` = identity, `mismatch rate` = mismatch_rate,
                  `insertion rate` = insertion_rate, `deletion rate` = deletion_rate)],
           id.vars = c("data", "version"), variable.name = "what")[!is.na(value)]
er[, `:=`(lo = quantile(value, 0.01), hi = quantile(value, 0.99)), by = what]
er_sum <- qc[, lapply(.(identity, mismatch_rate, insertion_rate, deletion_rate), function(x) 100 * mean(x, na.rm = TRUE)),
             by = .(data, version)]
p7 <- ggplot(er[value >= lo & value <= hi], aes(x = value, colour = lvl(data))) +
    geom_freqpoly(aes(y = after_stat(density)), bins = 60, linewidth = 0.6) + colour_scale +
    facet_wrap(version ~ what, scales = "free", nrow = 2) +
    labs(x = "per read (central 98% shown)", y = "density", title = "Kinnex per-read identity and error types (from NM and CIGAR)",
         subtitle = paste(er_sum[, sprintf("%s %s: identity %.2f%%, mismatch %.3f%%, insertion %.3f%%, deletion %.3f%%",
                                           data, version, V1, V2, V3, V4)], collapse = "\n"))

p8 <- ggplot(qc[!is.na(mean_quality)], aes(x = mean_quality, colour = lvl(data))) +
    geom_freqpoly(aes(y = after_stat(density)), binwidth = 0.5, linewidth = 0.6) + colour_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "mean base quality of the read (Phred)", y = "density", title = "Kinnex base quality per read",
         subtitle = paste(kq[, sprintf("%s %s: mean %.1f", data, version, `mean base quality`)], collapse = "\n"))

na <- qc[, .(reads = .N), by = .(data, version, n = pmin(n_alignments, 11))][, frac := reads / sum(reads), by = .(data, version)]
p9 <- ggplot(na, aes(x = factor(n, levels = 1:11, labels = c(1:10, ">10")), y = frac, fill = lvl(data))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) + fill_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "alignments per read (primary + secondary + supplementary)", y = "fraction of reads",
         title = "Kinnex multi-mapping: how many transcripts a read aligns to",
         subtitle = paste(qc[, .(m = mean(n_alignments), u = 100 * mean(n_alignments == 1), pl = 100 * mean(strand == "+")),
                              by = .(data, version)][, sprintf("%s %s: mean %.2f alignments, %.0f%% aligned once, %.0f%% on + strand",
                                                              data, version, m, u, pl)], collapse = "\n"))

ksl <- melt(ks, id.vars = c("data", "version"), variable.name = "metric")
p10 <- ggplot(ksl, aes(x = version, y = value, colour = lvl(data))) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, position = dodge, show.legend = FALSE) +
    geom_point(position = dodge, size = 1.6, alpha = 0.8) + colour_scale +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    labs(x = "GENCODE version (reads aligned to the same version)", y = NULL,
         title = "Kinnex alignment summary (samtools stats, all reads)", subtitle = "dots: simulated samples / the real library; bar: mean")


# 11-16. Illumina
insn <- ins[insert <= 800, .(pairs = sum(pairs)), by = .(data, version, insert)][, frac := pairs / sum(pairs), by = .(data, version)]
p11 <- ggplot(insn, aes(x = insert, y = frac, colour = lvl(data))) +
    geom_line(linewidth = 0.6) + colour_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "insert (fragment) size (bp)", y = "fraction of pairs", title = "Illumina insert size",
         subtitle = paste(is_[, sprintf("%s %s: mean %.0f bp", data, version, `mean insert size`)], collapse = "\n"))

p12 <- ggplot(ffq, aes(x = cycle, y = value, colour = lvl(data))) +
    geom_line(linewidth = 0.6) + colour_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "position in read 1 (cycle)", y = "mean base quality (Phred)", title = "Illumina base quality along the read",
         subtitle = "simulated: constant quality (our polyester wrapper writes 'I' = Q40 everywhere)")

p13 <- if (nrow(mpc)) {
    ggplot(merge(mpc, reads_n, by = c("data", "version")), aes(x = cycle, y = 100 * value / reads, colour = lvl(data))) +
        geom_line(linewidth = 0.6) + colour_scale +
        facet_wrap(~version, nrow = 1) +
        labs(x = "position in the read (cycle)", y = "mismatches per 100 mapped reads",
             title = "Illumina mismatches along the read",
             subtitle = "simulated: uniform substitution errors at the configured rate")
} else ggplot() + labs(title = "Illumina mismatches along the read: not available (no MPC section)")

gcn <- gcf[, .(n = sum(n)), by = .(data, version, gc)][, frac := n / sum(n), by = .(data, version)]
p14 <- ggplot(gcn, aes(x = gc, y = frac, colour = lvl(data))) +
    geom_line(linewidth = 0.6) + colour_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "GC content of the read (%)", y = "fraction of reads", title = "Illumina GC content",
         subtitle = "simulated reads come from the expressed transcripts without a GC-bias model")

pp <- pos[transcript_length >= 1000][, rel := ((ref_start + ref_end) / 2) / transcript_length]
p15 <- ggplot(pp, aes(x = pmin(pmax(rel, 0), 1), colour = lvl(data))) +
    geom_freqpoly(aes(y = after_stat(density)), binwidth = 0.02, boundary = 0, linewidth = 0.6) + colour_scale +
    facet_wrap(~version, nrow = 1) +
    labs(x = "middle of read 1 along its transcript (0 = 5' end, 1 = 3' end)", y = "density",
         title = "Illumina read position along the transcript (transcripts >= 1 kb)",
         subtitle = "simulated: uniform fragment positions (polyester without positional bias)")

isl <- melt(is_[, .(data, version, `% reads mapped`, `error rate (%)`, `read length`)], id.vars = c("data", "version"),
            variable.name = "metric")
p16 <- ggplot(isl, aes(x = version, y = value, fill = lvl(data))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = signif(value, 3), group = lvl(data)), position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
    fill_scale + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    labs(x = "GENCODE version (reads aligned to the same version)", y = NULL,
         title = "Illumina alignment and error rate (minimap2 short-read mode, first 500k pairs)",
         subtitle = "simulated: control_rep1; real: H69 (101 bp). Real reads include intronic, intergenic and rRNA fragments")


# 17-18. both
ab_paths <- grep("/parametric/v[0-9]+/control_rep1_v[0-9]+\\.tsv$", paths, value = TRUE)
samp <- rbindlist(lapply(grep("true_counts\\.tsv$", paths, value = TRUE), function(f) {
    v <- version_of(f)
    ab <- fread(grep(paste0("/", v, "/"), ab_paths, value = TRUE))[tpm > 0]
    plat <- sub(".*/(kinnex|illumina)_([^/]+)/.*", "\\1", f)
    d <- sub(".*/(kinnex|illumina)_([^/]+)/.*", "\\2", f)
    w <- if (plat == "kinnex") ab$tpm else ab$est_counts
    e <- data.table(transcript_id = ab$target_id, expected = w / sum(w) * depth_num(d))
    e[fread(f), on = "transcript_id", true := i.true_reads][is.na(true), true := 0]
    e[, `:=`(version = v, platform = sprintf("%s (%s %s)", ifelse(plat == "kinnex", "Kinnex", "Illumina"), d,
                                              ifelse(plat == "kinnex", "reads", "pairs")))]
}))
r17 <- samp[, .(r = cor(log1p(expected), log1p(true))), by = .(version, platform)]
p17 <- ggplot(samp, aes(x = expected + 1, y = true + 1)) +
    geom_abline(linewidth = 0.3, colour = "grey60") +
    rasterise(geom_point(size = 0.3, alpha = 0.25, colour = VCOL[1]), dpi = 250) +
    geom_text(data = r17, aes(x = 1, y = Inf, label = sprintf("r (log) = %.3f", r)), hjust = 0, vjust = 1.5, size = 3.2) +
    scale_x_log10() + scale_y_log10() +
    facet_grid(version ~ platform) +
    labs(x = "expected reads from the true abundance + 1", y = "reads actually simulated + 1",
         title = "Reads per transcript: simulated vs expected (control_rep1)",
         subtitle = "scatter around the diagonal = sampling noise of sequencing")

oar <- rbindlist(lapply(grep("run_oarfish/.*\\.quant$", paths, value = TRUE), function(f) {
    dl <- sub("^sim[0-9]+_kinnex([^_]+)_.*", "\\1", basename(f))
    v <- sub(".*\\.(v[0-9]+)\\.quant$", "\\1", f)
    data.table(data = paste("simulated", v), version = v, platform = "Kinnex (oarfish)",
               pct = 100 * sum(fread(f, select = "num_reads")$num_reads) / depth_num(dl))
}))
sal <- rbindlist(lapply(grep("meta_info\\.json$", paths, value = TRUE), function(f) {
    j <- paste(readLines(f), collapse = "")
    data.table(data = data_label(f), version = version_of(f), platform = "Illumina (salmon)",
               pct = as.numeric(sub('.*"percent_mapped": ([0-9.]+).*', "\\1", j)))
}))
use <- rbind(oar, sal, ks[data == "real H69", .(data, version, platform = "Kinnex (oarfish)", pct = `% reads mapped`)])
use[, platform := factor(platform, levels = c("Kinnex (oarfish)", "Illumina (salmon)"))]
p18 <- ggplot(use, aes(x = version, y = pct, colour = lvl(data))) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, position = dodge, show.legend = FALSE) +
    geom_point(position = dodge, size = 1.8, alpha = 0.8) + colour_scale +
    facet_wrap(~platform, nrow = 1, scales = "free_y") +
    labs(x = "GENCODE version (reads quantified against the same version)", y = "% of reads used",
         title = "Share of reads the quantifiers use",
         subtitle = "simulated Kinnex: reads assigned by oarfish / reads simulated; real Kinnex: % mapped (minimap2); Illumina: salmon mapping rate")


# log + pages
cat("\n== overview ==\n"); print(dcast(ov, metric ~ data + version, value.var = "value"), digits = 3)
cat("\n== truncation by transcript length ==\n"); print(tb, digits = 3)
cat("\n== sampling (r, log) ==\n"); print(r17, digits = 4)
pdf(args[[2]], width = 12, height = 6.5)
for (p in list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18)) print(p)
dev.off()
