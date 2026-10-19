suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrastr)
})


# simulated truth (one version): expression, control vs treatment, replicates, isoform shares

paths <- args[[1]]
labels <- fread(grep("de_labels_v", paths, value = TRUE))
abundance_paths <- grep("(control|treatment)_(abundance|rep[0-9]+)_v[0-9]+\\.tsv$", paths, value = TRUE)
version <- sub(".*_v([0-9]+)\\.tsv$", "\\1", abundance_paths[1])

read_tpm <- function(f) {
    d <- fread(f, select = c("target_id", "tpm"))
    name <- sub("_v[0-9]+\\.tsv$", "", basename(f))
    d[, `:=`(condition = sub("_.*", "", name),
             file = fifelse(grepl("abundance", name), "truth", sub(".*_", "", name)))]
    d
}
tpm <- rbindlist(lapply(abundance_paths, read_tpm))
tpm[, file := factor(file, levels = c("truth", sort(setdiff(unique(file), "truth"))))]
tpm[, condition := factor(condition, levels = c("control", "treatment"))]

BLUE <- "#2a78d6"
GREY <- "#8f8e88"
ORANGE <- "#eb6834"
CLASS_LEVELS <- c("null", "uniform", "dtu_swap", "partial_dte")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                           panel.grid.major = element_line(linewidth = 0.2),
                                           strip.background = element_blank()))


# 1. expression distribution of every file -------------------------------------
expr <- tpm[tpm > 0]
xlim <- quantile(log10(expr$tpm), c(0.0005, 1))
p1 <- ggplot(expr, aes(x = log10(tpm), colour = condition)) +
    geom_freqpoly(bins = 60, linewidth = 0.7) +
    scale_colour_manual(values = c(control = BLUE, treatment = ORANGE), name = NULL) +
    facet_wrap(~file, nrow = 1) +
    coord_cartesian(xlim = xlim) +
    labs(x = "log10 TPM (expressed transcripts)", y = "transcripts",
         title = "Transcript expression in the simulated truth and replicate files",
         subtitle = sprintf("GENCODE v%s; %d expressed transcripts in %d genes; 5%% of genes differ between conditions",
                            version, nrow(labels), uniqueN(labels$gene_id))) +
    theme(legend.position = "top")


# 2. control vs treatment truth, by change type --------------------------------
truth <- dcast(tpm[file == "truth"], target_id ~ condition, value.var = "tpm")
truth <- merge(truth, labels, by.x = "target_id", by.y = "transcript_id")
truth[, class := factor(class, levels = CLASS_LEVELS)]
n_genes <- truth[, .(n = uniqueN(gene_id)), by = class]
truth[n_genes, on = "class", facet := sprintf("%s (%d genes)", class, i.n)]
truth[, facet := factor(facet, levels = n_genes[order(match(class, CLASS_LEVELS)), sprintf("%s (%d genes)", class, n)])]
truth[, changed := factor(fifelse(is_de, "changed", "unchanged"), levels = c("unchanged", "changed"))]
setorder(truth, changed)
p2 <- ggplot(truth, aes(x = log10(control), y = log10(treatment), colour = changed)) +
    geom_abline(linewidth = 0.3, colour = "grey60") +
    rasterise(geom_point(size = 0.6, alpha = 0.6), dpi = 300) +
    scale_colour_manual(values = c(unchanged = GREY, changed = BLUE), name = NULL) +
    facet_wrap(~facet, nrow = 1) +
    coord_equal() +
    labs(x = "control truth, log10 TPM", y = "treatment truth, log10 TPM",
         title = "Control vs treatment truth by change type",
         subtitle = "uniform: whole gene shifts; dtu_swap: rank 1 and 2 trade places; partial_dte: 1-3 isoforms shift") +
    theme(legend.position = "top")


# 2b. true log2FC of the transcripts that change, per change type
de <- truth[is_de == TRUE & class != "null"]
de[, facet := factor(as.character(facet), levels = levels(truth$facet)[-1])]
de_n <- de[, .N, by = facet]
p2b <- ggplot(de, aes(x = true_log2fc)) +
    geom_histogram(binwidth = 0.1, fill = BLUE, colour = NA) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_text(data = de_n, aes(x = -Inf, y = Inf, label = sprintf("%d transcripts", N)),
              hjust = -0.1, vjust = 1.5, size = 3.4, colour = "grey20") +
    facet_wrap(~facet, nrow = 1) +
    labs(x = "true log2 fold change (treatment vs control)", y = "transcripts",
         title = "True log2FC of the changed transcripts, per change type",
         subtitle = paste("uniform / partial_dte: |log2FC| drawn from 0.58-2.58 (1.5- to 6-fold);",
                          "dtu_swap: +/- log2(rank-1 share / rank-2 share)"))


# 3. replicates vs their truth --------------------------------------------------
reps <- merge(tpm[file != "truth"], tpm[file == "truth", .(target_id, condition, truth = tpm)],
              by = c("target_id", "condition"))[truth > 0]
reps <- merge(reps, labels[, .(target_id = transcript_id, rank)], by = "target_id")
reps[, isoform := factor(fifelse(rank == 1, "rank 1 (most expressed)", "rank 2+"),
                         levels = c("rank 1 (most expressed)", "rank 2+"))]
p3 <- ggplot(reps, aes(x = log10(truth), y = log10(tpm))) +
    geom_abline(linewidth = 0.3, colour = "grey60") +
    rasterise(geom_point(size = 0.4, alpha = 0.3, colour = BLUE), dpi = 300) +
    facet_grid(isoform ~ condition + file) +
    labs(x = "truth, log10 TPM", y = "replicate, log10 TPM",
         title = "Replicates scatter around their condition truth",
         subtitle = "gene noise (bcv) moves whole genes; isoform-ratio noise (share_precision) mostly moves minor isoforms")


# 4. isoform shares by rank and isoforms per gene -------------------------------
ctrl <- truth[, .(gene_id, rank, control)]
ctrl[, share := control / sum(control), by = gene_id]
ctrl[, rank_bin := factor(fifelse(rank >= 5, "5+", as.character(rank)), levels = c("1", "2", "3", "4", "5+"))]
p4a <- ggplot(ctrl, aes(x = rank_bin, y = share)) +
    geom_boxplot(fill = NA, colour = BLUE, outlier.size = 0.3, linewidth = 0.4) +
    labs(x = "isoform rank within gene", y = "share of the gene (control truth)",
         title = "Isoform shares by rank (Dirichlet, sorted)")
k <- ctrl[, .(k = .N), by = gene_id]
p4b <- ggplot(k, aes(x = factor(pmin(k, 8), levels = 1:8, labels = c(1:7, "8+")))) +
    geom_bar(fill = BLUE, width = 0.7) +
    labs(x = "expressed isoforms per gene", y = "genes",
         title = sprintf("Expressed isoforms per gene (mean %.2f)", mean(k$k)))


# 5. the control truth itself ----------------------------------------------------
ctrl_tx <- merge(fread(grep("control_abundance_v", abundance_paths, value = TRUE),
                       select = c("target_id", "length", "tpm"))[tpm > 0],
                 labels[, .(target_id = transcript_id, gene_id)], by = "target_id")
gene_tpm <- ctrl_tx[, .(tpm = sum(tpm)), by = gene_id][order(-tpm)]
gene_tpm[, `:=`(gene_rank = seq_len(.N), cum_share = cumsum(tpm) / sum(tpm))]
top <- gene_tpm[gene_rank %in% c(10, 100, 1000)]
p5a <- ggplot(gene_tpm, aes(x = gene_rank, y = cum_share)) +
    geom_line(colour = BLUE, linewidth = 0.7) +
    geom_point(data = top, colour = BLUE, size = 2.5) +
    geom_text(data = top, aes(label = sprintf("top %d genes: %.0f%%", gene_rank, 100 * cum_share)),
              hjust = -0.1, vjust = 1.4, size = 3.2, colour = "grey20") +
    scale_x_log10() +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
    labs(x = "genes ranked by expression (log scale)", y = "cumulative share of all TPM",
         title = "How concentrated is the control truth? (gene level)",
         subtitle = sprintf("%d expressed genes; median gene TPM %.1f", nrow(gene_tpm), median(gene_tpm$tpm)))

ctrl_tx[, size := factor(fifelse(length < 250, "< 250 bp (shorter than an Illumina fragment)", ">= 250 bp"),
                         levels = c(">= 250 bp", "< 250 bp (shorter than an Illumina fragment)"))]
p5b <- ggplot(ctrl_tx, aes(x = length, y = tpm, colour = size)) +
    rasterise(geom_point(size = 0.5, alpha = 0.4), dpi = 300) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(BLUE, ORANGE), name = NULL) +
    labs(x = "transcript length (bp, log scale)", y = "control truth TPM (log scale)",
         title = "Expression vs transcript length (control truth)",
         subtitle = sprintf("%d of %d expressed transcripts are < 250 bp", sum(ctrl_tx$length < 250), nrow(ctrl_tx))) +
    theme(legend.position = "top")

thr <- data.table(threshold = c(0.01, 0.1, 1, 10))
thr[, n := sapply(threshold, function(t) sum(ctrl_tx$tpm >= t))]
thr[, label := factor(paste("TPM >=", threshold), levels = paste("TPM >=", threshold))]
p5c <- ggplot(thr, aes(x = label, y = n)) +
    geom_col(fill = BLUE, width = 0.6) +
    geom_text(aes(label = sprintf("%d (%.0f%%)", n, 100 * n / nrow(ctrl_tx))), vjust = -0.4, size = 3.4,
              colour = "grey20") +
    labs(x = NULL, y = "expressed transcripts",
         title = "Expressed transcripts passing common TPM thresholds (control truth)",
         subtitle = "miniQuant evaluates TPM >= 0.01 and uses TPM > 1 for its K-value analyses")


pdf(args[[2]], width = 12, height = 6)
print(p1)
print(p2)
print(p2b)
print(p3)
print(p4a)
print(p4b)
print(p5a)
print(p5b)
print(p5c)
dev.off()
