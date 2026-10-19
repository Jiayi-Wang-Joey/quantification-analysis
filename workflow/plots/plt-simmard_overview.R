suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})


# MARD of simulated samples against their true TPM (miniQuant), all transcripts and >= min_length

paths <- args[[1]]
summ <- rbindlist(lapply(grep("\\.summary\\.tsv$", paths, value = TRUE), fread))
tx <- rbindlist(lapply(grep("\\.transcripts\\.tsv$", paths, value = TRUE), fread), fill = TRUE)
min_len <- summ$min_length[1]

# one colour per platform x depth
depth_num <- function(x) as.numeric(sub("M$", "", x)) * fifelse(grepl("M$", x), 1e6, 1)
SERIES <- unique(tx[, .(platform, depth)])[order(platform != "kinnex", depth_num(depth))]
SERIES[, label := sprintf("%s %s (%s)", fifelse(platform == "kinnex", "Kinnex", "Illumina"),
                          fifelse(platform == "kinnex", paste(depth, "reads"), paste(depth, "pairs")),
                          fifelse(platform == "kinnex", "oarfish", "salmon"))]
BLUES <- c("#86b6ef", "#2a78d6", "#184f95", "#0d366b")
ORANGES <- c("#eb6834", "#b8461d", "#f2a37f")
SERIES[platform == "kinnex", colour := BLUES[seq_len(.N)]]
SERIES[platform == "illumina", colour := ORANGES[seq_len(.N)]]
COLOURS <- setNames(SERIES$colour, SERIES$label)
tx[SERIES, on = c("platform", "depth"), series := i.label]
tx[, platform := factor(series, levels = SERIES$label)]
REFS <- unique(tx$reference)
REFS <- REFS[order(-as.integer(sub("v", "", REFS)))]
tx[, reference := factor(reference, levels = REFS)]
REF_LABEL <- setNames(sprintf("GENCODE %s (truth and reference)", REFS), REFS)
N_SAMPLES <- uniqueN(tx$sample)
of_samples <- if (N_SAMPLES == 1) "1 sample" else sprintf("mean of %d samples", N_SAMPLES)
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                           panel.grid.major = element_line(linewidth = 0.2),
                                           strip.background = element_blank(),
                                           legend.position = "top"))
dodge <- position_dodge(width = 0.5)
K_LABELS <- c("=1", "1-2", "2-3", "3-5", "5-8", "8-13", "13-25", "25-90", ">=90")
k_bin <- function(k) factor(fcase(k == 1, "=1", k <= 2, "1-2", k <= 3, "2-3", k <= 5, "3-5", k <= 8, "5-8",
                                  k <= 13, "8-13", k <= 25, "13-25", k < 90, "25-90", default = ">=90"),
                            levels = K_LABELS)
E_LABELS <- c("< 1", "1-10", "10-100", "100-1000", ">= 1000")
e_bin <- function(x) factor(fcase(x < 1, "< 1", x < 10, "1-10", x < 100, "10-100", x < 1000, "100-1000",
                                  default = ">= 1000"), levels = E_LABELS)

SETS <- list(
    list(name = "all transcripts (no length cut-off)", rows = rep(TRUE, nrow(tx)),
         ard = "ard", true_tpm = "tpm_true", est_tpm = "tpm_est"),
    list(name = sprintf("transcripts >= %d bp", min_len), rows = tx$in_eval,
         ard = "ard_eval", true_tpm = "tpm_true_eval", est_tpm = "tpm_est_eval")
)


pages_for <- function(set) {
    t <- tx[set$rows]
    t[, `:=`(a = get(set$ard), ts = get(set$true_tpm), es = get(set$est_tpm))]

    g <- t[ts > 0, .(gene_mard = mean(a, na.rm = TRUE), n_iso = .N, gene_tpm = sum(ts),
                       kvalue = kvalue[1]), by = .(platform, reference, sample, gene)]
    g[, kbin := k_bin(kvalue)]
    # false signal
    fs <- t[, .(est_unexpressed = 100 * sum(es[status == "not_expressed"], na.rm = TRUE) / sum(es, na.rm = TRUE),
                truth_missing = 100 * sum(ts[status == "missing_from_reference"], na.rm = TRUE) / sum(ts, na.rm = TRUE)),
            by = .(platform, reference, sample)]

    head_tab <- merge(g[, .(median_gene_mard = median(gene_mard, na.rm = TRUE), genes = .N),
                        by = .(platform, reference, sample)],
                      fs, by = c("platform", "reference", "sample"))
    cat("\n==", set$name, "== (", of_samples, ")\n")
    print(head_tab[, lapply(.SD, mean), by = .(platform, reference), .SDcols = -"sample"], digits = 3)

    # 1. headline
    h <- melt(head_tab[, .(platform, reference, sample,
                           `accuracy: median gene MARD\n(expressed isoforms)` = median_gene_mard,
                           `false signal: % of estimate\non isoforms with no true reads` = est_unexpressed,
                           `% of truth on isoforms\nmissing from the reference` = truth_missing)],
              id.vars = c("platform", "reference", "sample"), variable.name = "metric")
    p1 <- ggplot(h, aes(x = reference, y = value, colour = platform)) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.3, position = dodge, show.legend = FALSE) +
        geom_point(position = dodge, size = 1.6, alpha = 0.8) +
        stat_summary(fun = mean, geom = "text", aes(label = signif(after_stat(y), 2)), vjust = -1.2,
                     size = 3, position = dodge, show.legend = FALSE) +
        scale_colour_manual(values = COLOURS, name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        facet_wrap(~metric, nrow = 1, scales = "free_y") +
        labs(x = "GENCODE version (truth and reference)", y = NULL,
             title = sprintf("Accuracy and false signal of the simulated samples: %s", set$name),
             subtitle = paste(sprintf("truth: each sample's true TPM (same for both platforms); dots: %d sample(s), bar and label: mean.", N_SAMPLES),
                              "\nAccuracy as in miniQuant: per gene, mean ARD over its expressed isoforms (TPM >= 0.01); median over genes"))

    # 2. K-value group sizes
    kc <- g[!is.na(kbin) & platform == levels(platform)[1],
            .(genes = .N, `expressed isoforms` = sum(n_iso)), by = .(reference, sample, kbin)]
    kc <- kc[, .(genes = mean(genes), `expressed isoforms` = mean(`expressed isoforms`)), by = .(reference, kbin)]
    kc <- melt(kc, id.vars = c("reference", "kbin"), variable.name = "what", value.name = "n")
    p2 <- ggplot(kc, aes(x = kbin, y = n, fill = reference)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.75) +
        geom_text(aes(label = round(n), group = reference), position = position_dodge(width = 0.8),
                  vjust = -0.3, size = 2.6, colour = "grey30") +
        scale_fill_manual(values = setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7", "#8f8e88")[seq_along(REFS)], REFS), name = "GENCODE") +
        facet_wrap(~what, nrow = 1, scales = "free_y") +
        labs(x = "K-value group of the gene", y = sprintf("count (%s)", of_samples),
             title = sprintf("Expressed genes and isoforms per K-value group: %s", set$name),
             subtitle = "the same expression design in both versions; v50 annotates more similar isoforms, so more genes fall into high K-value groups")

    # 3. gene MARD per K-value group
    gk <- g[!is.na(kbin), .(gene_mard = mean(gene_mard, na.rm = TRUE)), by = .(platform, reference, gene, kbin)]
    p3 <- ggplot(gk, aes(x = kbin, y = gene_mard, colour = platform)) +
        geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.4, position = position_dodge(width = 0.7)) +
        scale_colour_manual(values = COLOURS, name = NULL) +
        facet_wrap(~reference, nrow = 1, labeller = as_labeller(REF_LABEL)) +
        labs(x = "K-value group of the gene (in the reference quantified against)", y = "gene MARD (expressed isoforms)",
             title = sprintf("Gene MARD per K-value group: %s", set$name),
             subtitle = sprintf("each gene once (%s); box: quartiles, line: median; outliers not drawn", of_samples))

    # 3b. genes with true TPM > 1
    gk1 <- g[!is.na(kbin), .(gene_mard = mean(gene_mard, na.rm = TRUE), gene_tpm = mean(gene_tpm)),
             by = .(platform, reference, gene, kbin)][gene_tpm > 1]
    gk1_n <- gk1[platform == levels(platform)[1], .N, by = .(reference, kbin)]
    p3b <- ggplot(gk1, aes(x = kbin, y = gene_mard, colour = platform)) +
        geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.4, position = position_dodge(width = 0.7)) +
        geom_text(data = gk1_n, aes(x = kbin, y = -0.04, label = N), inherit.aes = FALSE, size = 2.8, colour = "grey40") +
        scale_colour_manual(values = COLOURS, name = NULL) +
        facet_wrap(~reference, nrow = 1, labeller = as_labeller(REF_LABEL)) +
        labs(x = "K-value group of the gene (in the reference quantified against)", y = "gene MARD (expressed isoforms)",
             title = sprintf("Gene MARD per K-value group, genes with true TPM > 1: %s", set$name),
             subtitle = paste(sprintf("as miniQuant, to reduce sampling error; gene TPM = sum over its isoforms (%s);", of_samples),
                              "grey: number of genes; outliers not drawn"))

    # 4. by expression
    ge <- g[, .(gene_mard = mean(gene_mard, na.rm = TRUE), gene_tpm = mean(gene_tpm)), by = .(platform, reference, gene)]
    ge[, ebin := e_bin(gene_tpm)]
    ge_med <- ge[, .(med = median(gene_mard, na.rm = TRUE)), by = .(platform, reference, ebin)]
    p4 <- ggplot(ge_med, aes(x = ebin, y = med, colour = platform, group = platform)) +
        geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
        scale_colour_manual(values = COLOURS, name = NULL) +
        scale_y_continuous(limits = c(0, NA)) +
        facet_wrap(~reference, nrow = 1, labeller = as_labeller(REF_LABEL)) +
        labs(x = "true expression of the gene in the sample (TPM, sum over its isoforms)", y = "median gene MARD",
             title = sprintf("Median gene MARD by expression level: %s", set$name),
             subtitle = sprintf("each gene once (%s)", of_samples))

    list(p1, p2, p3, p3b, p4)
}

pages <- lapply(SETS, pages_for)
pdf(args[[2]], width = 12, height = 6)
for (i in seq_along(pages[[1]])) for (s in seq_along(SETS)) print(pages[[s]][[i]])
dev.off()
