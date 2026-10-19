suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

# rough fit of the truth's distributions (expression, isoforms per gene, shares, gene types) to real quantifications

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
p <- snakemake@params
set.seed(1)
say <- function(...) message(sprintf(...))

SMALL <- c("snRNA", "snoRNA", "miRNA", "misc_RNA", "scaRNA", "scRNA", "sRNA", "vault_RNA", "Y_RNA", "rRNA",
           "rRNA_pseudogene", "Mt_tRNA", "Mt_rRNA", "ribozyme")
type_group <- function(t) fcase(t == "protein_coding", "protein_coding", t == "lncRNA", "lncRNA",
                                t %in% SMALL, "small RNA", grepl("pseudogene", t), "pseudogene", default = "other")
TYPES <- c("protein_coding", "lncRNA", "pseudogene", "small RNA", "other")

genes <- fread(cmd = paste("zcat", snakemake@input[["gtf"]], "| awk -F'\\t' '$3 == \"transcript\"' |",
                           "grep -oE 'gene_id \"[^\"]+\"|transcript_id \"[^\"]+\"|gene_type \"[^\"]+\"' |",
                           "cut -d'\"' -f2 | paste - - -"), header = FALSE, col.names = c("gene_id", "transcript_id", "gene_type"))
genes[, group := type_group(gene_type)]
n_annot <- genes[, .(n_annot = .N), by = gene_id]
say("%d transcripts in %d genes", nrow(genes), nrow(n_annot))

fit_lnorm_trunc <- function(x, t) {
    nll <- function(q) -sum(dnorm(log(x), q[1], exp(q[2]), log = TRUE) - pnorm(log(t), q[1], exp(q[2]), lower.tail = FALSE, log.p = TRUE))
    o <- optim(c(mean(log(x)), log(sd(log(x)))), nll)
    list(meanlog = o$par[1], sdlog = exp(o$par[2]), p_detect = pnorm(log(t), o$par[1], exp(o$par[2]), lower.tail = FALSE))
}
# k - 1 censored at n_annot - 1
fit_count <- function(k, n, d, p_) {
    ll <- function(par) sum(ifelse(k < n, d(k - 1, par, log = TRUE), log(pmax(1 - p_(n - 2, par), 1e-300))))
    if (identical(d, dpois)) optimize(ll, c(0.01, 20), maximum = TRUE) else optimize(ll, c(0.01, 0.99), maximum = TRUE)
}
fit_dirichlet <- function(shares) {
    ll <- function(a) sum(sapply(shares, function(s) lgamma(length(s) * a) - length(s) * lgamma(a) + (a - 1) * sum(log(s))))
    optimize(ll, c(0.01, 20), maximum = TRUE)$maximum
}
draw_shares <- function(k, alpha) {
    s <- rgamma(k, alpha); s <- pmax(s / sum(s), p$min_share); sort(s / sum(s), decreasing = TRUE)
}

fits <- list(); real <- list()
for (f in unlist(snakemake@input[["quant"]])) {
    if (grepl("quant\\.sf$", f)) {
        smp <- paste(basename(dirname(f)), "short (salmon)")
        q <- fread(f, select = c("Name", "NumReads", "TPM"), col.names = c("transcript_id", "reads", "tpm"))
    } else {
        smp <- paste(sub("\\.v[0-9]+\\.quant$", "", basename(f)), "long (oarfish)")
        q <- fread(f, select = c("tname", "num_reads"), col.names = c("transcript_id", "reads"))
        q[, tpm := reads / sum(reads) * 1e6]
    }
    q[genes, on = "transcript_id", `:=`(gene_id = i.gene_id, group = i.group)]
    g <- q[, .(reads = sum(reads), tpm = sum(tpm), group = group[1]), by = gene_id][n_annot, on = "gene_id", nomatch = 0]

    t <- p$min_gene_reads / sum(q$reads) * 1e6
    det <- g[tpm >= t]
    ln <- fit_lnorm_trunc(det$tpm, t)

    iso <- q[gene_id %in% g[reads >= p$iso_gene_reads]$gene_id][, share := tpm / sum(tpm), by = gene_id][share >= p$min_share]
    kd <- iso[, .(k = .N, shares = list(sort(share / sum(share), decreasing = TRUE))), by = gene_id][n_annot, on = "gene_id", nomatch = 0]
    kd_info <- kd[n_annot >= 2]
    po <- fit_count(kd_info$k, kd_info$n_annot, dpois, ppois)
    ge <- fit_count(kd_info$k, kd_info$n_annot, dgeom, pgeom)
    alpha <- fit_dirichlet(kd[k >= 2]$shares)

    types <- det[, .N, by = group][, frac := N / sum(N)]
    type_frac <- setNames(as.data.table(as.list(types[match(TYPES, group), fifelse(is.na(frac), 0, frac)])), paste0("frac_", TYPES))
    # stricter 'expressed' thresholds
    strict <- sapply(c(0.05, 0.1), function(m) {
        kk <- kd_info[, .(k = sapply(shares, function(s) sum(s >= m)), n_annot)][k >= 1]
        fit_count(kk$k, kk$n_annot, dpois, ppois)$maximum
    })
    fits[[smp]] <- data.table(sample = smp, total_reads = sum(q$reads), genes_detected = nrow(det),
                              sdlog = ln$sdlog, implied_genes = round(nrow(det) / ln$p_detect),
                              genes_for_isoforms = nrow(kd), mean_k = mean(kd$k), pct_single = 100 * mean(kd$k == 1),
                              lambda_poisson = po$maximum, loglik_poisson = po$objective,
                              mean_geometric = (1 - ge$maximum) / ge$maximum, loglik_geometric = ge$objective,
                              lambda_poisson_5pct = strict[1], lambda_poisson_10pct = strict[2],
                              alpha = alpha, mean_top_share = mean(sapply(kd[k >= 2]$shares, `[`, 1)), type_frac)
    real[[smp]] <- list(det = det[, .(sample = smp, tpm, t)], kd = kd[, .(sample = smp, k, n_annot, top = sapply(shares, `[`, 1))],
                        ln = ln, lambda = po$maximum, alpha = alpha)
    say("%s: %d genes detected, sdlog %.2f, lambda %.2f, alpha %.2f", smp, nrow(det), ln$sdlog, po$maximum, alpha)
}
fit <- rbindlist(fits)
fwrite(fit, snakemake@output[["table"]], sep = "\t")
sink(type = "message")
print(fit, digits = 3)


# figure
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
                                           legend.position = "top", plot.subtitle = element_text(size = 8)))
COL <- c(real = "#8f8e88", fitted = "#2a78d6", current = "#eb6834")
det <- rbindlist(lapply(real, `[[`, "det"))
curves <- rbindlist(lapply(names(real), function(s) {
    r <- real[[s]]; x <- seq(log10(r$det$t[1]), log10(max(r$det$tpm)), length.out = 200)
    dens <- function(m, sd) dnorm(x * log(10), m, sd) * log(10) / pnorm(log(r$det$t[1]), m, sd, lower.tail = FALSE)
    rbind(data.table(sample = s, x = x, y = dens(r$ln$meanlog, r$ln$sdlog), what = "fitted"),
          data.table(sample = s, x = x, y = dens(r$ln$meanlog, p$current$sdlog), what = "current"))
}))
fit_lab <- fit[, .(sample, label = sprintf("sdlog %.2f\n%d detected\n~%d genes in total", sdlog, genes_detected, implied_genes))]
p1 <- ggplot(det, aes(x = log10(tpm))) +
    geom_histogram(aes(y = after_stat(density)), bins = 60, fill = COL["real"], colour = NA) +
    geom_line(data = curves, aes(x = x, y = y, colour = what), linewidth = 0.6) +
    geom_text(data = fit_lab, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.2, size = 2.6) +
    scale_colour_manual(values = COL, labels = c(fitted = "fitted log-normal", current = sprintf("current sdlog = %g", p$current$sdlog)), name = NULL) +
    facet_wrap(~sample) +
    labs(x = "gene TPM (log10), genes with >= min reads", y = "density", title = "Gene expression: real vs log-normal",
         subtitle = sprintf("log-normal fitted to genes with >= %d reads, accounting for the undetected genes below that (truncated fit; same mean for 'current')", p$min_gene_reads))

kd <- rbindlist(lapply(real, `[[`, "kd"))
kobs <- kd[, .N, by = .(sample, k = pmin(k, 8))][, frac := N / sum(N), by = sample]
ksim <- rbindlist(lapply(names(real), function(s) {
    n <- real[[s]]$kd$n_annot
    rbind(data.table(sample = s, k = pmin(pmin(1 + rpois(length(n), real[[s]]$lambda), n), 8), what = "fitted"),
          data.table(sample = s, k = pmin(pmin(1 + rpois(length(n), p$current$lambda), n), 8), what = "current"))
}))[, .N, by = .(sample, what, k)][, frac := N / sum(N), by = .(sample, what)]
p2 <- ggplot(kobs, aes(x = factor(k), y = frac)) +
    geom_col(fill = COL["real"], width = 0.8) +
    geom_point(data = ksim, aes(colour = what), position = position_dodge(width = 0.5), size = 1.8) +
    scale_colour_manual(values = COL, labels = c(fitted = "1 + Poisson(fitted lambda)", current = sprintf("current lambda = %g", p$current$lambda)), name = NULL) +
    scale_x_discrete(labels = c(1:7, ">=8")) +
    facet_wrap(~sample) +
    labs(x = "expressed isoforms per gene", y = "fraction of genes", title = "Expressed isoforms per gene: real vs 1 + Poisson",
         subtitle = sprintf("genes with >= %d reads; isoform expressed if >= %g%% of its gene; simulated with each gene's annotated isoform count as the cap",
                            p$iso_gene_reads, 100 * p$min_share))

top <- kd[k >= 2, .(real = mean(top)), by = .(sample, k = pmin(k, 6))]
tsim <- rbindlist(lapply(names(real), function(s) {
    r <- real[[s]]$kd[k >= 2]
    rbind(data.table(sample = s, k = pmin(r$k, 6), top = sapply(r$k, function(k) draw_shares(k, real[[s]]$alpha)[1]), what = "fitted"),
          data.table(sample = s, k = pmin(r$k, 6), top = sapply(r$k, function(k) draw_shares(k, p$current$alpha)[1]), what = "current"))
}))[, .(top = mean(top)), by = .(sample, what, k)]
p3 <- ggplot(top, aes(x = factor(k), y = real)) +
    geom_col(fill = COL["real"], width = 0.8) +
    geom_point(data = tsim, aes(y = top, colour = what), position = position_dodge(width = 0.5), size = 1.8) +
    scale_colour_manual(values = COL, labels = c(fitted = "Dirichlet(fitted alpha)", current = sprintf("current alpha = %g", p$current$alpha)), name = NULL) +
    scale_x_discrete(labels = c(2:5, ">=6")) +
    facet_wrap(~sample) +
    labs(x = "expressed isoforms in the gene", y = "mean share of the top isoform", title = "Isoform shares: how dominant is the top isoform?",
         subtitle = sprintf("bars: real; dots: shares drawn from Dirichlet(alpha) with the %g%% floor, for the same genes", 100 * p$min_share))

ty <- melt(fit[, c("sample", paste0("frac_", TYPES)), with = FALSE], id.vars = "sample")
ty[, type := factor(sub("frac_", "", variable), levels = TYPES)]
ann_ty <- genes[, .(group = group[1]), by = gene_id][, .N, by = group][, .(sample = "GENCODE (all genes)", type = factor(group, levels = TYPES), value = N / sum(N))]
p4 <- ggplot(rbind(ty[, .(sample, type, value)], ann_ty), aes(x = sample, y = value, fill = type)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = ifelse(value > 0.03, sprintf("%.0f%%", 100 * value), "")), position = position_stack(vjust = 0.5), size = 2.6, colour = "white") +
    scale_fill_manual(values = c(protein_coding = "#2a78d6", lncRNA = "#eb6834", pseudogene = "#1baf7a", `small RNA` = "#4a3aa7", other = "#8f8e88"), name = NULL) +
    labs(x = NULL, y = "share of genes", title = "Gene types of expressed genes",
         subtitle = sprintf("genes with >= %d reads; right: all annotated genes (what a random draw gives)", p$min_gene_reads)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(snakemake@output[["pdf"]], width = 12, height = 7)
for (pl in list(p1, p2, p3, p4)) print(pl)
dev.off()
