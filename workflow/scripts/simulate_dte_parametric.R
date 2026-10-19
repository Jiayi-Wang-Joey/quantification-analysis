
# parametric DTE truth: genes per gene type, log-normal expression, 1 + Poisson isoforms,
# Dirichlet shares; treatment classes uniform / dtu_swap / partial_dte; replicates with noise


suppressPackageStartupMessages(library(data.table))

p <- snakemake@params[["p"]]
p$abs_log2fc <- unlist(p$abs_log2fc)
p$versions <- unlist(p$versions)
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
set.seed(p$seed)
dir.create(p$outdir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf(...))
unver <- function(id) sub("\\.[0-9]+(_PAR_Y)?$", "", id)
draw <- function(x, n) x[sample.int(length(x), n)]


# transcript length = sum of exons
load_annotation <- function(gtf) {
    exons <- as.data.table(rtracklayer::import(gtf, feature.type = "exon"))
    a <- exons[, .(length = sum(width)), by = .(transcript_id, gene_id, gene_type)]
    setorder(a, transcript_id)
    a[, `:=`(gene_key = unver(gene_id), eff_length = pmax(length - p$frag_mean + 1, 1))]
    a
}
ann <- setNames(lapply(unlist(snakemake@input[["gtf"]]), load_annotation), p$versions)
for (v in names(ann)) say("v%s: %d transcripts in %d genes", v, nrow(ann[[v]]), uniqueN(ann[[v]]$gene_id))
usable <- Reduce(intersect, lapply(ann, function(a) unique(a$gene_key)))
say("%d genes exist in every version", length(usable))


# 1. control truth (genes x ranks)
# gene type from the newest version
gtype <- unique(ann[[length(ann)]][gene_key %in% usable, .(gene_key, gene_type)], by = "gene_key")
gtype <- gtype[!gene_type %in% unlist(p$exclude_gene_types)]
gtype[, group := fifelse(gene_type %in% names(p$gene_type_fraction), gene_type, "other")]
n_type <- round(unlist(p$gene_type_fraction) * p$n_genes)
genes <- data.table(gene_key = unlist(lapply(names(n_type), function(t) draw(gtype[group == t]$gene_key, n_type[[t]]))))[sample(.N)]
for (t in names(n_type)) say("%s: %d of %d genes", t, n_type[[t]], gtype[group == t, .N])
genes[, `:=`(amount = rlnorm(.N, 0, p$sdlog), k = 1 + rpois(.N, p$lambda))]
# at most the isoforms annotated in every version
max_k <- rbindlist(lapply(ann, function(a) a[, .(n = .N), by = gene_key]))[
    , .(n = min(n)), by = gene_key]
genes[max_k, on = "gene_key", k := pmin(k, i.n)]

ranks <- genes[, .(rank = seq_len(k)), by = gene_key]
ranks[, share := {
    s <- rgamma(.N, p$alpha)
    s <- pmax(s / sum(s), p$min_share)
    sort(s / sum(s), decreasing = TRUE)
}, by = gene_key]
ranks[genes, on = "gene_key", control := i.amount * share]
say("%d genes, %d expressed isoforms (mean %.2f per gene)", nrow(genes), nrow(ranks), nrow(ranks) / nrow(genes))


# 2. treatment truth
genes[, class := "null"]
# dtu_swap eligibility as in Soneson et al. 2016
genes[, gene_tpm := amount / sum(amount) * 1e6]
genes[ranks[rank == 2], on = "gene_key", second_share := i.share]
pools <- list(
    dtu_swap = genes[k >= 2 & second_share > p$swap_min_share & gene_tpm >= p$swap_min_gene_tpm, gene_key],
    partial_dte = genes[k >= 2, gene_key],
    uniform = genes$gene_key
)
say("%d genes eligible for dtu_swap", length(pools$dtu_swap))
for (cl in names(pools)) {
    avail <- setdiff(pools[[cl]], genes[class != "null", gene_key])
    genes[gene_key %in% draw(avail, round(p$de_fraction[[cl]] * nrow(genes))), class := cl]
}
ranks[genes, on = "gene_key", class := i.class]
ranks[, target := FALSE]
ranks[class == "uniform", target := TRUE]
ranks[class == "dtu_swap", target := rank <= 2]
ranks[class == "partial_dte", target := seq_len(.N) %in% sample.int(.N, sample.int(min(p$max_targets, .N - 1), 1)),
      by = gene_key]

# direction chosen to keep the total RNA about constant
genes[class %in% c("uniform", "partial_dte"), abs_lfc := runif(.N, p$abs_log2fc[1], p$abs_log2fc[2])]
chg <- ranks[class %in% c("uniform", "partial_dte") & target, .(a = sum(control)), by = gene_key]
chg[genes, on = "gene_key", l := i.abs_lfc]
chg <- chg[sample.int(.N)]
running <- 0
chg[, sign := 0]
for (i in seq_len(nrow(chg))) {
    up <- chg$a[i] * (2^chg$l[i] - 1)
    down <- chg$a[i] * (2^-chg$l[i] - 1)
    s <- if (abs(running + up) <= abs(running + down)) 1 else -1
    set(chg, i, "sign", s)
    running <- running + if (s > 0) up else down
}
genes[chg, on = "gene_key", design_log2fc := i.l * i.sign]
ranks[genes, on = "gene_key", design_log2fc := i.design_log2fc]

ranks[, treatment := control]
ranks[class %in% c("uniform", "partial_dte") & target, treatment := control * 2^design_log2fc]
ranks[class == "dtu_swap" & target, treatment := rev(control), by = gene_key]
ranks[, true_log2fc := log2(treatment / control)]
ranks[, is_de := abs(true_log2fc) > 1e-8]
say("%s; %d up, %d down; %d DE isoforms; composition offset %.4f",
    paste(genes[, .N, by = class][, sprintf("%s=%d", class, N)], collapse = ", "),
    sum(chg$sign > 0), sum(chg$sign < 0), sum(ranks$is_de),
    log2(sum(ranks$treatment) / sum(ranks$control)))


# 3. replicates (genes x ranks)
samples <- data.table(condition = rep(c("control", "treatment"), each = p$n_reps), rep = rep(seq_len(p$n_reps), 2))
samples[, sample := sprintf("%s_rep%d", condition, rep)]
for (s in seq_len(nrow(samples))) {
    m <- samples$condition[s]
    g <- genes[, .(gene_key, b = rgamma(.N, shape = 1 / p$bcv^2, scale = p$bcv^2))]
    d <- ranks[, .(gene_key, a = get(m))]
    d[, share := a / sum(a), by = gene_key]
    d[, q := share]
    if (!is.null(p$share_precision)) {
        d[, q := { x <- rgamma(.N, shape = p$share_precision * share); x / sum(x) }, by = gene_key]
    }
    d[, total := sum(a), by = gene_key]
    ranks[, (samples$sample[s]) := d$total * g$b[match(d$gene_key, g$gene_key)] * d$q]
}


# per version: ranks to transcripts, write files
write_abundance <- function(tx, abund, path) {
    frag <- abund * tx$eff_length
    fwrite(data.table(target_id = tx$transcript_id, length = tx$length, eff_length = tx$eff_length,
                      est_counts = frag / sum(frag) * 1e7, tpm = abund / sum(abund) * 1e6),
           path, sep = "\t")
}
fwrite(genes, file.path(p$outdir, "de_design.tsv"), sep = "\t")
value_cols <- c("control", "treatment", samples$sample)
for (v in names(ann)) {
    tx <- ann[[v]]
    pick <- tx[gene_key %in% genes$gene_key]
    pick <- pick[, .(transcript_id = draw(transcript_id, genes$k[match(gene_key[1], genes$gene_key)])),
                 by = gene_key]
    pick[, rank := seq_len(.N), by = gene_key]
    pick <- merge(pick, ranks, by = c("gene_key", "rank"))
    tx <- merge(tx, pick[, c("transcript_id", "rank", "class", "target", value_cols, "design_log2fc",
                             "true_log2fc", "is_de"), with = FALSE], by = "transcript_id", all.x = TRUE)
    for (col in value_cols) set(tx, which(is.na(tx[[col]])), col, 0)

    d <- file.path(p$outdir, paste0("v", v))
    dir.create(d, showWarnings = FALSE)
    write_abundance(tx, tx$control, file.path(d, sprintf("control_abundance_v%s.tsv", v)))
    write_abundance(tx, tx$treatment, file.path(d, sprintf("treatment_abundance_v%s.tsv", v)))
    for (s in samples$sample) write_abundance(tx, tx[[s]], file.path(d, sprintf("%s_v%s.tsv", s, v)))
    fwrite(tx[!is.na(rank), .(transcript_id, gene_id, class, rank, target, true_log2fc, is_de)],
           file.path(d, sprintf("de_labels_v%s.tsv", v)), sep = "\t")
    say("v%s: %d transcripts in the annotation, %d expressed", v, nrow(tx), sum(tx$control > 0))
}
say("done: %s/v{%s}/", p$outdir, paste(p$versions, collapse = ","))
