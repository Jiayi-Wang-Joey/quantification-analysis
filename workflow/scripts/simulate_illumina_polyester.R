suppressPackageStartupMessages({
    library(polyester)
    library(Biostrings)
    library(data.table)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
p <- snakemake@params
set.seed(p$seed)

ab <- fread(snakemake@input[["abundance"]], select = c("target_id", "est_counts"))
tx <- readDNAStringSet(snakemake@input[["transcriptome"]])
names(tx) <- sub("\\s.*", "", names(tx))

counts <- as.vector(rmultinom(1, p$depth, ab$est_counts / sum(ab$est_counts)))
keep <- counts > 0
ab <- ab[keep][, n := counts[keep]]
missing <- setdiff(ab$target_id, names(tx))
if (length(missing)) stop(length(missing), " transcripts with fragments are not in the transcriptome FASTA")
message(sprintf("%d fragments over %d transcripts", sum(ab$n), nrow(ab)))

n_chunks <- ceiling(sum(ab$n) / p$chunk_pairs)
ab[, chunk := floor((cumsum(n) - 1) * n_chunks / sum(n))]
r1 <- gzfile(snakemake@output[["r1"]], "wt")
r2 <- gzfile(snakemake@output[["r2"]], "wt")
write_fastq <- function(fa, con, chunk) {
    reads <- readDNAStringSet(fa)
    seqs <- as.character(reads)
    writeLines(paste0("@c", chunk, "_", names(reads), "\n", seqs, "\n+\n",
                      strrep("I", nchar(seqs))), con)
}
for (ch in unique(ab$chunk)) {
    d <- ab[chunk == ch]
    tmp <- file.path(p$tmpdir, paste0("chunk", ch))
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    fa <- file.path(tmp, "transcripts.fa")
    writeXStringSet(tx[d$target_id], fa)
    simulate_experiment_countmat(fasta = fa, readmat = matrix(d$n, ncol = 1), outdir = tmp,
                                 readlen = p$read_length, paired = TRUE, error_rate = p$error_rate,
                                 seed = p$seed + ch, gzip = TRUE)
    write_fastq(file.path(tmp, "sample_01_1.fasta.gz"), r1, ch)
    write_fastq(file.path(tmp, "sample_01_2.fasta.gz"), r2, ch)
    unlink(tmp, recursive = TRUE)
    message(sprintf("chunk %d: %d fragments", ch, sum(d$n)))
}
close(r1)
close(r2)
unlink(p$tmpdir, recursive = TRUE)
