suppressPackageStartupMessages({
    library(tximport)
    library(data.table)
})

quant_sf <- snakemake@input[["quant_sf"]]
output_path <- snakemake@output[[1]]
version <- snakemake@wildcards[["version"]]
sample <- snakemake@params[["sample"]]

# fishpond::computeInfRV defaults
PC <- 5
SHIFT <- 0.01

txi <- tximport(quant_sf, type = "salmon", txOut = TRUE)

reps <- txi$infReps[[1]]

mu <- rowMeans(reps)
inf_var <- apply(reps, 1, var)
inf_rv <- pmax(inf_var - mu, 0) / (mu + PC) + SHIFT

dt <- data.table(
    transcript_id = rownames(txi$counts),
    mean_count = mu,
    inf_var = inf_var,
    InfRV = inf_rv,
    platform = "salmon",
    sample = sample,
    gencode_version = version
)

fwrite(dt, output_path, sep = "\t")
