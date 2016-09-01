library(ggplot2)
library(reshape2)
library(plyr)
library(qvalue)

# Import gene name <-> Uniprot accession key for IPA data ----------------------
key <- read.delim("data/IPA/genemap_HL.txt", skip = 1, stringsAsFactors = F)
key <- key[ , c("Symbol", grep("GenPept", names(key), value = T))]
key <- melt(key, id.vars = "Symbol", variable.name = "discard", value.name = "accession")
key <- key[ , c("Symbol", "accession")]
key <- unique(key)
key <- key[key$accession != "--", ]

# Import tgfb genes ------------------------------------------------------------
tgfb1 <- read.delim("data/IPA/tgfb1_HL.txt", skip = 1, stringsAsFactors = F)
tgfb2 <- read.delim("data/IPA/tgfb2_HL.txt", skip = 1, stringsAsFactors = F)
tgfb3 <- read.delim("data/IPA/tgfb3_HL.txt", skip = 1, stringsAsFactors = F)
smad7 <- read.delim("data/IPA/smad7_HL.txt", skip = 1, stringsAsFactors = F)

# Import Upstream Regulators ---------------------------------------------------
# read file with p.values for each upstream regulator
pvals <- read.delim("data/IPA/upstreamRegulators_pval_HL.txt", skip = 1, stringsAsFactors = F)
pvals <- pvals[ , 1:4]
names(pvals) <- c("regulator", "age_15w", "age_1y", "age_12d")
pvals <- melt(pvals, id.vars = "regulator", variable.name = "age", value.name = "log10_p.value")

# read file with Z-scores for each upstream regulator
zscores <- read.delim("data/IPA/upstreamRegulators_HL.txt", skip = 1, stringsAsFactors = F)
zscores <- zscores[ , 1:4]
names(zscores) <- c("regulator", "age_15w", "age_1y", "age_12d")
zscores <- melt(zscores, id.vars = "regulator", variable.name = "age", value.name = "z.score")

# Merge into one df
uReg <- merge(pvals, zscores, all = T)

uReg$log10_p.value[uReg$log10_p.value == "N/A"] <- NA
uReg$z.score[uReg$z.score == "N/A"] <- NA

# Calculate FDR for upstream regulators using Benjamini-Hochberg
uReg$log10_p.value <- as.numeric(uReg$log10_p.value)
uReg$z.score <- as.numeric(uReg$z.score)
uReg$p.value <- 10^(-uReg$log10_p.value)
uReg <- ddply(uReg, "age", function(p) {
  p$fdr <- p.adjust(p$p.value, method = "BH")
  p
})

# Pull out TGFB related regulators
tgfbReg <- uReg[grep("(tgfb|smad)", uReg$regulator, ignore.case = T), ]

# Look at all TGF-B related genes ----------------------------------------------

