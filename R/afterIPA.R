library(ggplot2)
library(reshape2)
library(plyr)
library(qvalue)
set.seed(15678563)

fdrCutoff <- 0.08  # maximum acceptable fdr in our analysis
logRatCutoff <- 1  # minimum log2(fold change) required to be declared significant

################################################################################
# Set up #######################################################################
################################################################################

# Load dataframes created by R/ratioAnalysis.R
load("temp/protTab.rda") # contains all protein info in wide format
load("temp/stat.rda") # contains some protein info in long format

# Theme for ggplot2 bar plots
mytheme <- theme_bw()+
    theme(legend.position = "none",
          legend.key = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black"),
          #panel.margin = unit(2, "mm"),
          legend.title = element_blank())


################################################################################
# Import individual IPA upstream regulator results #############################
################################################################################

# Files were exported from IPA. Each contains the genes in our analysis mapped to
# upstream regulators (either SMAD or TGFB) 
files <- list.files("data/IPA", full.names = T)
smadFiles <- grep("smad", files, value = T) # SMAD3 and SMAD7 are most interesting
tgfbFiles <- grep("tgfb", files, value = T)

# Reads IPA-exported file f and cleans columns to only have gene symbol and 
# upstream regulator.
readFilesIPA <- function(f) {
    dat <- read.delim(f, skip = 1, stringsAsFactors = F)
    dat$upReg <- gsub("^.*/", "", f)
    dat$upReg <- gsub(".txt$", "",  dat$upReg)
    return(data.frame(Symbol = dat[ , 1],
                      upReg = dat$upReg))
}

# Read SMAD and TGFB upstream regulator files
smad <- ldply(smadFiles, readFilesIPA)
tgfb <- ldply(tgfbFiles, readFilesIPA)

# Read Uniprot Accession <-> Gene Symbol key -----------------------------------
# The "Molecules" tab in IPA was exported for this analysis. It's useful 
# because it contains the precise mapping that IPA used for UniProt Accession
# to gene name (which is displayed in the analysis results).

key <- read.delim("data/IPA/key.txt", skip = 1, stringsAsFactors = F)
keep <- c("Symbol", grep("GenPept.UniProt", names(key), value = T))
key <- key[, keep]

key <- melt(key, id.vars = "Symbol", value.name = "accession")
key <- key[ , c("Symbol", "accession")]
key <- unique(key[key$accession != "--", ])

# Filter Protein Results -------------------------------------------------------

# This function merges and filters smad and tgfb dfs with the proteomics protein 
# info. Keeps only proteins quantified accross all 3 ages.
mergeTables <- function(df) {
    df <- merge(df, key, all = F)
    df <- merge(df, stat, all = F)
    df$upReg <- as.factor(df$upReg)
    
    df <- ddply(df, c("accession","upReg"), function(x) {
        # T when protein x is quantified in all 3 ages, F otherwise.
        if(nrow(x) >= 3) { x$k <- T } else { x$k <- F }
        return(x)
    })
    
    df <- df[df$k, ]
    return(df[ , names(df) != "k"])
}

# Filter SMAD and TGFB results
smadTab <- mergeTables(smad)
tgfbTab <- mergeTables(tgfb)

write.table(smadTab, "temp/smadTab.txt", sep = "\t", quote = F, row.names = F)
write.table(tgfbTab, "temp/tgfbTab.txt", sep = "\t", quote = F, row.names = F)



################################################################################
# Upstream Regulator Z-Score and p-value #######################################
################################################################################
# Read upstream regulator table exported form IPA.
regs <- read.delim("data/IPA/upReg.txt", skip = 1, stringsAsFactors = F)

# Filter for TGFB and SMAD upstream regulators, specifically
tgfbRegs <- regs[grep("TGFB[123]$", regs$Upstream.regulators), ]
smadRegs <- regs[grep("^SMAD[347]$", regs$Upstream.regulators), ]
otherRegs <- regs[grep("(PDGF BB|^TP53$|APOE)", regs$Upstream.regulators), ]

# Function Takes an upstream regulator list and puts in long format, based on age.
# Finishes by cleaning up the age column.
meltRegs <- function(dfRegs) {
    df <- dfRegs[ , 1:(ncol(dfRegs)-1)]
    df <- melt(df, id.vars = "Upstream.regulators", 
               value.name = "z.score", 
               variable.name = "age")
    df$age <- gsub("^X", "", df$age)
    return(df)
}

# Melt tgfbReg and smadReg
tgfbRegT <- meltRegs(tgfbRegs)
tgfbRegT$z.score <- as.numeric(tgfbRegT$z.score)

smadRegT <- meltRegs(smadRegs)
smadRegT$z.score <- as.numeric(smadRegT$z.score)

otherRegT <- meltRegs(otherRegs)
otherRegT$z.score <- as.numeric(otherRegT$z.score)

# Make pretty plots of Upstream Regulator Z-Scores -----------------------------

# TGFB
lvl <- levels(factor(tgfbRegT$Upstream.regulators))
annTxt <- data.frame(age = "12d", 
                     z.score = 0, 
                     Upstream.regulators = factor("TGFB2", levels = lvl))

ggplot(tgfbRegT, aes(x = age, y = z.score, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.75) +
    mytheme +
    geom_hline(yintercept = 0, size = 0.5, color = "black") + 
    facet_grid(. ~ Upstream.regulators) +
    ylab("Activation Z-Score") +
    xlab("Age") +
    geom_text(data = annTxt, label = "ND", vjust = -1, size = 2)
    
ggsave("results/tgfbRegs.pdf", width = 87, height = 40, units = "mm", useDingbats = F)
ggsave("results/tgfbRegs.tiff", width = 87, height = 40, units = "mm")


# SMAD
lvl <- levels(factor(smadRegT$Upstream.regulators))
annTxt <- data.frame(age = "12d", 
                     z.score = 0, 
                     Upstream.regulators = factor("SMAD4", levels = lvl))

ggplot(smadRegT, aes(x = age, y = z.score, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.75) +
    mytheme +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    facet_grid(. ~ Upstream.regulators) +
    ylab("Activation Z-Score") +
    xlab("Age") +
    geom_text(data = annTxt, label = "ND", vjust = -1, size = 2)

ggsave("results/smadRegs.pdf", width = 87, height = 40, units = "mm", useDingbats = F)
ggsave("results/smadRegs.tiff", width = 87, height = 40, units = "mm")

#APOE, PDGF and TP53
ggplot(otherRegT, aes(x = age, y = z.score, fill = age)) +
  geom_bar(stat = "identity", color = "black", width = 0.75) +
  mytheme +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  facet_grid(. ~ Upstream.regulators) +
  ylab("Activation Z-Score") +
  xlab("Age") 

ggsave("results/otherRegs.pdf", width = 87, height = 40, units = "mm", useDingbats = F)


################################################################################
# Individual Protein Plots #####################################################
################################################################################

# TGFB2 --------------------------------------------------------------
# Select TGFB2
tgfb2 <- protTab$accession[grep("Tgfb2", protTab[ , "Gene.names"])]
tgfb2Row <- stat[stat$accession == tgfb2, ]

# plot
ggplot(tgfb2Row, aes(x = age, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    theme_bw()+
    theme(legend.position = "none",
          legend.key = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black")) +
    ylim(c(-0.65, 2.5)) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("TGFB2")

ggsave("results/tgfb2.pdf", width = 30, height = 40, units = "mm", useDingbats = F)
ggsave("results/tgfb2.tiff", width = 30, height = 40, units = "mm")

# LTBP2 ------------------------------------------------------------------------
ltbp2 <- protTab$accession[grep("Ltbp2", protTab[ , "Gene.names"])]
ltbp2Row <- stat[stat$accession == ltbp2, ]

ggplot(ltbp2Row, aes(x = age, y = logFC, fill = age)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5) + 
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  theme_bw()+
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black")) +
  ylim(c(-0.8, 3)) +
  ylab(expression("Log"[2]~"KO/wt Ratio")) +
  xlab("LTBP2")

ggsave("results/ltbp2.pdf", width = 30, height = 40, units = "mm", useDingbats = F)
ggsave("results/ltbp2.tiff", width = 30, height = 40, units = "mm")

# LTBP4 ------------------------------------------------------------------------
ltbp4 <- protTab$accession[grep("Ltbp4", protTab[ , "Gene.names"])]
ltbp4Row <- stat[stat$accession == ltbp4[2], ]

ggplot(ltbp4Row, aes(x = age, y = logFC, fill = age)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5) + 
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  theme_bw()+
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black")) +
  #ylim(c(-0.65, 2.5)) +
  ylab(expression("Log"[2]~"KO/wt Ratio")) +
  xlab("LTBP4")

ggsave("results/ltbp4.pdf", width = 30, height = 40, units = "mm", useDingbats = F)
ggsave("results/ltbp4.tiff", width = 30, height = 40, units = "mm")

################################################################################
# Panels of Protein Plots ######################################################
################################################################################

ltbp2 <- protTab$accession[grep("Ltbp2", protTab[ , "Gene.names"])]
ltbp2Row <- stat[stat$accession == ltbp2, ]

ltbp4 <- protTab$accession[grep("Ltbp4", protTab[ , "Gene.names"])]
ltbp4Row <- stat[stat$accession == ltbp4[2], ]

ltbp <- rbind(tgfb2Row, ltbp2Row, ltbp4Row)
ltbp$lab <- factor(c(rep("Tgfb2",3), rep("Ltbp2",3), rep("Ltbp4",3)), 
                   levels = c("Tgfb2", "Ltbp2", "Ltbp4"))

dodge <- position_dodge(width=0.7)
ggplot(ltbp, aes(x = lab, y = logFC, fill = age)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  mytheme +
  #ylim(c(-3.5, 8)) +
  ylab(expression("Log"[2]~"KO/wt Ratio")) +
  xlab("Protein") +
  theme(axis.title.x = element_blank())

ggsave("results/ltbp.pdf", width = 60, height = 40, units = "mm", useDingbats = F)
ggsave("results/ltbp.tiff", width = 60, height = 40, units = "mm")


# Differentially expressed in all heatmap --------------------------------------
protTab2 <- ddply(protTab, "accession", function(p) {
  p$allSig <- all((p$fdr.ST_12d <= fdrCutoff & p$logFC_12d^2 >= logRatCutoff^2),
                  (p$fdr.ST_15w <= fdrCutoff & p$logFC_15w^2 >= logRatCutoff^2),
                  (p$fdr.ST_1y <= fdrCutoff & p$logFC_15w^2 >= logRatCutoff^2), 
                  na.rm = F)
  return(p)
})

protTab2 <- protTab2[protTab2$allSig & !is.na(protTab2$allSig), ]

allSig <- stat[stat$accession %in% protTab2$accession, ]
allSig$gn <- mapvalues(allSig$accession, 
                       from = protTab2$accession, 
                       to = protTab2$Gene.names)

write.table(protTab2, file = "temp/alldiff.txt", sep = "\t", row.names = F, quote = F)

# Categories
protease <- c("A0A0R4J0I9", "Q54AE5", "Q8R054", "Q9R118")
stressResponse <- c("A2A813", "P08228", "A0A0R4J139")
glycan <- c("A6H6K1", "A6MDD3", "Q8C253", "Q9CQ60", "Q9R045", "Q60675")
nucMet <- c("P08030", "Q4FK28", "Q9R0Y5")
signaling <- c("P29268","Q3TX21", "Q5EBQ2", "Q9JJU8")
cytoSkel <- c("Q4FK36", "Q9CRB6")

# Proteases
dodge <- position_dodge(width=0.7)
ggplot(allSig[allSig$accession %in% protease, ], aes(x = gn, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    mytheme +
    #ylim(c(-3.5, )) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("Protein") +
    theme(axis.title.x = element_blank())

ggsave("results/proteaseBar.pdf", height = 60, width = 100, units = "mm", useDingbats = F)

# Stress Response
ggplot(allSig[allSig$accession %in% stressResponse, ], aes(x = gn, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    mytheme +
    ylim(c(-3.75, 0)) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("Protein") +
    theme(axis.title.x = element_blank())

ggsave("results/stressBar.pdf", height = 30, width = 80, units = "mm", useDingbats = F)

# Glycans and ECM
ggplot(allSig[allSig$accession %in% glycan, ], aes(x = gn, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    mytheme +
    ylim(c(-2.5, 3.75)) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("Protein") +
    theme(axis.title.x = element_blank())

ggsave("results/ecmBar.pdf", height = 60, width = 150, units = "mm", useDingbats = F)

# Nucleotide Metabolism
ggplot(allSig[allSig$accession %in% nucMet, ], aes(x = gn, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    mytheme +
    ylim(c(-3.25, 0)) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("Protein") +
    theme(axis.title.x = element_blank())

ggsave("results/nucleotideBar.pdf", height = 30, width = 80, units = "mm", useDingbats = F)

# Signaling
ggplot(allSig[allSig$accession %in% signaling, ], aes(x = gn, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    mytheme +
    ylim(c(-3.25, 8)) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("Protein") +
    theme(axis.title.x = element_blank())

ggsave("results/signalingBar.pdf", height = 60, width = 100, units = "mm", useDingbats = F)

# Cytoskeleton
ggplot(allSig[allSig$accession %in% cytoSkel, ], aes(x = gn, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black", width = 0.7, position = dodge) +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.H), width = .2, size = 0.5, position = dodge) + 
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    mytheme +
    ylim(c(-5, 0)) +
    ylab(expression("Log"[2]~"KO/wt Ratio")) +
    xlab("Protein") +
    theme(axis.title.x = element_blank())

ggsave("results/cytoBar.pdf", height = 30, width = 60, units = "mm", useDingbats = F)
