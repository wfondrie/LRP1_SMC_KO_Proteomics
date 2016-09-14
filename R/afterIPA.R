library(ggplot2)
library(reshape2)
library(plyr)
library(qvalue)
set.seed(15678563)

################################################################################
# Set up #######################################################################
################################################################################

# Load dataframes created by R/ratioAnalysis.R
load("temp/protTab.rda") # contains all protein info in wide format
load("temp/stat.rda") # contains some protein info in long format

# Theme for ggplot2 bar plots
mytheme <- theme_bw()+
    theme(legend.key = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black"),
          legend.title = element_blank())


################################################################################
# Import individual IPA upstream regulator results #############################
################################################################################

# Files were exported from IPA. Each contains the genes in our analysis mapped to
# upstream regulators (either SMAD or TGFB) 
files <- list.files("data/IPA", full.names = T)
smadFiles <- grep("smad[37]", files, value = T) # SMAD3 and SMAD7 are most interesting
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
tgfbRegs <- regs[grep("TGFB[123]", regs$Upstream.regulators), ]
smadRegs <- regs[grep("^SMAD[37]$", regs$Upstream.regulators), ]

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

# Make pretty plots of Upstream Regulator Z-Scores -----------------------------

# TGFB
ggplot(tgfbRegT, aes(x = Upstream.regulators, y = z.score, fill = age)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    mytheme +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    ylab("Activation Z-Score") +
    xlab("Upstream Regulator")
    
ggsave("results/tgfbRegs.pdf", width = 60, height = 35, units = "mm", useDingbats = F)
ggsave("results/tgfbRegs.tiff", width = 60, height = 35, units = "mm")


# SMAD
ggplot(smadRegT, aes(x = Upstream.regulators, y = z.score, fill = age)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    mytheme +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    ylab("Activation Z-Score") +
    xlab("Upstream Regulator")

ggsave("results/smadRegs.pdf", width = 50, height = 35, units = "mm", useDingbats = F)
ggsave("results/smadRegs.tiff", width = 50, height = 35, units = "mm")


################################################################################
# Individual Protein Plots #####################################################
################################################################################

# TGFB2 ------------------------------------------------------------------------
# Select TGFB2
tgfb2 <- protTab$accession[grep("Tgfb2", protTab[ , "Gene.names"])]
tgfb2Row <- stat[stat$accession == tgfb2, ]

# plot
ggplot(tgfb2Row, aes(x = age, y = logFC, fill = age)) +
    geom_bar(stat = "identity", color = "black") +
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

ggsave("results/tgfb2.pdf", width = 35, height = 35, units = "mm", useDingbats = F)
ggsave("results/tgfb2.tiff", width = 35, height = 35, units = "mm")