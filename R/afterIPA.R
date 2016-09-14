library(ggplot2)
library(reshape2)
library(plyr)
library(qvalue)
load("temp/protTab.rda")
load("temp/stat.rda")

# Theme for ggplot2 bar plots
mytheme <- theme_bw()+
           theme(legend.key = element_blank(),
                 legend.key.size = unit(0.5, "lines"),
                 text = element_text(size = 8),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(color = "black"),
                 legend.title = element_blank())


# Proteins ---------------------------------------------------------------------
files <- list.files("data/IPA", full.names = T)
smadFiles <- grep("smad[37]", files, value = T)
tgfbFiles <- grep("tgfb", files, value = T)



# reads file f and formats it nicely
readFilesIPA <- function(f) {
  dat <- read.delim(f, skip = 1, stringsAsFactors = F)
  dat$upReg <- gsub("^.*/", "", f)
  dat$upReg <- gsub(".txt$", "",  dat$upReg)
  return(data.frame(Symbol = dat[ , 1],
                    upReg = dat$upReg))
}

smad <- ldply(smadFiles, readFilesIPA)
tgfb <- ldply(tgfbFiles, readFilesIPA)

key <- read.delim("data/IPA/key.txt", skip = 1, stringsAsFactors = F)
keep <- c("Symbol", grep("GenPept.UniProt", names(key), value = T))
key <- key[, keep]

key <- melt(key, id.vars = "Symbol", value.name = "accession")
key <- key[ , c("Symbol", "accession")]
key <- unique(key[key$accession != "--", ])


# Merge and filter smad and tgfb dfs with the proteomics protein info
mergeTables <- function(df) {
  df <- merge(df, key, all = F)
  df <- merge(df, stat, all = F)
  df$upReg <- as.factor(df$upReg)
  df <- ddply(df, c("accession","upReg"), function(x) {
    if(nrow(x) >= 3) { x$k <- T } else { x$k <- F }
    return(x)
  })
  df <- df[df$k, ]
  return(df[ , names(df) != "k"])
}


smadTab <- mergeTables(smad)
tgfbTab <- mergeTables(tgfb)

write.table(smadTab, "temp/smadTab.txt", sep = "\t", quote = F, row.names = F)
write.table(tgfbTab, "temp/tgfbTab.txt", sep = "\t", quote = F, row.names = F)

# Regulators -------------------------------------------------------------------
regs <- read.delim("data/IPA/upReg.txt", skip = 1, stringsAsFactors = F)
tgfbRegs <- regs[grep("TGFB[123]", regs$Upstream.regulators), ]
smadRegs <- regs[grep("^SMAD[37]$", regs$Upstream.regulators), ]

meltRegs <- function(dfRegs) {
  df <- dfRegs[ , 1:(ncol(dfRegs)-1)]
  df <- melt(df, id.vars = "Upstream.regulators", 
             value.name = "z.score", 
             variable.name = "age")
  df$age <- gsub("^X", "", df$age)
  return(df)
}

tgfbRegT <- meltRegs(tgfbRegs)
tgfbRegT$z.score <- as.numeric(tgfbRegT$z.score)
ggplot(tgfbRegT, aes(x = Upstream.regulators, y = z.score, fill = age)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  mytheme +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  ylab("Activation Z-Score") +
  xlab("Upstream Regulator") + 

ggsave("results/tgfbRegs.pdf", width = 60, height = 35, units = "mm", useDingbats = F)
ggsave("results/tgfbRegs.tiff", width = 60, height = 35, units = "mm")

smadRegT <- meltRegs(smadRegs)
smadRegT$z.score <- as.numeric(smadRegT$z.score)

ggplot(smadRegT, aes(x = Upstream.regulators, y = z.score, fill = age)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  mytheme +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  ylab("Activation Z-Score") +
  xlab("Upstream Regulator")

ggsave("results/smadRegs.pdf", width = 50, height = 35, units = "mm", useDingbats = F)
ggsave("results/smadRegs.tiff", width = 50, height = 35, units = "mm")

# TGFB2 Plot -------------------------------------------------------------------

tgfb2 <- protTab$accession[grep("Tgfb2", protTab[ , "Gene.names"])]
tgfb2Row <- stat[stat$accession == tgfb2, ]

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