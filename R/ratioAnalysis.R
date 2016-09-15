library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(qvalue)
library(ggrepel)
library(gplots)
library(qvalue)
library(data.table)
set.seed(45681)

# Parameters for Analysis
pValCutoff <- 0.01 # maximum acceptable p-value in our analysis
fdrCutoff <- 0.10  # maximum acceptable fdr in our analysis
logRatCutoff <- 1  # minimum log2(fold change) required to be declared significant
minRep <- 2        # minimum number of replicates quantified per condition
minPep <- 2        # minimum number of peptides quantified per run


################################################################################
# Import MaxQuant results ######################################################
################################################################################
prot <- read.delim("data/combined/txt/proteinGroups.txt", stringsAsFactors = F)
#prot <- read.delim("data/test.txt", stringsAsFactors = F)

# Remove Contaminants and Reverse Sequences
prot <- prot[!grepl("^(CON_|REV_)", prot$Protein.IDs), ]

# Need to rename columns to account for the bad 1y replicate that was thrown out
names(prot) <- gsub("1y_3", "1y_2", names(prot))
names(prot) <- gsub("1y_4", "1y_3", names(prot))

# Choose quantification columns
ratCols <- grep("^Ratio.H.L.normalized.1", names(prot), value = T)

# Remove proteins with less than minPep
pepCols <- grep("Peptides.1", names(prot), value = T)

# Remove ratios calculated from less than minPep peptides
prot <- ddply(prot, "Protein.IDs", function(p) {
    lowPep <- p[ , pepCols] < minPep
    p[ , ratCols[lowPep]] <- NaN
    return(p)
})


# Remove proteins with no quantified reps after filtering
prot <- prot[rowSums(prot[ , ratCols], na.rm = T) > 0, ]

# add a simplified accession column
prot$accession <- gsub(";.*$","",prot$Protein.IDs)


################################################################################
# H/L Differential Expression ##################################################
################################################################################
# Prepare expression matrix for limma
qProt <- prot[ , ratCols]
row.names(qProt) <- prot$accession
qProt <- as.matrix(log2(qProt))

# Retrieve names for 3 age groups from column names
ages <- unique(gsub("_.*$", "",colnames(qProt)))


# Emprical Bayes moderated t-test, split by age. -------------------------------
### If all were together, only proteins quantified in every age would be analyzed.
stat <- ldply(ages, function(age) {
    cond <- qProt[ , grep(age, colnames(qProt))]
    n <- rowSums(!is.na(cond))
    cond <- cond[n >= minRep, ]
    fit <- lmFit(cond)
    testFit <- eBayes(fit)
    tab <- topTable(testFit, number = nrow(cond), sort.by = "none", confint = T)
    fdrST <- qvalue(tab$P.Value)$qvalues
    
    # Output dataframe
    rdat <- data.frame(accession = row.names(cond),
                       age = gsub("^.*\\.", "", age),
                       CI.L = tab$CI.L,
                       CI.H = tab$CI.R,
                       logFC = tab$logFC,
                       p.value = tab$P.Value,
                       fdr.BH = tab$adj.P.Val,
                       fdr.ST = fdrST, 
                       n = n[n >= minRep])
    rdat
})

# Create Pretty Volcano Plots --------------------------------------------------
vProt <- stat[complete.cases(stat), ]
vProt$sig <- vProt$fdr.ST <= fdrCutoff & vProt$logFC^2 >= logRatCutoff^2
maxP <- ddply(vProt[vProt$sig, ], "age", function(x) {max(x$p.value)})
names(maxP)[2] <- "cutoff"

# make the pretty plot
volcano <- ggplot(vProt, aes(x=logFC, y=-log10(p.value))) +
    geom_point(aes(color = sig), size = 0.5) +
    scale_color_manual(values=c("gray","black")) +
    geom_hline(data = maxP,aes(yintercept = -log10(cutoff)), size = 0.5, color="black") +
    facet_grid(. ~ age) +
    geom_vline(xintercept = c(logRatCutoff, -logRatCutoff), size = 0.5) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black"),
          strip.background = element_rect(color = "black")) +
    xlab(expression("Normalized log"[2]~"KO/wt ratio (H/L)")) +
    ylab(expression("-log"[10]~"p-value")) +
    xlim(-6.6,6.6)
volcano
ggsave("results/volcano.pdf", width = 70, height = 35, units = "mm", useDingbats = F)
ggsave("results/volcano.tiff", width = 70, height = 35, units = "mm")



# Data needed for Venn Diagram -------------------------------------------------
total <- length(unique(vProt$accession[vProt$sig]))
sig12d <- vProt[vProt$age == "12d" & vProt$sig, ]
sig15w <- vProt[vProt$age == "15w" & vProt$sig, ]
sig1y <- vProt[vProt$age == "1y" & vProt$sig, ]

allShared <- length(unique(vProt$accession[
  (vProt$accession %in% sig12d$accession) &
  (vProt$accession %in% sig15w$accession) &
  (vProt$accession %in% sig1y$accession)
]))

btw12d15w <- length(sig12d$accession[sig12d$accession %in% sig15w$accession]) - allShared
btw12d1y <- length(sig12d$accession[sig12d$accession %in% sig1y$accession]) - allShared
btw15w1y <- length(sig12d$accession[sig15w$accession %in% sig1y$accession]) - allShared

only12d <- nrow(sig12d) - (btw12d1y + btw12d15w + allShared)
only15w <- nrow(sig15w) - (btw12d15w + btw15w1y + allShared)
only1y <- nrow(sig1y) - (btw12d1y + btw15w1y + allShared)

vennDat <- data.frame(total = total,
                      sig12d = nrow(sig12d),
                      sig15w = nrow(sig15w),
                      sig1y = nrow(sig1y),
                      allShared = allShared,
                      shared_12d_15w = btw12d15w,
                      shared_12d_1y = btw12d1y,
                      shared_15w_1y = btw15w1y,
                      only12d = only12d,
                      only15w = only15w,
                      only1y = only1y)
vennDat$sanityCheck <- sum(vennDat[1, 5:ncol(vennDat)])

################################################################################
# Clustering ###################################################################
################################################################################
# > Clustering by PCA and heirarchical clustering to identify if age makes a 
#   difference.

# PCA --------------------------------------------------------------------------
# Select proteins quantified in ALL replicates.
cProt <- qProt[rowSums(is.finite(qProt)) == ncol(qProt), ]

# Calculate and extract principal components.
comps <- prcomp(cProt)
compVar <- summary(comps)$importance[2, ]
pcaPlot <- data.frame(comps$rotation)

# Rename columns so they look pretty as labels
pcaPlot$lab <- gsub("_"," \\(",row.names(pcaPlot))
pcaPlot$lab <- paste0(pcaPlot$lab, ")")
pcaPlot$lab <- gsub("^.*\\.", "", pcaPlot$lab)
pcaPlot$age <- gsub("_.*", "", row.names(pcaPlot))

# Make a pretty plot of first 2 components
ggplot(pcaPlot, aes(x = PC1, y = PC2, label = lab)) + 
    geom_point() + 
    geom_label_repel(aes(fill = age), 
                     force = 10, 
                     segment.color = "black",
                     size = 2)+
    geom_point(aes(color = age), size = 2) + 
    theme(legend.position = "none") +
    #xlim(c(0.23,0.40)) +
    #ylim(c(-0.515,0.5)) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black")) +
    xlab(paste0("PC 1 (", round(compVar[1]*100, 0),"%)")) +
    ylab(paste0("PC 2 (", round(compVar[2]*100, 0),"%)"))

ggsave("results/pca.tiff", width = 70, height = 60, units = "mm")
ggsave("results/pca.pdf", width = 70, height = 60, units = "mm", useDingbats = F)


# Heirarchical Clustering and Heatmap ------------------------------------------
# Create a dataframe for labels and names
nameDat <- data.frame(colName = colnames(cProt))

# Extract age from column name
nameDat$age <- gsub("_.*", "", colnames(cProt))
nameDat$age <- factor(gsub("^.*\\.", "", nameDat$age))

# Assign age to color corresponding to base ggplot2 scheme, so plots match
nameDat$color <- as.character(revalue(nameDat$age, c("12d" = "#F8766D", 
                                                     "15w" = "#00BA38", 
                                                     "1y" = "#619CFF")))

# Create pretty labels fror columns
nameDat$colName <- gsub("_"," \\(", nameDat$colName)
nameDat$colName <- paste0(nameDat$colName, ")")
nameDat$colName <- factor(gsub("^.*\\.", "", nameDat$colName))

# Want to use a Red -> White -> Blue colorblind friendly color scheme.
myCols <- colorRampPalette(c("#b2182b",
                             "#d6604d",
                             "#f4a582",
                             "#fddbc7",
                             "#f7f7f7",
                             "#d1e5f0",
                             "#92c5de",
                             "#4393c3",
                             "#2166ac"))(n=64)

# Change gradient scale to highligh biologically significant changes.
b <- c(min(cProt),seq(-3,3, length = 63), max(cProt))

# Make a pretty heatmap
pdf(file = "results/heatmap.pdf", 
    width = 50 * 0.0392701, 
    height = 100 * 0.0392701, 
    useDingbats = F)
heatmap.2(cProt,
          trace = "none",
          dendrogram = "col",
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          labRow = F,
          col = myCols,
          breaks = b,
          density.info = "none",
          denscol = "black",
          key.title = "",
          key.xlab = expression("Log"[2]~"KO/wt Ratio (H/L)"),
          ColSideColors = nameDat$color,
          margins = c(4,1),
          srtCol = 45,
          cexCol = 0.85,
          keysize = 1,
          key.par = list(mar = c(5.5,2,0,2)),
          lmat = rbind(c(0,4),c(0,1),c(3,2),c(0,5)),
          lhei = c(1,0.25,3,1.25),
          labCol = nameDat$colName)
graphics.off()


################################################################################
# Age Comparisons ##############################################################
################################################################################
# Only want genes quantified in minRep replicates in all age groups
keepProt <- ddply(stat, "accession", function(prot){
    data.frame(fullReps = (length(prot$n) == 3))
})

# fqProt means fully quantified proteins = proteins quantified in all 3 ages.
keepAcc <- as.character(keepProt$accession[keepProt$fullReps])
k <- rownames(qProt) %in% keepAcc
fqProt <- qProt[ k, ]

# Setting up for limma ---------------------------------------------------------
# runDat contains info on age and run of each ratio column in fqProt
runDat <- data.frame(name = colnames(fqProt))
runDat$run <- gsub("^.*\\.", "", runDat$name)
runAge <- gsub("_.*$", "", runDat$run)

# limma ------------------------------------------------------------------------
design <- model.matrix(~ 0 + runAge)
fit <- lmFit(fqProt, design)

# Do pairwise contrasts between age groups.
contrastMat <- makeContrasts(YvM = runAge15w - runAge12d,
                             YvO = runAge1y - runAge12d,
                             MvO = runAge1y - runAge15w,
                             levels = design)
contrastFit <- contrasts.fit(fit, contrastMat)
testFit <- eBayes(contrastFit)

# Keep F statistics
anovaTab <- topTable(testFit, number = nrow(fqProt), sort.by = "none", confint = T)
anovaTab$F.fdrST <- qvalue(anovaTab$P.Value)$qvalues
anovaTab$sig <- anovaTab$F.fdrST <= 0.01
anovaTab$accession <- row.names(anovaTab)

# Age paiwise comparisions
ageStat <- ldply(1:3, function(contrast) {
    tab <- topTable(testFit, 
                    coef = contrast, 
                    nrow(fqProt), 
                    sort.by = "none", 
                    confint = T)
    fdrST <- qvalue(tab$P.Value)$qvalues
    
    # Output dataframe
    data.frame(accession = row.names(tab),
               contrast = colnames(testFit$coefficients)[contrast],
               CI.L = tab$CI.L,
               CI.H = tab$CI.R,
               logFC = tab$logFC,
               p.value = tab$P.Value,
               fdr.BH = tab$adj.P.Val,
               fdr.ST = fdrST)
    
})


# Make a pretty volcano plot ---------------------------------------------------
ageStat$sig <- ageStat$fdr.ST <= fdrCutoff & ageStat$logFC^2 >= logRatCutoff^2
ageStat$contrast <- revalue(ageStat$contrast, c("YvM" = "12d vs 15w",
                                                "YvO" = "12d vs 1y",
                                                "MvO" = "15w vs 1y"))
maxP <- ddply(ageStat[ageStat$sig, ], "contrast", function(x) {max(x$p.value)})
names(maxP)[2] <- "cutoff"


ageVolcano <- ggplot(ageStat, aes(x=logFC, y=-log10(p.value))) +
    geom_point(aes(color = sig), size = 0.5) +
    scale_color_manual(values=c("gray","black")) +
    geom_hline(data = maxP,aes(yintercept = -log10(cutoff)), size = 0.5, color="black") +
    facet_grid(~ contrast) +
    geom_vline(xintercept = c(logRatCutoff, -logRatCutoff), size = 0.5) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black"),
          strip.background = element_rect(color = "black")) +
    xlab(expression("Difference of log"[2]~"KO/wt ratios (H/L)")) +    
    ylab(expression("-log"[10]~"p-value")) +
    scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4.71, 4.71))

ageVolcano
ggsave("results/ageVolcano.pdf", width = 70, height = 35, units = "mm", useDingbats = F)
ggsave("results/ageVolcano.tiff", width = 70, height = 35, units = "mm")

# Save statistics data for after IPA
save(stat, file = "temp/stat.rda")


# Data needed for Venn Diagram -------------------------------------------------
total <- length(unique(ageStat$accession[ageStat$sig]))
comp1 <- ageStat[ageStat$contrast == "12d vs 15w" & ageStat$sig, ]
comp2 <- ageStat[ageStat$contrast == "12d vs 1y" & ageStat$sig, ]
comp3 <- ageStat[ageStat$contrast == "15w vs 1y" & ageStat$sig, ]

allShared <- length(unique(ageStat$accession[
  (ageStat$accession %in% comp1$accession) &
  (ageStat$accession %in% comp2$accession) &
  (ageStat$accession %in% comp3$accession)
]))

btw12 <- length(comp1$accession[comp1$accession %in% comp2$accession]) - allShared
btw13 <- length(comp1$accession[comp1$accession %in% comp3$accession]) - allShared
btw23 <- length(comp2$accession[comp2$accession %in% comp3$accession]) - allShared

only1 <- nrow(comp1) - (btw12 + btw13 + allShared)
only2 <- nrow(comp2) - (btw12 + btw23 + allShared)
only3 <- nrow(comp3) - (btw13 + btw23 + allShared)

vennDatAge <- data.frame(total = total,
                         allShared = allShared,
                         shared_comp1_comp2 = btw12,
                         shared_comp1_comp3 = btw13,
                         shared_comp2_comp3 = btw23,
                         only_comp1 = only1,
                         only_comp2 = only2,
                         only_comp3 = only3)
vennDatAge$sanityCheck <- sum(vennDatAge[1, 2:ncol(vennDatAge)])

middle <- ageStat[
  (ageStat$accession %in% comp1$accession) &
  (ageStat$accession %in% comp2$accession) &
  (ageStat$accession %in% comp3$accession),
]

################################################################################
# Write Tables #################################################################
################################################################################
# Protein Table (protTab)
# - summarize all of the protein level results

# Merge all of our different statistics tables
stat <- data.table(stat)
ageStat <- data.table(ageStat)

cols <- c("CI.L", "CI.H", "logFC", "p.value", "fdr.BH", "fdr.ST")

ageTable <- data.frame(dcast(ageStat, accession ~ contrast, value.var = c(cols, "sig")))
ratTable <- data.frame(dcast(stat, accession ~ age, value.var = c(cols, "n")))

protTab <- merge(merge(merge(prot, ratTable, all = T), ageTable, all = T), anovaTab, all = T)

save(protTab, file = "temp/protTab.rda")

# Keep only the columns we're interested in
keep <- c("accession",
          "Gene.names",
          "Protein.IDs", 
          "Protein.names",
          "Peptides",
          grep("^Peptides\\.", names(protTab), value = T),
          "Q.value",
          "Score",
          grep("^Ratio.H.L.1", names(protTab), value = T),
          grep("^Sequence.coverage", names(protTab), value = T),
          grep("^CI", names(protTab), value = T),
          grep("logFC", names(protTab), value = T),
          grep("p.value", names(protTab), value = T),
          grep("fdr.ST", names(protTab), value = T),
          "F",
          "P.Value",
          "F.fdrST")

protTab <- protTab[ , keep]

# Rename a few to clean it up
names(protTab) <- gsub("\\.+$", "", names(protTab))
names(protTab) <- gsub("ST", "", names(protTab))
names(protTab)[names(protTab) == "P.Value"] <- "F p-value"
names(protTab) <- gsub("p.value", "p-value", names(protTab))
names(protTab) <- gsub("Ratio.H\\.L", "H/L", names(protTab))
names(protTab) <- gsub("\\._"," ", names(protTab))
names(protTab) <- gsub("CI\\.", "CI-", names(protTab))
names(protTab) <- gsub("(\\.|_)", " ", names(protTab))

write.table(protTab, "results/proteinTable.txt", sep = "\t", row.names = F, quote = F)