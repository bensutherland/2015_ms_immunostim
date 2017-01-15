# Full analysis of the Ssal/Lsal immunostim project (v1.0)
# Overview:
#   01. Louse transcriptome analysis
#   02. Louse qPCR analysis
#   03. Salmon physiology and lice counts
#   04. Salmon qPCR analysis

# clean the space (remove #)
# rm(list=ls())

#
#######  01 - LICE TXOME ANALYSIS  ######
source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
library(limma)

# set working directory to the git repo main directory
setwd("~/Documents/koop/immunostim/03_analysis/02_lice_txome_analysis")


## input data ##
# define targets file (interpretation file for microarray data)
targets <- readTargets("Targets_2color_immunostim.csv", 
                       # path = "./02_raw_data/lice_txome_data/",
                       row.names="Name", sep = ",")

##read in datafiles (individual files are block and channel separated)
RG <- read.maimages(path = 
  "/Users/wayne/Documents/koop/immunostim/03_analysis/00_archive/louse-array-raw_data_and_sep/sep",
    files = targets[,c("FileNameCy3","FileNameCy5")], 
                    source="imagene",
                    columns = list(f="Signal Median", 
                                   b="Background Median"))
dim(RG$R) # number of probes and number samples
targets # see details on the samples

# view density plot of unprocessed data
plotDensities(RG, log=F)

#set control status of probes
spottypes <- readSpotTypes(file = "SpotTypes.txt")
RG$genes$Status <- controlStatus(types = spottypes, genes = RG)


## quality control plots ##
# plotMA3by2(RG) #MA plots per array saved to working directory
# imageplot3by2(RG) #imageplots per array saved to working directory

# boxplots of raw foreground/background for each channel
par(mfrow=c(2,2), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

units <- c("R", "G", "Rb", "Gb")
for(i in 1:4) {
  boxplot(log2(RG[[units[i]]]), xaxt = "n", ylab = "log2(fluor. units)", xlab = "samples", main = units[i])
}


## Data pre-processing ##
par(mfrow=c(3,1), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
## set low expressed genes statement in $RG$genes$isexpr column
RG$genes$isexpr.R <- rowSums(RG$R >= 500) >= 3 #here, value >= 500 in at least 3 arrays
                      ##isexpr > 500 in at least 3 arrays, lowest sample size = 4 for ld.
table(RG$genes$isexpr.R)
RG$genes$isexpr.G <- rowSums(RG$G >= 500) >= 3 ##isexpr > 500 in at least 3 arrays
table(RG$genes$isexpr.G)

## identify saturated probes
RG$genes$satspot.R <- rowSums(RG$R == 65535) >= 11 #saturated in all arrays
table(RG$genes$satspot.R) 
RG$genes$satspot.G <- rowSums(RG$G == 65535) >= 11
table(RG$genes$satspot.G) 

## background correction
MA.bm <- backgroundCorrect(RG, method="minimum") 
plotDensities(MA.bm)
## within array normalization (loess)
MA.bmW <- normalizeWithinArrays(MA.bm, method="loess", weights=NULL)
plotDensities(MA.bmW)

## between array normalization using Cy3 channel (Gquantile)
MA.bmWG <- normalizeBetweenArrays(MA.bmW, method="Gquantile")
plotDensities(MA.bmWG) 

## only retain custom probes, no control spots
MA.flt1 <- MA.bmWG[MA.bmWG$genes$Status == "Unknown",]
dim(MA.flt1) #38132 rows; 11 columns

## read in annotation file
annot <- read.table('koop_Lsal38K_2010-02-25b-short.txt',
                    header=TRUE, sep = "\t") #note sep = "\t" is important when empty cells
annot.df <- data.frame(annot)
dim(annot.df) #38132 rows; 5 columns - note: this should match the length of MA.flt1 given above

## attach annotation to the MA list
MA.flt1$genes$spotID <- annot.df[,1] 
MA.flt1$genes$description <- annot.df[,2] 
MA.flt1$genes$GeneID <- annot.df[,4] 

## Remove low expressed and saturated probes
MA.flt <- MA.flt1[MA.flt1$genes$isexpr.R == "TRUE" & MA.flt1$genes$isexpr.G == "TRUE"
                  & MA.flt1$genes$satspot.R != "TRUE" & MA.flt1$genes$satspot.G != "TRUE",]

## Save out full MAlist
write.csv(cbind(MA.flt$genes,MA.flt$M), "all_expr_data.csv", row.names=FALSE,
            #header=TRUE
          )


### ASIDE ### this dataframe will be required for correlation with qPCR in section 2 ###
# and as such, obtain probes of interest
probes.of.interest <- c("C042R126","C088R114","C263R087")
probes.for.cor <- as.data.frame(cbind(MA.flt$M[MA.flt$genes[12] == probes.of.interest[1]], MA.flt$M[MA.flt$genes[12] == probes.of.interest[2]] ,
      MA.flt$M[MA.flt$genes[12] == probes.of.interest[3]]))
rownames(probes.for.cor) <- colnames(MA.flt$M)
colnames(probes.for.cor) <- probes.of.interest
probes.for.cor
str(probes.for.cor)

## Save out probes.for.cor to the qPCR directory:
saveRDS(probes.for.cor, file = "~/Documents/koop/immunostim/03_analysis/03_lice_qPCR_analysis/probes.for.cor")

## Save out background list for GO enrichment testing
write.csv(MA.flt$genes, file = "background_immunostim.csv") #this can be used as the 'background list' for enrichment analysis

plotDensities(MA.flt) ##densitiy plot of expression values


### Differential Expression Analysis ## 
design <- modelMatrix(targets, ref = "ref")
attributes(design)
design

fit <- lmFit(object = MA.flt, design = design)
cont.matrix <- makeContrasts(ld-con, levels=design) ##makes contrasts
fit2 <- contrasts.fit(fit, cont.matrix) ###computes contrasts from lmod fit (outputs)
fit3 <- eBayes(fit2) ##given related parameter estim and SE, will comute mod t-stat, etc
output <- topTable(fit3, number=45000, adjust="none", p.value=0.05)
write.table(output, file="results.txt", sep="\t")
output.df <- data.frame(output)

## END ANALYSIS OF LOUSE TXOME



### CHECKED UP TO HERE ######

#
#######  02 - LICE QPCR  ########

# clean the space (remove #)
# rm(list=ls())

setwd("~/Documents/koop/immunostim/03_analysis/03_lice_qPCR_analysis")

# install required packages
source("http://bioconductor.org/biocLite.R")
biocLite("ReadqPCR")
biocLite("NormqPCR") # do not load from source requiring compilation.
library(ReadqPCR) # load the ReadqPCR library
library(NormqPCR)

immuno.data <- read.qPCR("new-together_formatted-notdouble.txt")

# inspection of batch object
dim(immuno.data) # 30 features, 11 samples
rownames(exprs(immuno.data)) # rows are genes (incl reps)
colnames(exprs(immuno.data)) # cols are samples
exprs(immuno.data) #lets you observe the raw Cts - note: there are two NA values for tech reps
which(is.na(exprs(immuno.data))) # currently two NAs

# combine technical replicates (and remove NAs)
combinedTechReps <- combineTechReps(immuno.data)
combinedTechReps #this combined gene names
exprs(combinedTechReps)  # gives Cts
featureNames(combinedTechReps) # gives gene names
which(is.na(exprs(combinedTechReps))) # in averaging, NAs are ignored; thus NAs removed

# remove genes that are not to be used in the analysis
row.names(combinedTechReps) # gives gene names
combinedTechReps <- combinedTechReps[-1,] # remove abcc1
row.names(combinedTechReps)

# geNORM
normalizers <- selectHKs(combinedTechReps, method = "geNorm",
                    Symbols = featureNames(combinedTechReps),
                    minNrHK = 2, log = T) #selectHKs is run with log = TRUE since we are using Cq
ranks <- data.frame(c(1, 1:13), normalizers$ranking)

# dCt
hkgs <- c("flna", "rps20")
dCt <- deltaCq(qPCRBatch =  combinedTechReps, hkgs = hkgs, calc="geo") # subtracts hkg values
head(exprs(dCt)) # in dCt format (i.e. still a log2 value)
# run statistics on log2 values, then transform to linear if desired

# create a dataframe using the dCt expression data, rows as samples, columns as genes
results <- as.data.frame(t(exprs(dCt)))
head(results)

### correlate microarray vs qPCR
## first obtain data from microarray
probes.for.cor <- readRDS(file = "~/Documents/koop/immunostim/03_analysis/03_lice_qPCR_analysis/probes.for.cor")
targets <- readTargets("~/Documents/koop/immunostim/03_analysis/02_lice_txome_analysis/Targets_2color_immunostim.csv", row.names="Name", sep = ",")

# load libraries for manipulating strings below:
library(stringi)
library(tidyr)

# split rownames into Cy5 and Cy3
probes.for.cor$sample <- rownames(probes.for.cor)
# provide a match column to match with qPCR targets file
probes.for.cor <- separate(data = probes.for.cor, col = sample, into = c("Cy3", "match"), sep = "\\.")
# provide a match column to match with probe file from array
targets$match <- stri_replace_all_fixed(str = targets$FileNameCy5,pattern = ".txt", replacement = "")

# need to pull in the targets file to obtain the sample ID's necessary to match with the qPCR sample IDs
probes.for.cor <- merge(x = targets, y = probes.for.cor, by = "match")
# remove extraneous columns
probes.for.cor.trimmed <- probes.for.cor[,c(2,1,6,9:11)]
# attach match names
probes.for.cor.trimmed$Name2match <- paste("Sample_",probes.for.cor.trimmed$Name, sep ="")
# order the two dataframes together
probes.for.cor.trimmed <- probes.for.cor.trimmed[with(probes.for.cor.trimmed, order(Name2match)), ]
probes.for.cor.trimmed

# in order:
# C042R126 Esterase SG1
# C088R114  Nuclear pore membrane glycoprotein 210
# C263R087 Histone-lysine N-methyltransferase SETD7

# models first to obtain R-squared
sg1 <- lm(probes.for.cor.trimmed$C042R126 ~ results$sg1)
pom210 <- lm(probes.for.cor.trimmed$C088R114 ~ results$pom210)
setd7 <- lm(probes.for.cor.trimmed$C263R087 ~ results$setd7)

sg1.adj.r.squared <- summary(sg1)$adj.r.squared
pom210.adj.r.squared <- summary(pom210)$adj.r.squared
setd7.adj.r.squared <- summary(setd7)$adj.r.squared

rsq <- c(sg1.adj.r.squared, pom210.adj.r.squared, setd7.adj.r.squared)
# names(rsq) <- c("sg1", "pom210", "setd7")
rsq

# plot the correlation
ylim.plot = c(-8,4)
xlim.plot = c(-8,4)
xtxt = -5.5
ytxt = 3.7
par(mfrow=c(1,3), mar = c(4,4,1,2))
plot(probes.for.cor.trimmed$C042R126 ~ {median(results$sg1)- results$sg1}, xlab = "ddCt sg1", ylab = "microarray log2(sg1)", las = 1
     , ylim = ylim.plot, xlim = xlim.plot
     )
text(x = xtxt, y = ytxt, labels = paste("R2 =", round(rsq[1], digits = 2)))
plot(probes.for.cor.trimmed$C088R114 ~ {median(results$pom21)- results$pom210}, xlab = "ddCt pom210", ylab = "microarray log2(pom210)", las = 1 
     , ylim = ylim.plot, xlim = xlim.plot
     )
text(x = xtxt, y = ytxt, labels = paste("R2 =", round(rsq[2], digits = 2)))
plot(probes.for.cor.trimmed$C263R087 ~ {median(results$setd7)- results$setd7}, xlab = "ddCt setd7", ylab = "microarray log2(setd7)", las = 1
     , ylim = ylim.plot, xlim = xlim.plot
     )
text(x = xtxt, y = ytxt, labels = paste("R2 =", round(rsq[3], digits = 2)))


### follows is ddCt calculation within NormqPCR, and significance testing among groups in qPCR (not really what is necessary)

# ## calculating ddCt within NormqPCR
# #need a contrast matrix (similar to limma, needs 0s and 1s (control and case))
# contM <- cbind(c(0,1,0,1,1,0,0,1,0,0,0),c(1,0,1,0,0,1,1,0,1,1,1))
# colnames(contM) <- c("low.dose","control")
# rownames(contM) <- colnames(combinedTechReps)
# contM
# hkgs<-c("flna", "rps20")
# ddCq.out <- deltaDeltaCq(qPCRBatch = combinedTechReps, 
#                         hkgs=hkgs, contrastM=contM, case="low.dose", 
#                         control="control", statCalc="arith", hkgCalc="geom")
# ddCq.out

# # Block G:  Plot of deltaDeltaCq values
# plotVals <- log2(as.numeric(as.character(ddCq.out[,6])))
# plotVals[as.character(ddCq.out[,6]) == "-                 "] <- min(plotVals, na.rm=TRUE)
# plotVals[as.character(ddCq.out[,6]) == "+                 "] <- max(plotVals, na.rm=TRUE)
# minusOneSD <- log2(as.numeric(as.character(ddCq.out[,7])))
# plusOneSD <- log2(as.numeric(as.character(ddCq.out[,8])))
# genes <- ddCq.out[,1]
# 
# colVec <- rep("blue", length(plotVals))
# colVec[as.character(ddCq.out[,6]) == "-                 " | as.character(ddCq.out[,6]) == "+                 "] <- "red"
# 
# library(gplots)
# barplot2(plotVals, col = colVec, ci.u = plusOneSD, ci.l = minusOneSD, plot.ci = TRUE, ylab="delta delta Cq", xlab = "genes")

## second method of analyzing:
# linear models of the values contained in 'results'
# {median(results$sg1)- results$sg1} #note the stats are the exact same whether you use dCt or ddCt

# dCt.df <- results
# dCt.df$group <- cbind(c("control","low-dose","control","low-dose","low-dose","control",
#                         "control","low-dose","control","control","control"))
# dCt.df
# sg1.mod <- aov(dCt.df$sg1
#   ~ dCt.df$group)
# summary(sg1.mod) # p = 0.04
# boxplot(2^(-dCt.df$sg1) ~ dCt.df$group, main = "sg1")
# 
# setd7.mod <- aov(dCt.df$setd7 ~ dCt.df$group)
# summary(setd7.mod) # p = 0.13
# boxplot(2^(-dCt.df$setd7) ~ dCt.df$group, main = "setd7 linear")
#                       
# pom210.mod <- aov(dCt.df$pom210 ~ dCt.df$group)
# summary(pom210.mod) # p = 0.08
# boxplot(2^(-dCt.df$pom210) ~ dCt.df$group, main = "pom210 linear")
# 
# 
# # re-check statistical significance of microarray data:
# t.test(x = probes.for.cor.trimmed$C042R126[probes.for.cor.trimmed$Cy5 == "con"], 
#        y = probes.for.cor.trimmed$C042R126[probes.for.cor.trimmed$Cy5 == "ld"])
# 
# 
# test <- aov(probes.for.cor.trimmed$C042R126 ~ probes.for.cor.trimmed$Cy5)
# summary(test)
# test2 <- aov(probes.for.cor.trimmed$C088R114 ~ probes.for.cor.trimmed$Cy5)
# summary(test2)
# test3 <- aov(probes.for.cor.trimmed$C263R087 ~ probes.for.cor.trimmed$Cy5)
# summary(test3)

#
#######  03 - LICE COUNTS & FISH PHYS   ############
# rm(list=ls())
setwd("~/Documents/koop/immunostim/03_analysis/04_salmon_physiol_analysis")

# import physiological and infection data
immuno.data <- read.csv(file = "immunostim_data-jan21-15.csv")
dim(immuno.data) # 313 rows, 15 columns
names(immuno.data)
immuno.data[1:10,]
levels(as.factor(immuno.data$time)) # confirms merge of time 3.1 and 3

# data setup, subset to remove the 'UNKN' treatment (tank control)
immuno.data <- subset(immuno.data, immuno.data$treatment != "UNKN")
immuno.data <- droplevels(immuno.data)
levels(immuno.data$treatment) # confirm the tank control has now been removed
immuno.data$treatment <- factor(immuno.data$treatment, 
                                levels = c("cont","LD","HD")) #relevel to put the right order (control, lowdose, highdose)

# weight at T(0) among tanks (confirm no difference at start) (use treatment groups instead of tanks for sample size)
# there are not enough replicates at T0 for independent tank evaluation, so do it with treatment groups
boxplot(immuno.data$weight.g[immuno.data$time==0] ~ 
         immuno.data$treatment[immuno.data$time==0], 
        na.rm =T, col = c("white","lightgrey","darkgrey"),
        xlab = "Tank ID", ylab = "Weight (g)")
points(immuno.data$weight.g[immuno.data$time==0] ~ immuno.data$treatment[immuno.data$time==0]) # n = 6 (total = 18 samples)
# model
pre.inf.mod = aov(immuno.data$weight.g[immuno.data$time==0] ~ immuno.data$treatment[immuno.data$time==0])
summary(pre.inf.mod) # no significant difference prior to infection (ANOVA p = 0.3)

# weight after start of feed groups
boxplot(immuno.data$weight.g[immuno.data$time!=0] ~ immuno.data$treatment[immuno.data$time!=0] * as.factor(immuno.data$time[immuno.data$time!=0])
        , col = c("white","lightgrey","darkgrey"))
weight.mod <- aov(immuno.data$weight.g[immuno.data$time!=0] ~ immuno.data$treatment[immuno.data$time!=0] * as.factor(immuno.data$time[immuno.data$time!=0]))
summary(weight.mod) # an effect of treatment ( p = 0.02)
TukeyHSD(weight.mod)

# note that weight data was not collected for indiv. 13-22 from tanks 6, 7, 8, 9, 10 (T4)
summary(is.na(immuno.data$weight.g[immuno.data$time==4])) # due to this, variable will be treated as LPF

# test without time 0 or 4
# weight after start of feed groups
boxplot(immuno.data$weight.g[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] ~ immuno.data$treatment[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] * 
          as.factor(immuno.data$time[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3])
        , col = c("white","lightgrey","darkgrey"))
weight.mod.noT0orT4 <- aov(immuno.data$weight.g[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] ~ immuno.data$treatment[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] * 
                             as.factor(immuno.data$time[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3]))
summary(weight.mod.noT0orT4) # an effect of treatment ( p = 0.02)
TukeyHSD(weight.mod.noT0orT4)




# Generate figure of lice per fish (T2-T4)
par(mar=c(2,2,2,2))
boxplot(immuno.data$per.fish.total[immuno.data$time > 1] ~ immuno.data$treatment[immuno.data$time > 1]
        * as.factor(immuno.data$time[immuno.data$time > 1]), na.rm =T, xlab = "Treatment",
        las = 1, ylab = "Lice per fish", xaxt = "n", 
        col = rep(c("white", "lightgrey", "darkgrey"), times = 3), ylim = c(0,12))
axis(1, at = seq(1:9), labels = c("CTR.T2", "LD.T2", "HD.T2", "CTR.T3", "LD.T3", "HD.T3", "CTR.T4", "LD.T4", "HD.T4"))
abline(v = c(3.5, 6.5), lty = 1)
text(x = seq(1:9), y = c(4,4,4, 9,9,9, 11.8,11.8,11.8), 
     label= c("a","a","a", "b","c","bc", "d","e","de"))

lpf.mod <- aov(immuno.data$per.fish.total[immuno.data$time > 1] ~ immuno.data$treatment[immuno.data$time > 1]
                    * as.factor(immuno.data$time[immuno.data$time > 1]))
anova(lpf.mod) 
#effect of treatment (p=6e-5); time (p = 3.5e-5); no interaction (p=0.4)
TukeyHSD(lpf.mod) # can do pair-wise comparisons to see within a time point which values are different

# fold change bw cont and ld
mean(immuno.data$per.fish.total[immuno.data$time == 3 & immuno.data$treatment ==  "cont"])
mean(immuno.data$per.fish.total[immuno.data$time == 3 & immuno.data$treatment ==  "LD"])
1/(1.4/2.9) # 2.1-fold lower at time 3
mean(immuno.data$per.fish.total[immuno.data$time == 4 & immuno.data$treatment ==  "cont"])
mean(immuno.data$per.fish.total[immuno.data$time == 4 & immuno.data$treatment ==  "LD"])
1/(1.56/3.48) # 2.2-fold lower at time 4


## are there significant differences in sex of lice between treatment groups?
names(immuno.data)

# collect proportions
LD.ad.m.4 <- sum(immuno.data$body.ad.m[immuno.data$time==4 & immuno.data$treatment=="LD"], na.rm =T)
LD.ad.fem.4 <-  sum(immuno.data$body.ad.fem[immuno.data$time==4 & immuno.data$treatment=="LD"], na.rm =T)
cont.ad.m.4 <- sum(immuno.data$body.ad.m[immuno.data$time==4 & immuno.data$treatment=="cont"], na.rm =T)
cont.ad.fem.4 <- sum(immuno.data$body.ad.fem[immuno.data$time==4 & immuno.data$treatment=="cont"], na.rm =T)
c(LD.ad.m.4, LD.ad.fem.4, cont.ad.m.4, cont.ad.fem.4)
# at time 4, more males than females in the LD, and opposite in the control

LD.ad.m.3 <- sum(immuno.data$body.ad.m[immuno.data$time==3 & immuno.data$treatment=="LD"], na.rm =T)
LD.ad.fem.3 <- sum(immuno.data$body.ad.fem[immuno.data$time==3 & immuno.data$treatment=="LD"], na.rm =T)
cont.ad.m.3 <- sum(immuno.data$body.ad.m[immuno.data$time==3 & immuno.data$treatment=="cont"], na.rm =T)
cont.ad.fem.3 <- sum(immuno.data$body.ad.fem[immuno.data$time==3 & immuno.data$treatment=="cont"], na.rm =T)
c(LD.ad.m.3, LD.ad.fem.3, cont.ad.m.3, cont.ad.fem.3)
# at time 3, there were equal proportions males and females in LD or control
# since these two time points were only one week apart, the inconsistencies viewed between T3 and T4 in the sex ratios suggest there is no effect.

#
#######  04 - SALMON QPCR #############
# rm(list=ls())
setwd("~/Documents/koop/immunostim/03_analysis/05_salmon_qPCR_analysis")

# read in data
salmon.gx <- read.csv(file = "salmon_qpcr-oct23-14-longform.csv")
names(salmon.gx)
levels(salmon.gx$gene) #genes are: il-1; il-8; mmp9; tlr9
levels(salmon.gx$tissue) #tissues are: hk; skin; spleen
str(salmon.gx) #6th column is numerical, and contains the expression values
min(salmon.gx$expr) # the expression data is in linear form currently (i.e. no negative values)
### 2^-ddCt?

#log2 transform expression data
salmon.gx[,6] <- log2(salmon.gx[,6])
str(salmon.gx) #confirm log2 transform of fold change values...

# obtain the correct order of the treatment
levels(salmon.gx$treatment)
salmon.gx$treatment <- factor(salmon.gx$treatment, 
                                levels = c("control","LD","HD"))
levels(salmon.gx$treatment)

#### Analysis and Generate Figures for Manuscript: 
## Plot the four genes ignoring time for visualization (i.e. there is no significant interaction effect),
    #but keep time in models (note: see full boxplots in supplemental)

# Generate salmon qPCR figure:
par(mfrow=c(2,2))
gene <- levels(salmon.gx$gene)

for (i in 1:4) {
  boxplot(salmon.gx$expr[salmon.gx$gene == gene[i] & salmon.gx$tissue == "skin"] ~ 
            salmon.gx$treatment[salmon.gx$gene == gene[i] & salmon.gx$tissue == "skin"]
          ,
          col = rep(c("white","lightgrey","darkgrey")),
          las = 1, ylab = paste("log2(",gene[i], ") in skin", sep = ""), xaxt = "n", ylim =c(-5,4.5))
  axis(1, at = seq(1:3), labels = rep(c("CTR", "LD", "HD")))
  if (gene[i] == "il-1") text(x = seq(1:3), y = 2.5 ,labels = c("a","b","a"))
  else text(x = seq(1:3), y = 2.7 ,labels = "a")  
}

# save 4x5 inch in landscape = figure qPCR salmon

# run models
mods <- as.list(NULL[[1]])
for (i in 1:4) {
  cat(gene[i])
  mods[[i]] <- aov(salmon.gx$expr[salmon.gx$gene == gene[i] & salmon.gx$tissue == "skin"] ~ 
        salmon.gx$treatment[salmon.gx$gene == gene[i] & salmon.gx$tissue == "skin"]
      * salmon.gx$time[salmon.gx$gene == gene[i] & salmon.gx$tissue == "skin"])  
}

# how to access the data:
gene
gene[1]
summary(mods[[1]]) #provides the summary of the aov for gene[i]
TukeyHSD(mods[[1]])



##### TO BE PUT INTO LOOP ####

# Supplemental plots
### Within gene, cross tissue for Supplemental file: ###
#### il-1
par(mfrow=c(3,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
boxplot(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"] ~ 
        salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"]
        * salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "skin log2(il1b)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "spleen"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "spleen"]
        * salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "spleen"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "spleen log2(il1b)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "hk"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "hk"]
        * salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "hk"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "hk log2(il1b)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)


####il-8
par(mfrow=c(3,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
boxplot(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"]
        * salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "skin log2(il-8)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"]
        * salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "spleen log2(il-8)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"]
        * salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "hk log2(il-8)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

####mmp9
par(mfrow=c(3,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
boxplot(salmon.gx$expr[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "skin"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "skin"]
        * salmon.gx$time[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "skin"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "skin log2(mmp9)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "spleen"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "spleen"]
        * salmon.gx$time[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "spleen"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "spleen log2(mmp9)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "hk"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "hk"]
        * salmon.gx$time[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "hk"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "hk log2(mmp9)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

####tlr9
par(mfrow=c(3,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
boxplot(salmon.gx$expr[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "skin"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "skin"]
        * salmon.gx$time[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "skin"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "skin log2(tlr9)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "spleen"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "spleen"]
        * salmon.gx$time[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "spleen"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "spleen log2(tlr9)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

boxplot(salmon.gx$expr[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "hk"] ~ 
          salmon.gx$treatment[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "hk"]
        * salmon.gx$time[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "hk"],
        col = rep(c("white","lightgrey","darkgrey")),
        las = 1, ylab = "hk log2(tlr9)", xaxt = "n")
axis(1, at = seq(1:12), labels = rep(c("CTR", "LD", "HD"), times = 4))
abline(v = c(3.5))
abline(v = c(6.5,9.5),lty=2)

# END PLOTS; START MODELS:


# Models #
## il-1 ## skin
mod.il1b.s <- aov(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"] ~ 
                  as.factor(salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"])
                * salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"])
summary(mod.il1b.s)
anova(mod.il1b.s)
TukeyHSD(mod.il1b.s)
#main treatment = 0.0018; no significant effect of time, and no significant interaction (p > 0.3)

#il-1 hk
mod.il1b.hk.1 <- aov(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "hk"] ~ 
                       as.factor(salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "hk"])
                     * salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "hk"])
summary(mod.il1b.hk.1)
#effect of time (p = 0.00027), no effect of treatment, no interaction

#il-1 spleen
mod.il1b.sp.1 <- aov(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "spleen"] ~ 
                       salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "spleen"]
                     * salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "spleen"])
summary(mod.il1b.sp.1)
TukeyHSD(mod.il1b.sp.1)
#marginal significant time effect (p = 0.03) and treatment (p = 0.0407) almost interact: 0.1385
# General conclusion: local production of il-1b, some indication of up-reg at spleen (no Tukey signif), not at HK


## il-8 ## skin
mod.il8.s <- aov(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"] ~ 
                    salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"]
                  * salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"])
summary(mod.il8.s)
TukeyHSD(mod.il8.s)
# Effect of time (decr slightly over time); no effect of treatment (p = 0.074)

#il-8 hk
mod.il8.hk.1 <- aov(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"] ~ 
                       salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"]
                     * salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"])
summary(mod.il8.hk.1)
#very strong effect of time (decr. over time), no effect of treatment. #marginal interact p =0.08

#il-8 spleen
mod.il8.sp.1 <- aov(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"] ~ 
                       salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"]
                     * salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"])
summary(mod.il8.sp.1)
#significant effect of time p = 0.0014; no other effect
# General conclusion: no DE for il-8 over treatments; but gradually decreases over time in all tissues


## mmp9 ## skin
mod.mmp9.s <- aov(salmon.gx$expr[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "skin"] ~ 
                    salmon.gx$time[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "skin"]
                  * salmon.gx$treatment[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "skin"])
summary(mod.mmp9.s)
# no significant effects (time p = 0.07)

#mmp9 hk
mod.mmp9.hk <- aov(salmon.gx$expr[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "hk"] ~ 
                       salmon.gx$time[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "hk"]
                     * salmon.gx$treatment[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "hk"])
summary(mod.mmp9.hk)
# Strong down-regulation at t2 in all groups (signif effect of time)

#mmp9 spleen
mod.mmp9.sp <- aov(salmon.gx$expr[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "spleen"] ~ 
                       salmon.gx$time[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "spleen"]
                     * salmon.gx$treatment[salmon.gx$gene == "mmp9" & salmon.gx$tissue == "spleen"])
summary(mod.mmp9.sp)
# No effect
# General conclusion: interesting down-regulation after first infection in hk, other than that no effect

## tlr9 ## skin
mod.tlr9.s <- aov(salmon.gx$expr[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "skin"] ~ 
                    salmon.gx$time[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "skin"]
                  * salmon.gx$treatment[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "skin"])
summary(mod.tlr9.s)
TukeyHSD(mod.tlr9.s)
# marginal effect of treatment ( p = 0.035) and effect size is very minor

#tlr9 hk
mod.tlr9.hk <- aov(salmon.gx$expr[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "hk"] ~ 
                       salmon.gx$time[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "hk"]
                     * salmon.gx$treatment[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "hk"])
summary(mod.tlr9.hk)
# effect of time (decr over time, slight interaction effect) 

#tlr9 spleen
mod.tlr9.sp <- aov(salmon.gx$expr[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "spleen"] ~ 
                       salmon.gx$time[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "spleen"]
                     * salmon.gx$treatment[salmon.gx$gene == "tlr9" & salmon.gx$tissue == "spleen"])
summary(mod.tlr9.sp)
# no effect
# General conclusion: no effect