# Lice qPCR analysis for Ssal/Lsal immunostim project

# rm(list=ls())

## Install packages
# source("http://www.bioconductor.org/biocLite.R")
# biocLite()
biocLite("ReadqPCR")
biocLite("NormqPCR") # do not load from source requiring compilation.
biocLite("limma")
library(ReadqPCR)
library(NormqPCR)
library(limma) # needed for reading targets
# libraries for manipulating strings
library(stringi)
library(tidyr)


#### 2. Lice qPCR analysis ####
# set working directory to the git repo main directory
setwd("~/Documents/koop/immunostim/2017_ms_immunostim")

#### 2.a. Input data and quality control ####
# Load qPCR data
immuno.data <- read.qPCR("02_raw_data/new-together_formatted-notdouble.txt")

# Load microarray data
probes.for.cor <- readRDS(file = "04_output/probes_for_cor.rds")

# Load interpretation file
targets <- readTargets("00_archive/Targets_2color_immunostim.csv", row.names="Name", sep = ",")


# Inspection of batch object
dim(immuno.data) # 30 features, 11 samples
rownames(exprs(immuno.data)) # rows are genes (incl reps)
colnames(exprs(immuno.data)) # cols are samples
exprs(immuno.data) #lets you observe the raw Cts - note: there are two NA values for tech reps
which(is.na(exprs(immuno.data))) # currently two NAs

# Combine technical replicates (and remove NAs)
combinedTechReps <- combineTechReps(immuno.data)
combinedTechReps #this combined gene names
exprs(combinedTechReps)  # gives Cts
featureNames(combinedTechReps) # gives gene names
which(is.na(exprs(combinedTechReps))) # in averaging, NAs are ignored; thus NAs removed

# Remove genes that are not to be used in the analysis
row.names(combinedTechReps) # gives gene names
combinedTechReps <- combinedTechReps[-1,] # remove abcc1
row.names(combinedTechReps)

# Calculate which are the best normalizers using geNORM
normalizers <- selectHKs(combinedTechReps, method = "geNorm",
                         Symbols = featureNames(combinedTechReps),
                         minNrHK = 2, log = T) #selectHKs is run with log = TRUE since we are using Cq
ranks <- data.frame(c(1, 1:13), normalizers$ranking)

# dCt
hkgs <- c("flna", "rps20")
dCt <- deltaCq(qPCRBatch =  combinedTechReps, hkgs = hkgs, calc="geo") # subtracts hkg values
head(exprs(dCt)) # in dCt format (i.e. still a log2 value)
# Note: run statistics on log2 values

# Create a dataframe using the dCt expression data, rows as samples, columns as genes
# Note: need to transform it to be negative as it is still higher = fewer transcripts
results <- as.data.frame(t(exprs(dCt) * -1))
head(results)

#### 2.b. Correlate microarray vs qPCR ####
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
# C088R114 Nuclear pore membrane glycoprotein 210
# C263R087 Histone-lysine N-methyltransferase SETD7

# models first to obtain R-squared
## Still need to bring in the two datasets #TODO#

# match together the two dataframes
results$Name2match <- rownames(results)
head(results)
head(probes.for.cor.trimmed)

collected.df <- merge(x = probes.for.cor.trimmed, y = results, by = "Name2match")
names(collected.df)

sg1 <- NULL; pom210 <- NULL; setd7 <- NULL
gois <- c("sg1","pom210","setd7")
probes <- c("C042R126","C088R114","C263R087")
lm.results <- list()
adj.r.sq.res <- NULL

for(i in 1:3){
  lm.results[[gois[i]]] <- lm(collected.df[,gois[i]] ~ collected.df[,probes[i]])
  adj.r.sq.res[i] <- summary(lm.results[[i]])$adj.r.squared
  names(adj.r.sq.res)[i] <- gois[i]
}

rsq <- adj.r.sq.res

# Plot
par(mfrow=c(1,3), mar = c(4,4,1,2))
xtxt = -5.5
ytxt = 3.7

for(i in 1:3){
  plot(x = collected.df[,probes[i]]
       , y = collected.df[,gois[i]]
       , xlab = paste("log2(", probes[i],")", sep = "")
       , ylab = paste("log2(",gois[i],")", sep = "")
       , las = 1)
  text(x = max(collected.df[,probes[i]])-2, 
       y = min(collected.df[,gois[i]])+0.5
       , labels = paste("R2 =", round(rsq[i], digits = 2)))
}

# save out as 8 x 4 plot


#### 2.c. ddCt calculations ####
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