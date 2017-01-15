# Transcriptome analysis for the Ssal/Lsal immunostim project

# rm(list=ls())

## Install packages
# source("http://www.bioconductor.org/biocLite.R")
# biocLite()
# biocLite("limma")
library(limma)


#### 1. Lice transcriptome analysis ####
setwd("~/Documents/koop/immunostim/2017_ms_immunostim")

#### 1.a. Input data and quality control ####
# Load targets file (i.e. interpretation file)
targets <- readTargets("00_archive/Targets_2color_immunostim.csv", row.names="Name", sep = ",")

# Load annotation file
annot <- read.table("00_archive/koop_Lsal38K_2010-02-25b-short.txt",
                    header=TRUE, sep = "\t") #note sep = "\t" is important when empty cells
annot.df <- as.data.frame(annot)
dim(annot.df) #38132 rows; 5 columns - note: this should match the length of MA.flt1 given above
str(annot.df)

# Load expression data (block and channel sep.)
RG <- read.maimages(path = "02_raw_data/",
                    files = targets[,c("FileNameCy3","FileNameCy5")], 
                    source="imagene",
                    columns = list(f="Signal Median", b="Background Median"))
dim(RG$R) # number of probes and number samples
targets # see details on the samples

# View density plot of unprocessed data
plotDensities(RG, log=F) # note that there is a bump at the saturation point
plotDensities(RG, log=T) # this makes it easier to view the distribution, will produce NAs

# Set probe control type
spottypes <- readSpotTypes(file = "00_archive/SpotTypes.txt")
RG$genes$Status <- controlStatus(types = spottypes, genes = RG)

# Quality Control plots ##
plotMA3by2(RG, path = "03_analysis/") #MA plots per array saved to working directory
imageplot3by2(RG, path = "03_analysis/") #imageplots per array saved to working directory

# Boxplots non-normalized foreground/background per channel
par(mfrow=c(2,2), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

units <- c("R", "G", "Rb", "Gb")
for(i in 1:4) {
  boxplot(log2(RG[[units[i]]]), 
          xaxt = "n", 
          ylab = "log2(fluor. units)", xlab = "samples", main = units[i])
}
print(dimnames(RG[["R"]])[[2]]) # this is the order of the plot


#### 1.b. Quality filtering and normalization ####
par(mfrow=c(3,1), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

## Set low expressed genes statement in $RG$genes$isexpr column
RG$genes$isexpr.R <- rowSums(RG$R >= 500) >= 3 #here, value >= 500 in at least 3 arrays
                      ##isexpr > 500 in at least 3 arrays, lowest sample size = 4 for ld.
table(RG$genes$isexpr.R)
RG$genes$isexpr.G <- rowSums(RG$G >= 500) >= 3 ##isexpr > 500 in at least 3 arrays
table(RG$genes$isexpr.G)

## Identify saturated probes
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

## Attach annotation to the MA list
MA.flt1$genes$spotID <- annot.df[,1] 
MA.flt1$genes$description <- annot.df[,2] 
MA.flt1$genes$GeneID <- annot.df[,4] 

## Remove low expressed and saturated probes
MA.flt <- MA.flt1[MA.flt1$genes$isexpr.R == "TRUE" & MA.flt1$genes$isexpr.G == "TRUE"
                  & MA.flt1$genes$satspot.R != "TRUE" & MA.flt1$genes$satspot.G != "TRUE",]
dim(MA.flt)

plotDensities(MA.flt) ##densitiy plot of expression values

#### 1.c. Data extraction ####
## Save out quality filtered MAlist
write.csv(cbind(MA.flt$genes,MA.flt$M), "04_output/all_expr_data.csv", row.names=FALSE)

## Select probe expression specifically for correlation with qPCR GOIs 
probes.of.interest <- c("C042R126","C088R114","C263R087")
probes.for.cor <- as.data.frame(cbind(
              MA.flt$M[MA.flt$genes[12] == probes.of.interest[1]]
            , MA.flt$M[MA.flt$genes[12] == probes.of.interest[2]] 
            , MA.flt$M[MA.flt$genes[12] == probes.of.interest[3]]))
rownames(probes.for.cor) <- colnames(MA.flt$M)
colnames(probes.for.cor) <- probes.of.interest
probes.for.cor
str(probes.for.cor)

## Save out probes.for.cor
saveRDS(probes.for.cor, file = "04_output/probes_for_cor.rds")

## Save out background list for GO enrichment testing
write.csv(MA.flt$genes, file = "04_output/background_immunostim.csv")


#### 1.c. Differential expression analysis ####
### Differential Expression Analysis ## 
design <- modelMatrix(targets, ref = "ref")
attributes(design)
design

fit <- lmFit(object = MA.flt, design = design)
cont.matrix <- makeContrasts(ld-con, levels=design) ##makes contrasts
fit2 <- contrasts.fit(fit, cont.matrix) ###computes contrasts from lmod fit (outputs)
fit3 <- eBayes(fit2) ##given related parameter estim and SE, will comute mod t-stat, etc
output <- topTable(fit3, number=45000, adjust="none", p.value=0.05) # number=45000 makes sure all will be output
write.table(output, file="04_output/transcriptome_DE_results.txt", sep="\t")
output.df <- data.frame(output)