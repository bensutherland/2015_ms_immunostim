# Salmon qPCR for the Ssal/Lsal immunostim project #

# rm(list=ls())

# Todo: Put suppl. plots into a for-loop

#### 4. Salmon qPCR ####
setwd("~/Documents/koop/immunostim/2017_ms_immunostim/")

#### 4.a. Import data ####
# read in data (note: the data was received as linearized ddCT data)
salmon.gx <- read.csv(file = "02_raw_data/salmon_qpcr-oct23-14-longform.csv")
names(salmon.gx)
levels(salmon.gx$gene) #genes are: il-1; il-8; mmp9; tlr9
levels(salmon.gx$tissue) #tissues are: hk; skin; spleen
levels(salmon.gx$time) # note that there is no T0 in the normalized expr data
str(salmon.gx) # expression data in 6th column 'expr', and is is numerical
min(salmon.gx$expr) # the expression data is in linear form currently (i.e. no negative values)


#### 4.b. Transform data ####
# log2 transformation
salmon.gx[,6] <- log2(salmon.gx[,6])
str(salmon.gx) #confirm log2 transform of fold change values...

# Set the correct order of treatments for plotting
levels(salmon.gx$treatment)
salmon.gx$treatment <- factor(salmon.gx$treatment, 
                              levels = c("control","LD","HD"))
levels(salmon.gx$treatment)


#### 4.c. Analyze data ####
## Plot the four genes ignoring time for visualization (i.e. there is no significant interaction effect),
# but keep time in models (note: see full boxplots in supplemental)

# Plot with time averaged (for manuscript):
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

# Run models
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
#note that ld-control in tukey is diff = 1.16 and p = 0.0056
2^1.16 # linear = 2.2-fold overexpressed when not considering time (but see suppl. fig.)


# double check model for il1b
par(mfrow=c(1,1))
boxplot(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"] ~ 
                  salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"]
                * salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"])  

mod.test <- aov(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"] ~ 
                  salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"]
                * salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"])  
summary(mod.test)


#### 4.d. Create supplemental plots ####
# Within gene, cross tissue for Supplemental file:
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


#### 4.e. Create supplemental models
## il-1 ## 
# il-1 skin
mod.il1b.s <- aov(salmon.gx$expr[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"] ~ 
                    as.factor(salmon.gx$time[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"])
                  * salmon.gx$treatment[salmon.gx$gene == "il-1" & salmon.gx$tissue == "skin"])
summary(mod.il1b.s)
anova(mod.il1b.s)
TukeyHSD(mod.il1b.s)
#main treatment = 0.0018; no significant effect of time, and no significant interaction (p > 0.3)

# il-1 hk
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
# Marginal significant time effect (p = 0.03) and treatment (p = 0.0407) almost interact: 0.1385
# General conclusion: local production of il-1b, some indication of up-reg from feed at spleen (no Tukey signif), not at HK


## il-8 ## 
#il-8 skin
mod.il8.s <- aov(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"] ~ 
                   salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"]
                 * salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "skin"])
summary(mod.il8.s)
TukeyHSD(mod.il8.s)
# Effect of time (decr slightly over time); no effect of treatment (p = 0.074)

# il-8 hk
mod.il8.hk.1 <- aov(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"] ~ 
                      salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"]
                    * salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "hk"])
summary(mod.il8.hk.1)
#very strong effect of time (decr. over time), no effect of treatment. #marginal interact p =0.08

# il-8 spleen
mod.il8.sp.1 <- aov(salmon.gx$expr[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"] ~ 
                      salmon.gx$time[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"]
                    * salmon.gx$treatment[salmon.gx$gene == "il-8" & salmon.gx$tissue == "spleen"])
summary(mod.il8.sp.1)
#significant effect of time p = 0.0014; no other effect
# General conclusion: no DE for il-8 over treatments; but gradually decreases over time in all tissues


## mmp9 ## 
#mmp9 skin
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

## tlr9 ## 
#tlr9 skin
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