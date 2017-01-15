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
