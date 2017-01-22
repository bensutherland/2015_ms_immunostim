# Lice counts and fish physiology for the Ssal/Lsal immunostim project

# rm(list=ls())

#### 3. Lice counts and fish physiology ####
setwd("~/Documents/koop/immunostim/2017_ms_immunostim/")

#### 3.a. Import physiological and infection data ####
immuno.data <- read.csv(file = "02_raw_data/immunostim_data-jan21-15.csv")
dim(immuno.data) # 313 rows, 15 columns
names(immuno.data)
immuno.data[1:10,]
levels(as.factor(immuno.data$time)) # Note that data collection time 3.1 and 3.2 were merged

# Set up data and remove the tank control (UNKN)
immuno.data <- subset(immuno.data, immuno.data$treatment != "UNKN")
immuno.data <- droplevels(immuno.data)
levels(immuno.data$treatment) # Confirm UNKN has now been removed
immuno.data$treatment <- factor(immuno.data$treatment, levels = c("cont","LD","HD")) # Relevel to put the right order (control, lowdose, highdose)
                                

#### 3.b. Pre-experiment analysis ####
# Eval weight at T0 to confirm no diff at start
# Note: not enough reps for indep. tank eval, so do by group instead of by tank
# Plot with individual records
boxplot(immuno.data$weight.g[immuno.data$time==0] ~ 
          immuno.data$treatment[immuno.data$time==0], 
        na.rm =T, col = c("white","lightgrey","darkgrey"),
        xlab = "Tank ID", ylab = "Weight (g)")
points(immuno.data$weight.g[immuno.data$time==0] ~ immuno.data$treatment[immuno.data$time==0]) # n = 6 (total = 18 samples)

# model
pre.inf.mod = aov(immuno.data$weight.g[immuno.data$time==0] ~ immuno.data$treatment[immuno.data$time==0])
summary(pre.inf.mod) # no significant difference prior to infection (ANOVA p > 0.3)


#### 3.b. Weight after start feed (not incl. T0) ####
# Plot
boxplot(immuno.data$weight.g[immuno.data$time!=0] 
        ~ immuno.data$treatment[immuno.data$time!=0] * as.factor(immuno.data$time[immuno.data$time!=0])
        , col = c("white","lightgrey","darkgrey"))

# Model
weight.mod <- aov(immuno.data$weight.g[immuno.data$time!=0]
                  ~ immuno.data$treatment[immuno.data$time!=0] * as.factor(immuno.data$time[immuno.data$time!=0]))
summary(weight.mod) # Sig. fx of treatment (p = 0.02) and of time (p < 0.00001)
TukeyHSD(weight.mod)

# Note: weight data missing in data collection for indiv. 13-22 from T4 (tanks 6, 7, 8, 9, 10) (noted in manuscript)
summary(is.na(immuno.data$weight.g[immuno.data$time==4])) # due to this, variable will be treated as LPF, not using weight

# Eval w/o T0 and T4 (due to missing data)
boxplot(immuno.data$weight.g[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] ~ immuno.data$treatment[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] * 
          as.factor(immuno.data$time[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3])
        , col = c("white","lightgrey","darkgrey"))
weight.mod.noT0orT4 <- aov(immuno.data$weight.g[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3] 
                          ~ immuno.data$treatment[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3]
                          * as.factor(immuno.data$time[immuno.data$time==1 | immuno.data$time == 2 | immuno.data$time == 3]))
summary(weight.mod.noT0orT4) # marginal effect of treatment ( p = 0.076)
TukeyHSD(weight.mod.noT0orT4) # no sig pairwise b/w HD and contr at T1, T2, or T3 # but marginal overall HD-cont paired (p = 0.06)


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
TukeyHSD(lpf.mod) # pair-wise comparisons to see within a time point which values are different
# LD3 v cont3 (p = 0.027)
# LD4 v cont4 (p = 0.016)

# Evaluate fold diff bw cont and ld
mean(immuno.data$per.fish.total[immuno.data$time == 3 & immuno.data$treatment ==  "cont"])
mean(immuno.data$per.fish.total[immuno.data$time == 3 & immuno.data$treatment ==  "LD"])
1/(1.4/2.9) # 2.1-fold lower at time 3
mean(immuno.data$per.fish.total[immuno.data$time == 4 & immuno.data$treatment ==  "cont"])
mean(immuno.data$per.fish.total[immuno.data$time == 4 & immuno.data$treatment ==  "LD"])
1/(1.56/3.48) # 2.2-fold lower at time 4

## Effect of lice sex bw treatment groups?
names(immuno.data)

# collect proportions
LD.ad.m.4 <- sum(immuno.data$body.ad.m[immuno.data$time==4 & immuno.data$treatment=="LD"], na.rm =T)
LD.ad.fem.4 <-  sum(immuno.data$body.ad.fem[immuno.data$time==4 & immuno.data$treatment=="LD"], na.rm =T)
cont.ad.m.4 <- sum(immuno.data$body.ad.m[immuno.data$time==4 & immuno.data$treatment=="cont"], na.rm =T)
cont.ad.fem.4 <- sum(immuno.data$body.ad.fem[immuno.data$time==4 & immuno.data$treatment=="cont"], na.rm =T)
c(LD.ad.m.4, LD.ad.fem.4, cont.ad.m.4, cont.ad.fem.4)
# At time 4, more males than females in the LD, and opposite in the control
# Numbers are pretty low to determine this.

LD.ad.m.3 <- sum(immuno.data$body.ad.m[immuno.data$time==3 & immuno.data$treatment=="LD"], na.rm =T)
LD.ad.fem.3 <- sum(immuno.data$body.ad.fem[immuno.data$time==3 & immuno.data$treatment=="LD"], na.rm =T)
cont.ad.m.3 <- sum(immuno.data$body.ad.m[immuno.data$time==3 & immuno.data$treatment=="cont"], na.rm =T)
cont.ad.fem.3 <- sum(immuno.data$body.ad.fem[immuno.data$time==3 & immuno.data$treatment=="cont"], na.rm =T)
c(LD.ad.m.3, LD.ad.fem.3, cont.ad.m.3, cont.ad.fem.3)
# At time 3, there were equal proportions males and females in LD or control
# Since these two time points were only one week apart, the inconsistencies viewed between T3 and T4 in the sex ratios suggest there is no effect.
