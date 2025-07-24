library(reshape2)
library(lmerTest)
library(dplyr)
library(plyr)
library(survival)
library(vegan)
library(ape)
library(emmeans)

##### 1 - Discovery in the MetaSalt study ===========================================================================

### 1) Gut-microbial diversity --------------------------------------------------------

## α-diversity --------------
species_alpha <- alpha_diversity(species.abundance)
species_alpha$labid <- rownames(species_alpha)
species_alpha <- merge(species_alpha, alldata, by="labid")

# Difference across study phases 
dat1 <- subset(species_alpha, phaseid!=3)
fit <- lmer(data=dat1, Simpson~phaseid+(1|famid)+(1|labid2))
p.12 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(Simpson, phaseid!=2)
fit <- lmer(data=dat1, Simpson~phaseid+(1|famid)+(1|labid2))
p.13 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(Simpson, phaseid!=1)
fit <- lmer(data=dat1, Simpson~phaseid+(1|famid)+(1|labid2))
p.23 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
alpha_phase.p <- data.frame(p.12, p.13, p.23)

# Difference across SSBP group in each study phase
dat1 <- subset(species_alpha, phaseid==1)
fit <- lmer(data=dat1, Simpson~ss3+age+gender+bmi+fcc+smoking+htn+chol+(1|famid))
p.1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(species_alpha, phaseid==2)
fit <- lmer(data=dat1, Simpson~ss3+age+gender+bmi+fcc+smoking+htn+chol+(1|famid))
p.2 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(species_alpha, phaseid==3)
fit <- lmer(data=dat1, Simpson~ss3+age+gender+bmi+fcc+smoking+htn+chol+(1|famid))
p.3 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
alpha_phase_ss.p <- data.frame(p.1, p.2, p.3)

## β-diversity ----------------
species_beta <- vegdist(species.abundance, method = "bray")
species_beta <- as.data.frame(as.matrix(species_beta))
species_beta.pcoa <- pcoa(species_beta)

temp.1 <- species_beta
temp.1 <- temp.1[sort(rownames(temp.1)), sort(colnames(temp.1))]
temp.1 <- as.dist(temp.1)
temp.2 <- alldata[sort(rownames(alldata)), -1]
set.seed(20230601)
species_beta.p <- adonis2(temp.1~phaseid, method="bray", data = temp.2, permutations = 999, parallel = 16)

## Median Bray-Curtis dissimilarity -------------
temp.1 <- species_beta[paste0("B", labid2), paste0("B", labid2)]
temp.1 <- betadis(temp.1)
temp.2 <- species_beta[paste0("L", labid2), paste0("L", labid2)]
temp.2 <- betadis(temp.2)
temp.3 <- species_beta[paste0("H", labid2), paste0("H", labid2)]
temp.3 <- betadis(temp.3)
temp.4 <- rbind(temp.1, temp.2, temp.3)
median_bray_dis <- merge(alldata, temp.4, by="labid")

# Difference across study phase
dat1 <- subset(median_bray_dis, phaseid!=3)
fit <- lmer(data=dat1, log(betadis)~phaseid+(1|famid)+(1|labid2))
p.12 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(median_bray_dis, phaseid!=2)
fit <- lmer(data=dat1, log(betadis)~phaseid+(1|famid)+(1|labid2))
p.13 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(median_bray_dis, phaseid!=1)
fit <- lmer(data=dat1, log(betadis)~phaseid+(1|famid)+(1|labid2))
p.23 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
median_phase.p <- data.frame(p.12, p.13, p.23)

# Difference across SSBP group in each study phase
dat1 <- subset(median_bray_dis, phaseid==1)
fit <- lmer(data=dat1, log(betadis)~ss3+age+gender+bmi+fcc+smoking+htn+chol+(1|famid))
p.1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(median_bray_dis, phaseid==2)
fit <- lmer(data=dat1, log(betadis)~ss3+age+gender+bmi+fcc+smoking+htn+chol+(1|famid))
p.2 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(median_bray_dis, phaseid==3)
fit <- lmer(data=dat1, log(betadis)~ss3+age+gender+bmi+fcc+smoking+htn+chol+(1|famid))
p.3 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
median_phase_ss.p <- data.frame(p.1, p.2, p.3)

### 2) - Identification of salt-related biomarkers ------------------------------------

# salt-related gut-microbial species ----------------
temp <- melt(alldata, c("labid2","famid","phaseid"), speid)
metasalt.sr.species <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.species <- subset(metasalt.sr.species, p.h<0.05/531)

# salt-related gut-microbial functions ----------------
temp <- melt(alldata, c("labid2","famid","phaseid"), keggid)
metasalt.sr.kegg <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.kegg <- subset(metasalt.sr.kegg, p.h<0.05/190)

temp <- melt(alldata, c("labid2","famid","phaseid"), cazyid)
metasalt.sr.cazy <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.cazy <- subset(metasalt.sr.cazy, p.h<0.05/133)

# salt-related metabolites ---------------------------
temp <- melt(alldata, c("labid2","famid","phaseid"), metid)
metasalt.sr.metabol <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.metabol <- subset(metasalt.sr.metabol, p.h<0.05/221)

### 2) - Identification of SSBP-related biomarkers ------------------------------------

# SSBP-related gut-microbial species ------------------------
temp <- melt(chgdata, c("labid2","famid","ss","age","gender","bmi","fcc","htn","smoking","chol"), metasalt.sr.species$speid)
metasalt.ss.species <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
  se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
  se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
  p2 = data.frame(summary(fit)$coefficients)$Pr...t..[3]
  
  output = data.frame(coef.ssbp=coef, se.ssbp=se, p.ssbp=p,
                      coef.ssbp.1=coef1, se.ssbp.1=se1, p.ssbp.1=p1,
                      coef.ssbp.2=coef2, se.ssbp.2=se2, p.ssbp.2=p2)
  print(output)
})

metasalt.ss.species <- subset(metasalt.ss.species, p.ssbp<0.05)
species.lsmeans <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp.1, value~as.factor(ss)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="y") %>% data.frame()
})

# SSBP-related gut-microbial function ----------------------
temp <- melt(chgdata, c("labid2","famid","ss","age","gender","bmi","fcc","htn","smoking","chol"), metasalt.sr.kegg$keggid)
metasalt.ss.species <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
  se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
  se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
  p2 = data.frame(summary(fit)$coefficients)$Pr...t..[3]
  
  output = data.frame(coef.ssbp=coef, se.ssbp=se, p.ssbp=p,
                      coef.ssbp.1=coef1, se.ssbp.1=se1, p.ssbp.1=p1,
                      coef.ssbp.2=coef2, se.ssbp.2=se2, p.ssbp.2=p2)
  print(output)
})

metasalt.ss.kegg <- subset(metasalt.ss.kegg, p.ssbp<0.05)
kegg.lsmeans <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp.1, value~as.factor(ss)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="y") %>% data.frame()
})

temp <- melt(chgdata, c("labid2","famid","ss","age","gender","bmi","fcc","htn","smoking","chol"), metasalt.sr.cazy$cazyid)
metasalt.ss.species <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
  se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
  se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
  p2 = data.frame(summary(fit)$coefficients)$Pr...t..[3]
  
  output = data.frame(coef.ssbp=coef, se.ssbp=se, p.ssbp=p,
                      coef.ssbp.1=coef1, se.ssbp.1=se1, p.ssbp.1=p1,
                      coef.ssbp.2=coef2, se.ssbp.2=se2, p.ssbp.2=p2)
  print(output)
})

metasalt.ss.cazy <- subset(metasalt.ss.cazy, p.ssbp<0.05)
cazy.lsmeans <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp.1, value~as.factor(ss)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="y") %>% data.frame()
})

# SSBP-related metabolites ------------------------------------
temp <- melt(chgdata, c("labid2","famid","ss","age","gender","bmi","fcc","htn","smoking","chol"), metasalt.sr.metabol$metid)
metasalt.ss.metabol <-  ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
  se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
  se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
  p2 = data.frame(summary(fit)$coefficients)$Pr...t..[3]
  
  output = data.frame(coef.ssbp=coef, se.ssbp=se, p.ssbp=p,
                      coef.ssbp.1=coef1, se.ssbp.1=se1, p.ssbp.1=p1,
                      coef.ssbp.2=coef2, se.ssbp.2=se2, p.ssbp.2=p2)
  print(output)
})

metasalt.ss.metabol <- subset(metasalt.ss.metabol, p.ssbp<0.05)
metabol.lsmeans <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp.1, value~as.factor(ss)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="y") %>% data.frame()
})

# Associations between SSBP-related species and metabolites
temp <- melt(chgdata, c("labid2","famid","ss","age","gender","bmi","fcc","htn","smoking","chol", metasalt.ss.metabol$metid), metasalt.ss.species$speid)
names(temp)[which(names(temp) %in% c("variable", "value"))] = c("speid", "spe.value")
temp <- melt(temp, c("labid2","famid","ss","age","gender","bmi","fcc","htn","smoking","chol", "speid", "spe.value"), metasalt.ss.metabol$metid)
names(temp)[which(names(temp) %in% c("variable", "value"))] = c("metid", "met.value")
assoc.spe.met <- ddply(temp, c("speid","metid"), function(temp){
  fit <- lmer(data=temp, met.chg~spe.chg+age+gender+bmi+fcc+chol+smoking+htn+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  res <- data.frame(coef, se, p)
  print(res)
})
assoc.spe.met <- subset(assoc.spe.met, p<0.05)


##### 2 - Replication in the GenSalt study =============================================================================

## 1) Replication of salt-related metabolites -------------------------------------
temp <- melt(alldata.gst, c("labid2","phaseid"), metasalt.sr.metabol$metid)
gensalt.sr.metabol <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  res <- data.frame(coef.h=coef, se.h=se, p.h=p)
  print(res)
})
gensalt.sr.metabol <- subset(gensalt.sr.metabol, p.h<0.05)

## 2) Replication of SSBP-related metabolites --------------------------------------
temp <- melt(chgdata.gst, c("labid2","famid","ss","bmi","fcc","smoking","chol"), metasalt.ss.metabol$metid)
gensalt.ss.metabol <-  ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss+bmi+fcc+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss)+bmi+fcc+smoking+chol+(1|famid))
  coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
  se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p1 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
  se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
  p2 = data.frame(summary(fit)$coefficients)$Pr...t..[3]
  
  output = data.frame(coef.ssbp=coef, se.ssbp=se, p.ssbp=p,
                      coef.ssbp.1=coef1, se.ssbp.1=se1, p.ssbp.1=p1,
                      coef.ssbp.2=coef2, se.ssbp.2=se2, p.ssbp.2=p2)
  print(output)
})

gensalt.ss.metabol <- subset(gensalt.ss.metabol, p.ssbp<0.05)
metabol.lsmeans.gst <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp.1, value~as.factor(ss)+fcc+bmi+smoking+chol+(1|famid))
  ls <- emmeans(fit, specs="y") %>% data.frame()
})

##### 3 - Validation in Rats ===========================================================================================

### 1) Wistar Rats --------------------------------------------------------------------------------

# Differences in BP levels and isovalerylcarnitine
aov(data=alldata.wistar, SBP~group*time+Error(ratid/time))
aov(data=alldata.wistar, DBP~group*time+Error(ratid/time))
t.test(data=alldata.wistar, isoval~group, var.equal=T)

### 2) Dahl-SS rats ------------------------------------------------------------------------------

# Differences in BP levels
aov(data=alldata.dahl, SBP~group*Time+Error(num/Time))
aov(data=subset(alldata.dahl, group %in% c("HSD","HSD_isoval")),  SBP~group*Time+Error(num/Time))
t.test(data=subset(alldata.dahl, Time=="28"), SBP~group, var.equal=T)

aov(data=alldata.dahl, DBP~group*Time+Error(num/Time))
aov(data=subset(alldata.dahl, group %in% c("HSD","HSD_isoval")),  DBP~group*Time+Error(num/Time))
t.test(data=subset(alldata.dahl, Time=="28"), DBP~group, var.equal=T)

# Differences in isovalerylcarnitine
t.test(data=subset(temp, Group %in% c("Control","HSD")), iso~Group, var.equal=T)
t.test(data=subset(temp, Group %in% c("Control","HSD_isoval")), iso~Group, var.equal=T)
t.test(data=subset(temp, Group %in% c("HSD","HSD_isoval")), iso~Group, var.equal=T)

# Differences in the relaxation of mesenteric arteries
aov(data=alldata.dahl, acetyl~group*con+Error(num/con))
aov(data=subset(alldata.dahl, group %in% c("HSD","HSD_isoval")), acetyl~group*con+Error(num/con))
aov(data=alldata.dahl, nitrop~group*con+Error(num/con))
aov(data=subset(alldata.dahl, group %in% c("HSD","HSD_isoval")), nitrop~group*con+Error(num/con))

### 3) SD rats -----------------------------------------------------------------------------------

# Differences in ASV richness
t.test(data=subset(alldata.sd, group %in% c("Control","ACT")), ASV_num~group, var.equal=T)
t.test(data=subset(alldata.sd, group %in% c("Control","ACT to FMT")), ASV_num~group, var.equal=T)
t.test(data=subset(alldata.sd, group %in% c("ACT","ACT to FMT")), ASV_num~group, var.equal=T)

# Differences in isovalerylcarnitine
t.test(data=subset(alldata.sd, group %in% c("Control","ACT")), isoval~group, var.equal=T)
t.test(data=subset(alldata.sd, group %in% c("Control","ACT to FMT")), isoval~group, var.equal=T)
t.test(data=subset(alldata.sd, group %in% c("ACT","ACT to FMT")), isoval~group, var.equal=T)

##### 4 - Validation in the cohort study ===============================================================================

### 1) Prevalent hypertension -------------------------------------------------
fit <- glm(data=cohort.hyp1, hyp~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine, family=binomial(link="logit"))
coef = data.frame(summary(fit)$coefficients)$Estimate[2]
se = data.frame(summary(fit)$coefficients)$Std..Error[2]
p = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- glm(data=cohort.hyp1, hyp~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine, family=binomial(link="logit"))
coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[2]
coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[3]

fit <- glm(data=cohort.hyp1, hyp~isoval.group+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine, family=binomial(link="logit"))
ptrend <- summary(fit)$coefficients[2,4]

AC.hyp1 = data.frame(outcome="Prevalent hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 2) Ideal BP status - incident prehypertension/hypertension --------------
fit <- coxph(data=cohort.hyp2, Surv(time, prehyp_hyp)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$Estimate[1]
se = data.frame(summary(fit)$coefficients)$Std..Error[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort.hyp2, Surv(time, prehyp_hyp)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$Estimate[1]
se1 = data.frame(summary(fit)$coefficients)$Std..Error[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$Estimate[2]
se2 = data.frame(summary(fit)$coefficients)$Std..Error[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort.hyp2, Surv(time, prehyp_hyp)~isoval.group+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp2 = data.frame(outcome="Ideal BP status - Prehypertension/hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 3) Prehypertension - incident hypertension ------------------------------
fit <- coxph(data=cohort.hyp3, Surv(time, inchyp)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$Estimate[1]
se = data.frame(summary(fit)$coefficients)$Std..Error[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort.hyp3, Surv(time, inchyp)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$Estimate[1]
se1 = data.frame(summary(fit)$coefficients)$Std..Error[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$Estimate[2]
se2 = data.frame(summary(fit)$coefficients)$Std..Error[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort.hyp3, Surv(time, inchyp)~isoval.group+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp3 = data.frame(outcome="Prehypertension - Incident hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 4) Normotension - incident hypertension ---------------------------------
fit <- coxph(data=cohort.hyp4, Surv(time, inchyp2)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$Estimate[1]
se = data.frame(summary(fit)$coefficients)$Std..Error[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort.hyp4, Surv(time, inchyp2)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$Estimate[1]
se1 = data.frame(summary(fit)$coefficients)$Std..Error[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$Estimate[2]
se2 = data.frame(summary(fit)$coefficients)$Std..Error[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort.hyp4, Surv(time, inchyp2)~isoval.group+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp4 = data.frame(outcome="Normotension - Incident hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 5) BP status increases >= 1 step ----------------------------------------
fit <- coxph(data=cohort.hyp5, Surv(time, probp)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$Estimate[1]
se = data.frame(summary(fit)$coefficients)$Std..Error[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort.hyp5, Surv(time, probp)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$Estimate[1]
se1 = data.frame(summary(fit)$coefficients)$Std..Error[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$Estimate[2]
se2 = data.frame(summary(fit)$coefficients)$Std..Error[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort.hyp5, Surv(time, probp)~isoval.group+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp5 = data.frame(outcome="BP status increases >= 1 step", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

AC.hyp <- rbind(AC.hyp1, AC.hyp2, AC.hyp3, AC.hyp4, AC.hyp5)


