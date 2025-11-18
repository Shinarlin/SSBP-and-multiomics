library(reshape2)
library(lmerTest)
library(dplyr)
library(plyr)
library(survival)
library(vegan)
library(ape)
library(emmeans)

files <- list.files(pattern="*.rds")
for (i in files){
  assign(gsub(".rds","",i),readRDS(i))
}

##### 1 - Discovery in the MetaSalt study ===========================================================================

### 1) Gut-microbial diversity --------------------------------------------------------

## α-diversity --------------
alpha_diversity <- function(x, tree = NULL){
  observed_species <- estimateR(x)[1,]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
  return(result)
}

species_alpha <- alpha_diversity(species.abundance)
species_alpha$labid <- rownames(species_alpha)
species_alpha <- merge(species_alpha, alldata.mst, by="labid")
species_alpha$Simpson <- as.numeric(species_alpha$Simpson)

# Difference across study phases 
dat1 <- subset(species_alpha, phaseid!=3)
fit <- lmer(data=dat1, Simpson~phaseid+(1|famid)+(1|labid2))
p.12 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(species_alpha, phaseid!=2)
fit <- lmer(data=dat1, Simpson~phaseid+(1|famid)+(1|labid2))
p.13 = data.frame(summary(fit)$coefficients)$Pr...t..[2]
dat1 <- subset(species_alpha, phaseid!=1)
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
temp.2 <- alldata.mst[sort(rownames(alldata.mst)), -1]
set.seed(20230601)
species_beta.p <- adonis2(temp.1~phaseid, method="bray", data = temp.2, permutations = 999, parallel = 16)

## Median Bray-Curtis dissimilarity -------------
betadis <- function(indata=data){
  temp.1 <- indata
  for (i in 1:nrow(temp.1)){
    for (j in 1:ncol(temp.1)){
      if (i==j) {
        temp.1[i,j] = NA
      }
    }
  }
  temp.1$labid <- rownames(temp.1)
  temp.2 <- melt(temp.1, "labid") %>% na.omit()
  temp.3 <- aggregate(data=temp.2, value~labid, "median")
  names(temp.3)[2] = "betadis"
  return(temp.3)
}

labid2 <- unique(alldata.mst$labid2)
temp.1 <- species_beta[paste0("B", labid2), paste0("B", labid2)]
temp.1 <- betadis(temp.1)
temp.2 <- species_beta[paste0("L", labid2), paste0("L", labid2)]
temp.2 <- betadis(temp.2)
temp.3 <- species_beta[paste0("H", labid2), paste0("H", labid2)]
temp.3 <- betadis(temp.3)
temp.4 <- rbind(temp.1, temp.2, temp.3)
median_bray_dis <- merge(alldata.mst, temp.4, by="labid")

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
speid <- names(alldata.mst)[grepl("spe", names(alldata.mst))]
temp <- melt(alldata.mst, c("labid2","famid","phaseid"), speid)
temp <- subset(temp, phaseid!=1)
metasalt.sr.species <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.species <- subset(metasalt.sr.species, p.sr<0.05/531)

# salt-related gut-microbial functions ----------------
keggid <- names(alldata.mst)[grepl("kegg", names(alldata.mst))]
temp <- melt(alldata.mst, c("labid2","famid","phaseid"), keggid)
temp <- subset(temp, phaseid!=1)
metasalt.sr.kegg <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.kegg <- subset(metasalt.sr.kegg, p.sr<0.05/190)

cazyid <- names(alldata.mst)[grepl("cazy", names(alldata.mst))]
temp <- melt(alldata.mst, c("labid2","famid","phaseid"), cazyid)
temp <- subset(temp, phaseid!=1)
metasalt.sr.cazy <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.cazy <- subset(metasalt.sr.cazy, p.sr<0.05/133)

# salt-related metabolites ---------------------------
metid <- names(alldata.mst)[grepl("met", names(alldata.mst))]
temp <- melt(alldata.mst, c("labid2","famid","phaseid"), metid)
temp <- subset(temp, phaseid!=1)
metasalt.sr.metabol <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~phaseid+(1|labid2)+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  output = data.frame(coef.sr=coef, se.sr=se, p.sr=p)
  print(output)
})
metasalt.sr.metabol <- subset(metasalt.sr.metabol, p.sr<0.05/221)

### 2) - Identification of SSBP-related biomarkers ------------------------------------

# SSBP-related gut-microbial species ------------------------
speid <- names(chgdata.mst)[grepl("spe", names(chgdata.mst))]
temp <- melt(chgdata.mst, c("labid2","famid","ss3","age","gender","bmi","fcc","htn","smoking","chol"), speid)
metasalt.ss.species <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss3+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
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
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="ss3") %>% data.frame()
})

# SSBP-related gut-microbial function ----------------------
keggid <- names(chgdata.mst)[grepl("kegg", names(chgdata.mst))]
temp <- melt(chgdata.mst, c("labid2","famid","ss3","age","gender","bmi","fcc","htn","smoking","chol"), keggid)
metasalt.ss.kegg <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss3+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
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
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="ss3") %>% data.frame()
})

cazyid <- names(chgdata.mst)[grepl("cazy", names(chgdata.mst))]
temp <- melt(chgdata.mst, c("labid2","famid","ss3","age","gender","bmi","fcc","htn","smoking","chol"), cazyid)
metasalt.ss.cazy <- ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss3+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
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
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="ss3") %>% data.frame()
})

# SSBP-related metabolites ------------------------------------
metid <- names(chgdata.mst)[grepl("met", names(chgdata.mst))]
temp <- melt(chgdata.mst, c("labid2","famid","ss3","age","gender","bmi","fcc","htn","smoking","chol"), metid)
metasalt.ss.metabol <-  ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss3+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+bmi+fcc+htn+smoking+chol+(1|famid))
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
  fit <- lmer(data=temp, value~as.factor(ss3)+age+gender+fcc+bmi+smoking+chol+htn+(1|famid))
  ls <- emmeans(fit, specs="ss3") %>% data.frame()
})

# Associations between SSBP-related species and metabolites
speid <- names(chgdata.mst)[grepl("spe", names(chgdata.mst))]
metid <- names(chgdata.mst)[grepl("met", names(chgdata.mst))]
temp <- melt(chgdata.mst, c("labid2","famid","ss3","age","gender","bmi","fcc","htn","smoking","chol", metid), speid)
names(temp)[which(names(temp) %in% c("variable", "value"))] = c("speid", "spe.chg")
temp <- melt(temp, c("labid2","famid","ss3","age","gender","bmi","fcc","htn","smoking","chol", "speid", "spe.chg"), metid)
names(temp)[which(names(temp) %in% c("variable", "value"))] = c("metid", "met.chg")
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

## 1) Replication of salt-related metabolite ----------------------------
metid <- names(alldata.gst)[grepl("met", names(alldata.gst))]
temp <- melt(alldata.gst, c("labid2","phaseid"), metid)
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
metid <- names(chgdata.gst)[grepl("met", names(chgdata.gst))]
temp <- melt(chgdata.gst, c("labid2","matchid","ss3","bmi","fcc","smoking","chol"), metid)
gensalt.ss.metabol <-  ddply(temp, .(variable), function(temp){
  fit <- lmer(data=temp, value~ss3+bmi+fcc+smoking+chol+(1|matchid))
  coef = data.frame(summary(fit)$coefficients)$Estimate[2]
  se = data.frame(summary(fit)$coefficients)$Std..Error[2]
  p = data.frame(summary(fit)$coefficients)$Pr...t..[2]
  fit <- lmer(data=temp, value~as.factor(ss3)+bmi+fcc+smoking+chol+(1|matchid))
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
  fit <- lmer(data=temp, value~as.factor(ss3)+fcc+bmi+smoking+chol+(1|matchid))
  ls <- emmeans(fit, specs="ss3") %>% data.frame()
})

##### 3 - Validation in Rats ===========================================================================================

### 1) Wistar Rats --------------------------------------------------------------------------------

# Differences in BP levels and isovalerylcarnitine
fit <- aov(data=alldata.wistar, SBP~group*time+Error(ratid/time))
summary(fit)
fit <- aov(data=alldata.wistar, DBP~group*time+Error(ratid/time))
summary(fit)
temp <- unique(alldata.wistar[,c("group","isoval")])
t.test(data=temp, isoval~group, var.equal=T)

### 2) Dahl-SS rats ------------------------------------------------------------------------------

# Differences in BP levels
fit <- aov(data=alldata.dahl.bp, SBP~group*Time+Error(ratid/Time))
summary(fit)
fit <- aov(data=subset(alldata.dahl.bp, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")),  SBP~group*Time+Error(ratid/Time))
summary(fit)
t.test(data=subset(alldata.dahl.bp, Time==28 & group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), SBP~group, var.equal=T)

fit <- aov(data=alldata.dahl.bp, DBP~group*Time+Error(ratid/Time))
summary(fit)
fit <- aov(data=subset(alldata.dahl.bp, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")),  DBP~group*Time+Error(ratid/Time))
summary(fit)
t.test(data=subset(alldata.dahl.bp, Time==28 & group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), DBP~group, var.equal=T)

# Differences in isovalerylcarnitine 
t.test(data=subset(alldata.dahl.others, group %in% c("A: Control","B: HSD")), isoval~group, var.equal=T)
t.test(data=subset(alldata.dahl.others, group %in% c("A: Control","C: HSD + Isovalerylcarnitine")), isoval~group, var.equal=T)
t.test(data=subset(alldata.dahl.others, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), isoval~group, var.equal=T)

# Differences in the relaxation of mesenteric arteries
fit <- aov(data=alldata.dahl.artery, acety~group*con+Error(num/con))
summary(fit)
fit <- aov(data=subset(alldata.dahl.artery, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), acety~group*con+Error(num/con))
summary(fit)
fit <- aov(data=alldata.dahl.artery, nitrop~group*con+Error(num/con))
summary(fit)
fit <- aov(data=subset(alldata.dahl.artery, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), nitrop~group*con+Error(num/con))
summary(fit)

# Differences in the Fluorescence intensity of eNOS
t.test(data=subset(alldata.dahl.others, group %in% c("A: Control","B: HSD")), eNOS~group, var.equal=T)
t.test(data=subset(alldata.dahl.others, group %in% c("A: Control","C: HSD + Isovalerylcarnitine")), eNOS~group, var.equal=T)
t.test(data=subset(alldata.dahl.others, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), eNOS~group, var.equal=T)

# Differences in the Fibrotic area of aorta and kidney
wilcox.test(data=subset(alldata.dahl.others, group %in% c("A: Control","B: HSD")), Masson~group)
wilcox.test(data=subset(alldata.dahl.others, group %in% c("A: Control","C: HSD + Isovalerylcarnitine")), Masson~group)
wilcox.test(data=subset(alldata.dahl.others, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), Masson~group)

wilcox.test(data=subset(alldata.dahl.others, group %in% c("A: Control","B: HSD")), sirius~group)
wilcox.test(data=subset(alldata.dahl.others, group %in% c("A: Control","C: HSD + Isovalerylcarnitine")), sirius~group)
wilcox.test(data=subset(alldata.dahl.others, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")), sirius~group)

# Difference in BP levels of Dahl SS rats in the experiment of oral isovalerylcarnitine
fit <- aov(data=pharmadata, value~group*Time+Error(variable/Time))
summary(fit)

fit <- aov(data=alldata.dahl.oral, SBP~group*Time+Error(num/Time))
summary(fit)
fit <- aov(data=subset(alldata.dahl.oral, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")),  SBP~group*Time+Error(num/Time))
summary(fit)
t.test(data=subset(alldata.dahl.oral, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine") & Time==28), SBP~group, var.equal=T)

fit <- aov(data=alldata.dahl.oral, DBP~group*Time+Error(num/Time))
summary(fit)
fit <- aov(data=subset(alldata.dahl.oral, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine")),  DBP~group*Time+Error(num/Time))
summary(fit)
t.test(data=subset(alldata.dahl.oral, group %in% c("B: HSD","C: HSD + Isovalerylcarnitine") & Time=="28"), DBP~group, var.equal=T)

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
fit <- glm(data=cohort1, hyp~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine, family=binomial(link="logit"))
coef = data.frame(summary(fit)$coefficients)$Estimate[2]
se = data.frame(summary(fit)$coefficients)$Std..Error[2]
p = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- glm(data=cohort1, hyp~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine, family=binomial(link="logit"))
coef1 = data.frame(summary(fit)$coefficients)$Estimate[2]
se1 = data.frame(summary(fit)$coefficients)$Std..Error[2]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[2]
coef2 = data.frame(summary(fit)$coefficients)$Estimate[3]
se2 = data.frame(summary(fit)$coefficients)$Std..Error[3]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[3]

fit <- glm(data=cohort1, hyp~as.numeric(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine, family=binomial(link="logit"))
ptrend <- summary(fit)$coefficients[2,4]

AC.hyp1 = data.frame(outcome="Prevalent hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 2) Ideal BP status - incident prehypertension/hypertension --------------
fit <- coxph(data=cohort2, Surv(time, prehyp_hyp)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$coef[1]
se = data.frame(summary(fit)$coefficients)$se.coef.[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort2, Surv(time, prehyp_hyp)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$coef[1]
se1 = data.frame(summary(fit)$coefficients)$se.coef[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$coef[2]
se2 = data.frame(summary(fit)$coefficients)$se.coef[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort2, Surv(time, prehyp_hyp)~as.numeric(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp2 = data.frame(outcome="Ideal BP status - Prehypertension/hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 3) Prehypertension - incident hypertension ------------------------------
fit <- coxph(data=cohort3, Surv(time, inchyp)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$coef[1]
se = data.frame(summary(fit)$coefficients)$se.coef[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort3, Surv(time, inchyp)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$coef[1]
se1 = data.frame(summary(fit)$coefficients)$se.coef[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$coef[2]
se2 = data.frame(summary(fit)$coefficients)$se.coef[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort3, Surv(time, inchyp)~as.numeric(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp3 = data.frame(outcome="Prehypertension - Incident hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 4) Normotension - incident hypertension ---------------------------------
fit <- coxph(data=cohort4, Surv(time, inchyp2)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$coef[1]
se = data.frame(summary(fit)$coefficients)$se.coef[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort4, Surv(time, inchyp2)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$coef[1]
se1 = data.frame(summary(fit)$coefficients)$se.coef[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$coef[2]
se2 = data.frame(summary(fit)$coefficients)$se.coef[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort4, Surv(time, inchyp2)~as.numeric(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp4 = data.frame(outcome="Normotension - Incident hypertension", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

### 5) BP status increases >= 1 step ----------------------------------------
fit <- coxph(data=cohort5, Surv(time, probp)~isoval+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef = data.frame(summary(fit)$coefficients)$coef[1]
se = data.frame(summary(fit)$coefficients)$se.coef[1]
p = data.frame(summary(fit)$coefficients)$Pr...z..[1]

fit <- coxph(data=cohort5, Surv(time, probp)~as.factor(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
coef1 = data.frame(summary(fit)$coefficients)$coef[1]
se1 = data.frame(summary(fit)$coefficients)$se.coef[1]
p1 = data.frame(summary(fit)$coefficients)$Pr...z..[1]
coef2 = data.frame(summary(fit)$coefficients)$coef[2]
se2 = data.frame(summary(fit)$coefficients)$se.coef[2]
p2 = data.frame(summary(fit)$coefficients)$Pr...z..[2]

fit <- coxph(data=cohort5, Surv(time, probp)~as.numeric(isoval.group)+age+sex+bmi+area+region+diet_gdline+as.factor(smoke)+drink+edu+dm+dyslipid+work_pha+carnitine)
ptrend <- summary(fit)$coefficients[1,5]

AC.hyp5 = data.frame(outcome="BP status increases >= 1 step", coef=coef, se=se, p=p, coef1=coef1, se1=se1, p1=p1, coef2=coef2, se2=se2, p2=p2, ptrend=ptrend)

AC.hyp <- rbind(AC.hyp1, AC.hyp2, AC.hyp3, AC.hyp4, AC.hyp5)


