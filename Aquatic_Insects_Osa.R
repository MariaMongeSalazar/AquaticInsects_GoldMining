#1 General ------------------------------------------------------
#1.1 Clean the environment --------------------------------------
rm(list=ls())

#1.2 Set working directory --------------------------------------
setwd("C:/Users/Amparo/Desktop/BC/Analisis")

#1.3 Install packages -------------------------------------------
install.packages("vegan") #For ecological analysis
install.packages("ggplot2") #To create aesthetic graphs
install.packages("ggrepel") #to repel text on ggplot
install.packages("ellipse") #To draw ellipses with 3 points
install.packages("lmer") #For glmm and lmm
install.packages("r2glmm") #To get the R2 for glmm

#1.4 Load packages ----------------------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(ellipse)
library(lme4)
library(r2glmm)

#1.5 Import data ------------------------------------------------
General<-read.csv("General.csv", header=TRUE)
Richness<-read.csv("Richness.csv", header=TRUE)
Abundance<-read.csv("Abundance.csv", header=TRUE)
Abundance$Stream<-as.factor(Abundance$Stream)

#2 NMDS ---------------------------------------------------------
#NMDS
nmds<-metaMDS(Abundance[,6:86],distance="bray", k=2, trymax=10000)

#Create table with coordinates
NMDS_table <- data.frame(MDS1 = nmds$points[,1],
                         MDS2 = nmds$points[,2], 
                         group=Abundance$Level)

#Create data frame with environmental loadings
env_var<- envfit(nmds, Abundance[,c(4,5)], na.rm=TRUE,
                 perm = 10000) 
env_scrs<-as.data.frame(scores(env_var,display="vectors"))
env_var #Check the significant variables

#Create hull values
Hull.B <-NMDS_table[NMDS_table$group == "Before", ][chull
  (NMDS_table[NMDS_table$group =="Before",c("MDS1", "MDS2")]), ]
Hull.GM <-NMDS_table[NMDS_table$group == "Gold mining", ][chull
 (NMDS_table[NMDS_table$group =="Gold mining",c("MDS1", "MDS2")]), ]
Hull.T <-NMDS_table[NMDS_table$group == "Tributary", ][chull
 (NMDS_table[NMDS_table$group =="Tributary",c("MDS1", "MDS2")]), ]
Hull.C <-NMDS_table[NMDS_table$group == "Confluence", ][chull
 (NMDS_table[NMDS_table$group =="Confluence",c("MDS1", "MDS2")]), ]

#Create Hull data frame
hull.data <- rbind(Hull.B, Hull.GM, Hull.T, Hull.C)  #combine hulls
hull.data

#Plot
ggplot(NMDS_table, aes(x=MDS1, y=MDS2, col=group)) +
  geom_polygon(data=hull.data,aes(x=MDS1,y=MDS2,fill=group,group=group),alpha=0.30)+
  geom_point() +
  theme_classic()

##NMDS According to Stream
#NMDS
nmds<-metaMDS(Abundance[,6:86],distance="bray", k=2, trymax=10000)

#Create table with coordinates
NMDS_table <- data.frame(MDS1 = nmds$points[,1],
                         MDS2 = nmds$points[,2], 
                         group=Abundance$Stream)

#Create data frame with environmental loadings
env_var<- envfit(nmds, Abundance[,c(4,5)], na.rm=TRUE,
                 perm = 10000) 
env_scrs<-as.data.frame(scores(env_var,display="vectors"))
env_var #Check the significant variables

#Create hull values
Hull.1 <-NMDS_table[NMDS_table$group == "1", ][chull
  (NMDS_table[NMDS_table$group =="1",c("MDS1", "MDS2")]), ]
Hull.2 <-NMDS_table[NMDS_table$group == "2", ][chull
  (NMDS_table[NMDS_table$group =="2",c("MDS1", "MDS2")]), ]

#Create Hull data frame
hull.data <- rbind(Hull.1, Hull.2)  #combine hulls
hull.data

#Plot
ggplot(NMDS_table, aes(x=MDS1, y=MDS2, col=group)) +
  geom_polygon(data=hull.data,aes(x=MDS1,y=MDS2,fill=group,group=group),alpha=0.30)+
  geom_point() +
  theme_classic()

#3 Permanova ----------------------------------------------------
#Dissimilarity matrix
Bray<-vegdist(Abundance[,6:86], 
              method="bray", binary=FALSE, diag=FALSE, 
              upper=FALSE, na.rm=FALSE)

#Permanova
Permanova<-adonis(Bray~Level+Stream, 
                  data=Abundance, 
                  na.rm=TRUE, permutations=10000, 
                  method="bray", strata=NULL)
densityplot(permustats(Permanova)) #Check plots
Permanova # Check results

#4 Simper -------------------------------------------------------
Simp<-simper(Abundance[,6:86], group=Abundance$Level, 
             permutations=0, trace=FALSE)
summary(Simp, ordered=TRUE) 
#Check the overall contribution in the list of results

Simp<-simper(Abundance[,6:86], group=Abundance$Stream, 
             permutations=0, trace=FALSE)
summary(Simp, ordered=TRUE)

#5 GLM ----------------------------------------------------------
# Check for ramdom variables
#Transect variance / (transect variance + residual varaince)*100
MOD1<-lmer(Richness~1+(1|Replicate), data=General)
summary(MOD1)
2.618/(2.618+53.351)*100 #4.67% of variance explained

MOD1<-lmer(Abundance~1+(1|Replicate), data=General)
summary(MOD1) #0% of variance explained

MOD1<-lmer(Richness~1+(1|Stream), data=General)
summary(MOD1)
9.601/(9.601+50.163)*100 #16% of variance explained

MOD1<-lmer(Abundance~1+(1|Stream), data=General)
summary(MOD1)
882.5/(882.5+5324.6)*100 #14% of variance explained

# Models
GLMA<-glmer(Abundance~Area+Microhabitats+ Level+ (1|Stream), 
            family=poisson, data=General)
summary(GLMA)

GLMR<-glmer(Richness~Area+Microhabitats+ Level+ (1|Stream), 
            family=poisson, data=General)
summary(GLMR)

# Model assumptions
##Normality of the residuals
qqnorm(residuals(GLMA))
qqnorm(residuals(GLMR))

##Homogeneity of variances
plot(GLMA)
plot(GLMR)

## Goodness-of-fit
r2beta(GLMA, partial = TRUE, method = "sgv", data = NULL)
r2beta(GLMR, partial = TRUE, method = "sgv", data = NULL)

#Effect sizes
##Abundance
#Stream= 4% of Variance
exp(0.04536) #Area = 1.046404
exp(0.02668) #Microhabitats = 1.027039
exp(-0.17272) #Confluence = 0.8413732
exp(-1.08605) #Gold Mining = 0.3375472
exp(-0.11674) #Tributary = 0.8898165

##Richness
#Stream=0% of variance
exp(0.13066) #Area = 1.13958
exp(0.09316) #Microhabitats = 1.097637
exp(0.10886) #Confluence = 1.115006
exp(-0.55664) #Gold Mining = 0.5731316
exp(0.10522) #Tributary = 1.110955

#Plots
## Abundance / Area
ggplot(General,aes(x=Area, y=Abundance))+ 
  geom_point(shape=19, color="royalblue2")+
  theme_classic()+
  labs(x="Area (m2)", y="Abundance")+
  geom_smooth(method=lm, linetype="solid", color="royalblue3", 
              fill="grey70")+ ggtitle("A")+
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size = 12),
        plot.title = element_text(colour="black", size=18))

## Abundance / Microhabitats
ggplot(General,aes(x=Microhabitats, y=Abundance))+ 
  geom_point(shape=19, color="royalblue2")+
  theme_classic()+
  labs(x="Number of microhabitats", y="Abundance")+
  geom_smooth(method=lm, linetype="solid", color="royalblue3", 
              fill="grey70")+ ggtitle("A")+
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size = 12),
        plot.title = element_text(colour="black", size=18))

## Richness / Area
ggplot(General,aes(x=Area, y=Richness))+ 
  geom_point(shape=19, color="royalblue2")+
  theme_classic()+
  labs(x="Area (m2)", y="Richness")+
  geom_smooth(method=lm, linetype="solid", color="royalblue3", 
              fill="grey70")+ ggtitle("A")+
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size = 12),
        plot.title = element_text(colour="black", size=18))

## Richness / Microhabitats
ggplot(General,aes(x=Microhabitats, y=Richness))+ 
  geom_point(shape=19, color="royalblue2")+
  theme_classic()+
  labs(x="Number of microhabitats", y="Richness")+
  geom_smooth(method=lm, linetype="solid", color="royalblue3", 
              fill="grey70")+ ggtitle("A")+
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size = 12),
        plot.title = element_text(colour="black", size=18))

# 6 Kruskal-Wallis
kruskal.test(Area ~ Level, data = General)
kruskal.test(Microhabitats ~ Level, data = General)
kruskal.test(Area ~ Stream, data = General)
kruskal.test(Microhabitats ~ Stream, data = General)
