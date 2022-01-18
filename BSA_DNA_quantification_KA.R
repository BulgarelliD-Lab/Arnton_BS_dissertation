#############################################################
# Clean-up the memory and start a new session
#############################################################

#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

rm(list=ls())
dev.off()

##working directory

setwd("~/Dundee Year 4/Semester 1/R Directory")
getwd()

library ("ggplot2")
library("dplyr")
library("MASS")
library("ggfortify")
library("ggthemes")
library("colorspace")
library("grDevices")
library("vegan")
library ("ape")
library("PMCMRplus")
library("readxl") 
library("lsmeans")



###BSA_quantification of bacterial and fungal DNA


##########################################################
###############Plotting Log Bacteria 16S copies by variety
##########################################################

#import without ct outliers (over 33, end of calibration graph)
BacData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Organised Results Minus", range = "A1:E61")
View(BacData)

#import
BacData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Organised Results Clean", range = "A1:E72")
View(BacData)

#order the factor
BacData$Variety <-  ordered(BacData$Variety, levels=c( "HID144","Barke", "Int17","Int52", "Int19", "Int56", "Morex", "RGT", "Hockett", "GP", "FTII", "Bulk"))


p <-ggplot(BacData, aes(x=Variety, y=lne_16S, fill=Variety)) + geom_boxplot() #+ ggtitle(" Natural Log Bacterial 16S copies")#+ylim(5,25)
#p + geom_jitter( size=5,shape=21, position=position_jitter(0))#+ scale_fill_manual(values = c ("#000000","#0072B2","#56B4E9", "#E69F00", "#D55E00","#D55E00"))
p + ggtitle("Natural log of 16S Copies by Variety") +
  xlab("Variety") + ylab("lne(16S)") +
  geom_jitter( size=3,shape=20, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

#Organise by P/A

BacData$PA <-  ordered(BacData$PA, levels=c( "P","A","N/A")) 

BacData$Variety <-  ordered(BacData$Variety, levels=c( "HID144","Barke", "Int17","Int52", "Hockett", "GP", "Int19", "Int56", "Morex", "RGT","FTII", "Bulk"))

p <-ggplot(BacData, aes(x=Variety, y=lne_16S, fill=PA)) + geom_boxplot() 
p + ggtitle("Natural log of 16S Copies by Variety & NLR Presence/Absence") +
  xlab("Variety") + ylab("lne(16S)") +
  geom_jitter( size=3,shape=20, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)+
  #scale_fill_viridis_d()+
  scale_fill_manual(values=c("#648FFF", "#DC267F", "#FFB000"))+
  labs(fill='NLR\nPresence/\nAbsence') 

##################################################################
#Stats Bacterial by variety
shapiro.test(BacData$lne_16S)

kruskal.test(lne_16S~Variety , data=BacData)

kwAllPairsDunnTest(x= BacData$lne_16S, g= BacData$Variety, p.adjust.method="BH")










###################################################################
###############Plotting Log Bacteria 16S copies by Presence/Absence
###################################################################

#import without ct outliers (over 33, end of calibration graph)
BacData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Organised Results Minus", range = "A1:E61")
View(BacData)

#import
BacData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Results Clean", range = "A1:E72")
View(BacData)

#remove bulk
Rhizo<-subset(BacData, Variety!= "Bulk")

#Order
Rhizo$PA <-  ordered(Rhizo$PA, levels=c( "P","A")) 

#plot
s <-ggplot(Rhizo, aes(x=PA, y=lne_16S, fill=PA)) + geom_boxplot() 
s + ggtitle("Natural log of 16S Copies by Presence/Absence") +
  xlab("Presence/Absence") + ylab("lne(16S)") +
  #scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values=c("#648FFF", "#DC267F", "#FFB000"))+
  labs(fill='NLR\nPresence/\nAbsence')+
  geom_jitter( size=3,shape=20, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)
##################################################################
#Stats Bacterial by P/A
shapiro.test(Rhizo$lne_16S)

#kruskal.test(lne_16S~PA , data=Rhizo)

#kwAllPairsDunnTest(x= Rhizo$lne_16S, g= Rhizo$PA, p.adjust.method="BH")

wilcox.test(lne_16S~PA , data=Rhizo)











#######################################################################
###############plotting concentrations Log Fungal ITS copies by variety
#######################################################################

#Import and view data
FungalData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Liam's Results Clean", range = "A1:C65")
View(FungalData)


FungalData$Variety <-  ordered(FungalData$Variety, levels=c( "HID144","Barke", "Int17","Int52", "Int19", "Int56", "Morex", "RGT", "Hockett", "GP", "Bulk"))


q <-ggplot(FungalData, aes(x=Variety, y=lne_ITS, fill=Variety)) + geom_boxplot() + ggtitle(" Natural Log Fungal ITS copies")#+ylim(5,25)
#q + geom_jitter( size=5,shape=21, position=position_jitter(0))#+ scale_fill_manual(values = c ("#000000","#0072B2","#56B4E9", "#E69F00", "#D55E00","#D55E00"))
q+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  geom_jitter( size=5,shape=20, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

##################################################################
#Stats Fungal by variety
shapiro.test(FungalData$lne_ITS)

kruskal.test(lne_ITS~Variety , data=FungalData)

kwAllPairsDunnTest(x= FungalData$lne_ITS, g= FungalData$Variety, p.adjust.method="BH")











#############################################################
##Dry Weight by Variety
#############################################################

#import
WeightData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Organised Results Clean", range = "A1:E72")
View(WeightData)

RhizoW<-subset(WeightData, Variety!= "Bulk" & Sample!=128)

#Order the levels according to a defined order 
RhizoW$Variety <-  ordered(RhizoW$Variety, levels=c( "HID144","Barke", "Int17","Int52", "Int19", "Int56", "Morex", "RGT", "Hockett", "GP", "FTII", "Bulk")) 


shapiro.test(RhizoW$DryWeight)

#If normal
#dry_w_results_aov<-aov(DryWeight~Variety, data=RhizoW)
#summary(dry_w_results_aov)


#dry_w_results_tukey.test <- TukeyHSD(dry_w_results_aov)
#dry_w_results_tukey.test


#If non-normal
kruskal.test(DryWeight~Variety , data=RhizoW)

kwAllPairsDunnTest(x= RhizoW$DryWeight, g= RhizoW$Variety, p.adjust.method="BH")



r<-ggplot(RhizoW, aes(x=Variety, y=DryWeight, fill=Variety))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


#r+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  #geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

r + ggtitle("Dry Weight by Variety") +
  xlab("Variety") + ylab("Dry Weight (g)") +
  scale_fill_viridis_d()


#Also by P/A
RhizoW$PA <-  ordered(RhizoW$PA, levels=c( "P","A")) 

RhizoW$Variety <-  ordered(RhizoW$Variety, levels=c( "HID144","FTII","Int17","Int52", "Int19", "Int56", "Barke", "Hockett", "GP", "Morex", "RGT", "Bulk"))

p <-ggplot(RhizoW, aes(x=Variety, y=DryWeight, fill=PA)) + geom_boxplot() 
p + ggtitle("Dry Weight by Variety & NLR Presence/Absence") +
  xlab("Variety") + ylab("Dry Weight (g)") +
  geom_jitter( size=3,shape=20, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)+
  scale_fill_manual(values=c("#648FFF", "#DC267F", "#FFB000"))+
  labs(fill='NLR\nPresence/\nAbsence') 





#############################################################
##Dry Weight by presence/absence
#############################################################

#import
WeightData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Organised Results Clean", range = "A1:E72")
View(WeightData)

RhizoW<-subset(WeightData, Variety!= "Bulk" & Sample!=128)

#Order the levels according to a defined order 
RhizoW$PA <-  ordered(RhizoW$PA, levels=c( "P","A")) 

shapiro.test(RhizoW$DryWeight)

#If normal
#PA_dry_w_results_aov<-aov(DryWeight~PA, data=RhizoW)
#summary(PA_dry_w_results_aov)


#PA_dry_w_results_tukey.test <- TukeyHSD(PA_dry_w_results_aov)
#PA_dry_w_results_tukey.test

#If non-normal
#kruskal.test(DryWeight~PA , data=RhizoW)

#kwAllPairsDunnTest(x= RhizoW$DryWeight, g= RhizoW$PA, p.adjust.method="BH")

#Wait, shouldn't this be a 2 factor?
wilcox.test(DryWeight~PA , data=RhizoW)


t<-ggplot(RhizoW, aes(x=PA, y=DryWeight, fill=PA))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


#t+#scale_fill_manual(values=c( "#56B4E9","#E69F00"))+  
  #geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

t + ggtitle("Dry Weight by Presence/Absence") +
  xlab("Presence/Absence") + ylab("Dry Weight (g)") +
  #geom_boxplot(alpha=0.3) +
  #theme(legend.position="none") +
  labs(fill='NLR\nPresence/\nAbsence')+
  #scale_fill_brewer(palette="Set2")
  scale_fill_manual(values=c("#648FFF", "#DC267F", "#FFB000"))



#############################################################
##Dry Weight by presence/absence minus wild varieties
#############################################################

#import
WeightData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Organised Results Clean", range = "A1:E72")
View(WeightData)

EliteRhizoW<-subset(WeightData, Variety!="HID144" & Variety!="Int17" & Variety!="Int52" & Variety!="Int19" & Variety!="Int56" & Variety!="FTII" & Variety!="Bulk" & Sample!=128)

#Order the levels according to a defined order 
EliteRhizoW$PA <-  ordered(EliteRhizoW$PA, levels=c( "P","A")) 

shapiro.test(EliteRhizoW$DryWeight)

#As non normal
wilcox.test(DryWeight~PA , data=EliteRhizoW)

#wild
WildRhizoW<-subset(WeightData, Variety!="Barke" & Variety!="Int17" & Variety!="Int52" & Variety!="Int19" & Variety!="Int56" & Variety!="Hockett" & Variety!="Bulk" & Variety!="GP" & Variety!="Morex" & Variety!="RGT"& Sample!=128)
WildRhizoW$PA <-  ordered(WildRhizoW$PA, levels=c( "P","A")) 
shapiro.test(WildRhizoW$DryWeight)
t.test(DryWeight~PA , data=WildRhizoW)

#Int
IntRhizoW<-subset(WeightData, Variety!="Barke" & Variety!="FTII" & Variety!="HID144" & Variety!="Hockett" & Variety!="Bulk" & Variety!="GP" & Variety!="Morex" & Variety!="RGT"& Sample!=128)
IntRhizoW$PA <-  ordered(IntRhizoW$PA, levels=c( "P","A")) 
shapiro.test(IntRhizoW$DryWeight)
t.test(DryWeight~PA , data=IntRhizoW)


b<-ggplot(EliteRhizoW, aes(x=PA, y=DryWeight, fill=PA))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


#b+#scale_fill_manual(values=c( "#56B4E9","#E69F00"))+  
 # geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

b + ggtitle("Dry Weight by Presence/Absence in Elite Varieties") +
  xlab("Presence/Absence") + ylab("Dry Weight (g)") +
  #geom_boxplot(alpha=0.3) +
  #theme(legend.position="none") +
  scale_fill_brewer(palette="Set2")





#############################################################
#Nitrogen Experiment
#############################################################
##Dry Weight by Variety
#############################################################

#import
NData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Nitrogen", range = "A1:D36")
View(NData)

#Order the levels according to a defined order 
NData$Variety <-  ordered(NData$Variety, levels=c( "HOR7552","HOR10350", "HOR13942","HOR3365", "HOR8148", "HOR9043")) 


shapiro.test(NData$DryWeight)

#If normal
#NV_dry_w_results_aov<-aov(DryWeight~Variety, data=NData)
#summary(NV_dry_w_results_aov)


#NV_dry_w_results_tukey.test <- TukeyHSD(NV_dry_w_results_aov)
#NV_dry_w_results_tukey.test

#If non-normal
kruskal.test(DryWeight~Variety , data=NData)

kwAllPairsDunnTest(x= NData$DryWeight, g= NData$Variety, p.adjust.method="BH")





m<-ggplot(NData, aes(x=Variety, y=DryWeight, fill=Variety))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


m+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

m + ggtitle("Dry Weight by Variety") +
  xlab("Variety") + ylab("Dry Weight (g)") +
  scale_fill_viridis_d()








#############################################################
##Dry Weight by N Treatment
#############################################################

#import
NData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Nitrogen", range = "A1:D36")
View(NData)

#Order the levels according to a defined order 
NData$Treatment <-  ordered(NData$Treatment, levels=c( "N100","N0", "DW")) 


shapiro.test(NData$DryWeight)

#If normal
#N_dry_w_results_aov<-aov(DryWeight~Treatment, data=NData)
#summary(N_dry_w_results_aov)


#N_dry_w_results_tukey.test <- TukeyHSD(N_dry_w_results_aov)
#N_dry_w_results_tukey.test

#If non-normal
kruskal.test(DryWeight~Treatment , data=NData)

kwAllPairsDunnTest(x= NData$DryWeight, g= NData$Treatment, p.adjust.method="BH")




n<-ggplot(NData, aes(x=Treatment, y=DryWeight, fill=Treatment))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


n+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

n + ggtitle("Dry Weight by Treatment") +
  xlab("Treatment") + ylab("Dry Weight (g)") +
  scale_fill_viridis_d()

###########################################################
#Separate by Variety

#Order
NData$Variety <-  ordered(NData$Variety, levels=c( "HOR7552","HOR10350", "HOR13942","HOR3365", "HOR8148", "HOR9043")) 
NData$Treatment <-  ordered(NData$Treatment, levels=c( "N100","N0", "DW")) 


w<-ggplot(NData, aes(x=Variety, y=DryWeight, fill=Treatment))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


#w+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  #geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

w + ggtitle("Dry Weight by Treatment & Variety") +
  xlab("Variety") + ylab("Dry Weight (g)") +
  scale_fill_viridis_d()

#Analysis

#Overall


#By subsets
ND7552 <- subset(NData, Variety== "HOR7552")
shapiro.test(ND7552$DryWeight)
#If non normal
#kruskal.test(DryWeight~Treatment , data=ND7552)
#kwAllPairsDunnTest(x= ND7552$DryWeight, g= ND7552$Treatment, p.adjust.method="BH")
#If normal
ND7552aov<-aov(DryWeight~Treatment, data=ND7552)
summary(ND7552aov)
TukeyHSD(ND7552aov)

ND10350 <- subset(NData, Variety== "HOR10350")
shapiro.test(ND10350$DryWeight)
#If non normal
#kruskal.test(DryWeight~Treatment , data=ND10350)
#kwAllPairsDunnTest(x= ND10350$DryWeight, g= ND10350$Treatment, p.adjust.method="BH")
#If normal
ND10350aov<-aov(DryWeight~Treatment, data=ND10350)
summary(ND10350aov)
TukeyHSD(ND10350aov)

ND13942 <- subset(NData, Variety== "HOR13942")
shapiro.test(ND13942$DryWeight)
#If non normal
#kruskal.test(DryWeight~Treatment , data=ND13942)
#kwAllPairsDunnTest(x= ND13942$DryWeight, g= ND13942$Treatment, p.adjust.method="BH")
#If normal
ND13942aov<-aov(DryWeight~Treatment, data=ND13942)
summary(ND13942aov)
TukeyHSD(ND13942aov)

ND3365 <- subset(NData, Variety== "HOR3365")
shapiro.test(ND3365$DryWeight)
#If non normal
#kruskal.test(DryWeight~Treatment , data=ND3365)
#kwAllPairsDunnTest(x= ND3365$DryWeight, g= ND3365$Treatment, p.adjust.method="BH")
#If normal
ND3365aov<-aov(DryWeight~Treatment, data=ND3365)
summary(ND3365aov)
TukeyHSD(ND3365aov)

ND8148 <- subset(NData, Variety== "HOR8148")
shapiro.test(ND8148$DryWeight)
#If non normal
#kruskal.test(DryWeight~Treatment , data=ND8148)
#kwAllPairsDunnTest(x= ND8148$DryWeight, g= ND8148$Treatment, p.adjust.method="BH")
#If normal
ND8148aov<-aov(DryWeight~Treatment, data=ND8148)
summary(ND8148aov)
TukeyHSD(ND8148aov)

ND9043 <- subset(NData, Variety== "HOR9043")
shapiro.test(ND9043$DryWeight)
kruskal.test(DryWeight~Treatment , data=ND9043)
kwAllPairsDunnTest(x= ND9043$DryWeight, g= ND9043$Treatment, p.adjust.method="BH")
#If normal
#ND9043aov<-aov(DryWeight~Treatment, data=ND9043)
#summary(ND9043aov)
#TukeyHSD(ND9043aov)



##############################################
#Dry Weight by Soil
##############################################

#import
SData <- read_excel("~/Dundee Year 4/Semester 1/BSA_bacterial_quantification_181121_data_organised.xlsx", sheet = "Soil", range = "A1:D72")
View(SData)

#Order the levels according to a defined order 
SData$Variety <-  ordered(SData$Variety, levels=c( "HID144","Barke", "Int17","Int52", "Int19", "Int56")) 
SData$Soil <-  ordered(SData$Soil, levels=c( "Q3","Q4")) 

shapiro.test(SData$DryWeight)

#As non-normal
wilcox.test(DryWeight~Soil , data=SData)

#Per Variety
SDH <- subset(SData, Variety== "HID144")
shapiro.test(SDH$DryWeight)
wilcox.test(DryWeight~Soil , data=SDH)

SDB <- subset(SData, Variety== "Barke")
shapiro.test(SDB$DryWeight)
wilcox.test(DryWeight~Soil , data=SDB)

SD17 <- subset(SData, Variety== "Int17")
shapiro.test(SD17$DryWeight)
wilcox.test(DryWeight~Soil , data=SD17)

SD52 <- subset(SData, Variety== "Int52")
shapiro.test(SD52$DryWeight)
wilcox.test(DryWeight~Soil , data=SD52)

SD19 <- subset(SData, Variety== "Int19")
shapiro.test(SD19$DryWeight)
wilcox.test(DryWeight~Soil , data=SD19)

SD56 <- subset(SData, Variety== "Int56")
shapiro.test(SD56$DryWeight)
t.test(DryWeight~Soil , data=SD56)



#For Soil and Variety
f<-ggplot(SData, aes(x=Variety, y=DryWeight, fill=Soil))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


f+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

f + ggtitle("Dry Weight by Variety & Soil") +
  xlab("Variety") + ylab("Dry Weight (g)") +
  scale_fill_viridis_d()


#For Soil alone
h<-ggplot(SData, aes(x=Soil, y=DryWeight, fill=Soil))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


h+#scale_fill_manual(values=c( "#999999", "#0072B2","#56B4E9","#E69F00","#D55E00","#D55E00"))+  
  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

h + ggtitle("Dry Weight by Soil") +
  xlab("Variety") + ylab("Dry Weight (g)") +
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set2")








