###############################################################################
#                             PACKAGES UTILES                                 #
###############################################################################
rm(list=ls(all=TRUE))
setwd("C:/Users/Cassandre/Documents/Stage/Stage/Transfert_APAD_IDESP_2021")

library(dplyr)
library(ggplot2)
library(nlme)
library(AICcmodavg)
library(stringr)
library(knitr)
library(lcmm)
library(gridExtra)
library(splines)
library(lmtest)

#--TABLE-----------------------------------------------------------------------
Data <- read.csv("QLQC30_score.csv", header = T, sep = ";")
df <- subset(Data, select = -c(q1:q30))
df <- df %>% relocate(QL2, .after = qsdtcd1)

#--PARTICIPANTES---------------------------------------------------------------
N <- length(unique(df$npat)) #nbr de patientes ITT
N_inclusion <- length(filter(df, visit==0)$npat) #nbr de patientes à l'inclusion

npat_ITT <- unique(df$npat) #nbr de patientes ITT
npat_mITT<- filter(df, visit==0)$npat
which(!npat_ITT %in% npat_mITT) #patiente qui n'est pas à l'inclusion (62)

###############################################################################
#                         CALCUL VRAIES DATES (mois)                          #
###############################################################################
df$qsdtcd1 <- as.Date(df$qsdtcd1, format = "%d/%m/%Y")

#--Patiente 62-----------------------------------------------------------------
npat_62 <- data.frame(matrix("",1,19))
npat_62 <- data.frame(npat=c(62), BRAS=c(0), visit=c(0), 
                      qsdtcd1=as.Date("2012-03-21"), QL2=c(NA), PF2=c(NA),
                      RF2=c(NA), EF=c(NA), CF=c(NA), SF=c(NA), FA=c(NA), 
                      NV=c(NA), PA=c(NA), DY=c(NA), SL=c(NA), AP=c(NA), 
                      CO=c(NA), DI=c(NA), FI=c(NA))
npat_62 <- rbind(npat_62, filter(df, npat=="62"))
time_p62 <- npat_62 %>% group_by(npat) %>% 
  mutate(LagDate = qsdtcd1[visit=="0"], # ci-dessous : conversion en mois
         Time = as.numeric((difftime(qsdtcd1, LagDate, units="days"))/30.4375))

#--Toutes les autres patientes-------------------------------------------------
tps <- df[-287,]
which(is.na(tps$qsdtcd1)) #val: 253 / 308 / 481 / 485 / 658
tps$qsdtcd1[253] <- as.Date("2012-05-24") #date trouvée fichier VS excel
tps$qsdtcd1[308] <- as.Date("2013-10-21") #date trouvée fichier VS excel
tps$qsdtcd1[481] <- as.Date("2014-04-07") #date trouvée fichier VS excel
tps$qsdtcd1[485] <- as.Date("2014-04-16") #date trouvée fichier VS excel 
tps$qsdtcd1[658] <- as.Date("2014-11-10") #date trouvée fichier VS excel
which(is.na(tps$qsdtcd1)) #verifier qu'il n'y a plus de NA

temps <- tps %>% group_by(npat) %>% 
  mutate(LagDate = qsdtcd1[visit=="0"], #ci-dessous : conversion en mois
         Time = as.numeric((difftime(qsdtcd1, LagDate, units="days"))/30.4375))

#--tableau total----------------------------------------------
df_tpreel <- bind_rows(temps, time_p62)
df_tpreel <- df_tpreel[order(df_tpreel$npat),]
df_tpreel <- df_tpreel[-287,] #enleve la ligne creee pour pat 62 (inclusion)
df_tpreel <- df_tpreel %>% relocate(LagDate, .after = qsdtcd1)
df_tpreel <- df_tpreel %>% relocate(Time, .after = LagDate)
df_tpreel <- as.data.frame(df_tpreel)
df_tpreel[,"BRAS"] <- as.factor(df_tpreel$BRAS)

# write.table(df_tpreel, file = "tableau_donnees_mois.csv", sep = ";",
#             dec = ",", row.names = FALSE)

###############################################################################
#                        CALCUL DES TAUX DE COMPLIANCE                        #
###############################################################################
#--Taleau sans doublons
df2 <- df[, 1:2] #tableau avec que numero et bras
doublons <- which(duplicated(df2$npat)) #cherche les doublons
df3 <- df2[-doublons,] #enleve les doublons

#--Effectif des 2 bras
bras0 <- sum(df3$BRAS==0) #nb de patientes en groupe controle
bras1 <- sum(df3$BRAS==1) #nb de patientes en groupe APAD

#--------BRAS 0 (controle) & chaque visites-------------------
arm0v0 <- length(which(df$visit==0 & df$BRAS==0)) #Baseline
arm0v1 <- length(which(df$visit==18 & df$BRAS==0)) #18 weeks
arm0v2 <- length(which(df$visit==27 & df$BRAS==0)) #27 weeks
arm0v3 <- length(which(df$visit==54 & df$BRAS==0)) #54 weeks
arm0v4 <- length(which(df$visit==81 & df$BRAS==0)) #81 weeks

#--------BRAS 1 (APAD) & chaque visites-----------------------
arm1v0 <- length(which(df$visit==0 & df$BRAS==1)) #Baseline
arm1v1 <- length(which(df$visit==18 & df$BRAS==1)) #18 weeks
arm1v2 <- length(which(df$visit==27 & df$BRAS==1)) #27 weeks
arm1v3 <- length(which(df$visit==54 & df$BRAS==1)) #54 weeks
arm1v4 <- length(which(df$visit==81 & df$BRAS==1)) #81 weeks

#--Taux compliances 
B0_baseline <- arm0v0/bras0
B0_W18 <- arm0v1/bras0
B0_W27 <- arm0v2/bras0
B0_W54 <- arm0v3/bras0
B0_W81 <- arm0v4/bras0

B1_baseline <- arm1v0/bras1
B1_W18 <- arm1v1/bras1
B1_W27 <- arm1v2/bras1
B1_W54 <- arm1v3/bras1
B1_W81 <- arm1v4/bras1

###############################################################################
#                TEST DE STUDENT POUR CHAQUE DIMENSION                        #
###############################################################################

#~~~~~~~~~~~~~~Statut de sante gobale~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t.test(df_tpreel$QL2[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$QL2[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#~~~~~~~~~~~~~~Echelles Fonctionnelles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Physique
t.test(df_tpreel$PF2[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$PF2[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Personnelle
t.test(df_tpreel$RF2[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$RF2[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Emotionnelle
t.test(df_tpreel$EF[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$EF[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Cognitive
t.test(df_tpreel$CF[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$CF[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Sociale
t.test(df_tpreel$SF[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$SF[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#~~~~~~~~~~~~~~Echelles Symptomatiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Fatigue
t.test(df_tpreel$FA[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$FA[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Nausée et vomissement
t.test(df_tpreel$NV[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$NV[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Douleur
t.test(df_tpreel$PA[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$PA[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Dyspnée
t.test(df_tpreel$DY[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$DY[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Insomnie
t.test(df_tpreel$SL[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$SL[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Perte d'appétit
t.test(df_tpreel$AP[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$AP[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Constipation
t.test(df_tpreel$CO[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$CO[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Diarrhée
t.test(df_tpreel$DI[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$DI[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

#Difficultés financière
t.test(df_tpreel$FI[which(df_tpreel$visit==0 & df_tpreel$BRAS==0)], 
       df_tpreel$FI[which(df_tpreel$visit==0 & df_tpreel$BRAS==1)])

###############################################################################
#                          LMMs (LINEAR MIXED MODELS)                         #
###############################################################################

#------------------------------------------------------------------------------
#         trajectoire linéaire: intercept et pente aléatoires
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~~Statut de sante gobale~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.lme.QL2 <- lme(QL2 ~ Time + Time:BRAS, 
                   random = ~ Time|npat, data=df_tpreel, 
                   na.action = na.omit, method="REML")
summary(fit.lme.QL2)

#~~~~~~~~~~~~~~Echelles Fonctionnelles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fonctions physiques
fit.lme.PF <- lme(PF2 ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.PF)

#Fonctions personnelles
fit.lme.RF <- lme(RF2 ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.RF)

#Fonctions emotionnelles
fit.lme.EF <- lme(EF ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML", 
                  control = lmeControl(opt="optim"))
summary(fit.lme.EF)

#Fonctions cognitives
fit.lme.CF <- lme(CF ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.CF)

#Fonctions sociales
fit.lme.SF <- lme(SF ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML",
                  control = lmeControl(opt="optim"))
summary(fit.lme.SF)

#~~~~~~~~~~~~~~Echelles Symptomatiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fatigue
fit.lme.FA <- lme(FA ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.FA)

#Nausées et vomissements
fit.lme.NV <- lme(NV ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML", 
                  control = lmeControl(opt="optim"))
summary(fit.lme.NV)

#Douleur
fit.lme.PA <- lme(PA ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.PA)

#Dyspnée
fit.lme.DY <- lme(DY ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.DY)

#Insomnie
fit.lme.SL <- lme(SL ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.SL)

#Perte d'appetit
fit.lme.AP <- lme(AP ~ Time + Time:BRAS,  
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML", 
                  control = lmeControl(opt="optim"))
summary(fit.lme.AP)

#Constipation
fit.lme.CO <- lme(CO ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML")
summary(fit.lme.CO)

#Diarrhées
fit.lme.DI <- lme(DI ~ Time + Time:BRAS, 
                  random = ~ Time|npat, data=df_tpreel, 
                  na.action = na.omit, method="REML",
                  control = lmeControl(opt="optim"))
summary(fit.lme.DI)

#Difficultés financières
fit.lme.FI <- lme(FI ~ Time + Time:BRAS,
                  random = ~ Time|npat, data=df_tpreel,
                  na.action = na.omit, method="REML",
                  control = lmeControl(opt="optim"))
summary(fit.lme.FI)

#------------------------------------------------------------------------------
#             Tableau des estimateurs beta1 et beta2 & p-value
#           Calcul des beta1+beta2 & Intervalles de confiance (95%)
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~~Statut de sante gobale~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Estimateurs
QL2_b1 <- summary(fit.lme.QL2)$tTable[2,1]
QL2_b2 <- summary(fit.lme.QL2)$tTable[3,1]
pval_QL2_b1 <- summary(fit.lme.QL2)$tTable[2,5]
pval_QL2_b2 <- summary(fit.lme.QL2)$tTable[3,5]

#Beta1 + Beta2
QL2_b1_b2 <- QL2_b1 + QL2_b2

#IC
binf_QL2_b1 <- QL2_b1 - summary(fit.lme.QL2)$tTable[2,2]*1.96
bsup_QL2_b1 <- QL2_b1 + summary(fit.lme.QL2)$tTable[2,2]*1.96
binf_QL2_b2 <- QL2_b2 - summary(fit.lme.QL2)$tTable[3,2]*1.96
bsup_QL2_b2 <- QL2_b2 + summary(fit.lme.QL2)$tTable[3,2]*1.96

IC_QL2_b1 <- str_c("[", round(binf_QL2_b1, 3), ";", round(bsup_QL2_b1, 3), "]", 
                   sep = " ") 
IC_QL2_b2 <- str_c("[", round(binf_QL2_b2, 3), ";", round(bsup_QL2_b2, 3), "]", 
                   sep = " ")


#~~~~~~~~~~~~~~Echelles Fonctionnelles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#------Fonctions physiques-----------------------------------------------------
#Estimateurs
PF_b1 <- summary(fit.lme.PF)$tTable[2,1]
PF_b2 <- summary(fit.lme.PF)$tTable[3,1]
pval_PF_b1 <- summary(fit.lme.PF)$tTable[2,5]
pval_PF_b2 <- summary(fit.lme.PF)$tTable[3,5]

#Beta1 + Beta2
PF_b1_b2 <- PF_b1 + PF_b2

#IC
binf_PF_b1 <- PF_b1 - summary(fit.lme.PF)$tTable[2,2]*1.96
bsup_PF_b1 <- PF_b1 + summary(fit.lme.PF)$tTable[2,2]*1.96
binf_PF_b2 <- PF_b2 - summary(fit.lme.PF)$tTable[3,2]*1.96
bsup_PF_b2 <- PF_b2 + summary(fit.lme.PF)$tTable[3,2]*1.96

IC_PF_b1 <- str_c("[", round(binf_PF_b1, 3), ";", round(bsup_PF_b1, 3), "]", 
                  sep = " ") 
IC_PF_b2 <- str_c("[", round(binf_PF_b2, 3), ";", round(bsup_PF_b2, 3), "]", 
                  sep = " ")


#------Fonctions personnelles--------------------------------------------------
#Estimateurs
RF_b1 <- summary(fit.lme.RF)$tTable[2,1]
RF_b2 <- summary(fit.lme.RF)$tTable[3,1]
pval_RF_b1 <- summary(fit.lme.RF)$tTable[2,5]
pval_RF_b2 <- summary(fit.lme.RF)$tTable[3,5]

#Beta1 + Beta2
RF_b1_b2 <- RF_b1 + RF_b2

#IC
binf_RF_b1 <- RF_b1 - summary(fit.lme.RF)$tTable[2,2]*1.96
bsup_RF_b1 <- RF_b1 + summary(fit.lme.RF)$tTable[2,2]*1.96
binf_RF_b2 <- RF_b2 - summary(fit.lme.RF)$tTable[3,2]*1.96
bsup_RF_b2 <- RF_b2 + summary(fit.lme.RF)$tTable[3,2]*1.96

IC_RF_b1 <- str_c("[", round(binf_RF_b1, 3), ";", round(bsup_RF_b1, 3), "]", 
                  sep = " ") 
IC_RF_b2 <- str_c("[", round(binf_RF_b2, 3), ";", round(bsup_RF_b2, 3), "]", 
                  sep = " ")


#------Fonctions emotionnelles-------------------------------------------------
#Estimateurs
EF_b1 <- summary(fit.lme.EF)$tTable[2,1]
EF_b2 <- summary(fit.lme.EF)$tTable[3,1]
pval_EF_b1 <- summary(fit.lme.EF)$tTable[2,5]
pval_EF_b2 <- summary(fit.lme.EF)$tTable[3,5]

#Beta1 + Beta2
EF_b1_b2 <- EF_b1 + EF_b2

#IC
binf_EF_b1 <- EF_b1 - summary(fit.lme.EF)$tTable[2,2]*1.96
bsup_EF_b1 <- EF_b1 + summary(fit.lme.EF)$tTable[2,2]*1.96
binf_EF_b2 <- EF_b2 - summary(fit.lme.EF)$tTable[3,2]*1.96
bsup_EF_b2 <- EF_b2 + summary(fit.lme.EF)$tTable[3,2]*1.96

IC_EF_b1 <- str_c("[", round(binf_EF_b1, 3), ";", round(bsup_EF_b1, 3), "]", 
                  sep = " ") 
IC_EF_b2 <- str_c("[", round(binf_EF_b2, 3), ";", round(bsup_EF_b2, 3), "]", 
                  sep = " ")


#------Fonctions cognitives----------------------------------------------------
#Estimateurs
CF_b1 <- summary(fit.lme.CF)$tTable[2,1]
CF_b2 <- summary(fit.lme.CF)$tTable[3,1]
pval_CF_b1 <- summary(fit.lme.CF)$tTable[2,5]
pval_CF_b2 <- summary(fit.lme.CF)$tTable[3,5]

#Beta1 + Beta2
CF_b1_b2 <- CF_b1 + CF_b2

#IC
binf_CF_b1 <- CF_b1 - summary(fit.lme.CF)$tTable[2,2]*1.96
bsup_CF_b1 <- CF_b1 + summary(fit.lme.CF)$tTable[2,2]*1.96
binf_CF_b2 <- CF_b2 - summary(fit.lme.CF)$tTable[3,2]*1.96
bsup_CF_b2 <- CF_b2 + summary(fit.lme.CF)$tTable[3,2]*1.96

IC_CF_b1 <- str_c("[", round(binf_CF_b1, 3), ";", round(bsup_CF_b1, 3), "]", 
                  sep = " ") 
IC_CF_b2 <- str_c("[", round(binf_CF_b2, 3), ";", round(bsup_CF_b2, 3), "]", 
                  sep = " ")


#------Fonctions sociales------------------------------------------------------
#Estimateurs
SF_b1 <- summary(fit.lme.SF)$tTable[2,1]
SF_b2 <- summary(fit.lme.SF)$tTable[3,1]
pval_SF_b1 <- summary(fit.lme.SF)$tTable[2,5]
pval_SF_b2 <- summary(fit.lme.SF)$tTable[3,5]

#Beta1 + Beta2
SF_b1_b2 <- SF_b1 + SF_b2

#IC
binf_SF_b1 <- SF_b1 - summary(fit.lme.SF)$tTable[2,2]*1.96
bsup_SF_b1 <- SF_b1 + summary(fit.lme.SF)$tTable[2,2]*1.96
binf_SF_b2 <- SF_b2 - summary(fit.lme.SF)$tTable[3,2]*1.96
bsup_SF_b2 <- SF_b2 + summary(fit.lme.SF)$tTable[3,2]*1.96

IC_SF_b1 <- str_c("[", round(binf_SF_b1, 3), ";", round(bsup_SF_b1, 3), "]", 
                  sep = " ") 
IC_SF_b2 <- str_c("[", round(binf_SF_b2, 3), ";", round(bsup_SF_b2, 3), "]", 
                  sep = " ")


#~~~~~~~~~~~~~~Echelles Symptomatiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#------Fatigue-----------------------------------------------------------------
#Estimateurs
FA_b1 <- summary(fit.lme.FA)$tTable[2,1]
FA_b2 <- summary(fit.lme.FA)$tTable[3,1]
pval_FA_b1 <- summary(fit.lme.FA)$tTable[2,5]
pval_FA_b2 <- summary(fit.lme.FA)$tTable[3,5]

#Beta1 + Beta2
FA_b1_b2 <- FA_b1 + FA_b2

#IC
binf_FA_b1 <- FA_b1 - summary(fit.lme.FA)$tTable[2,2]*1.96
bsup_FA_b1 <- FA_b1 + summary(fit.lme.FA)$tTable[2,2]*1.96
binf_FA_b2 <- FA_b2 - summary(fit.lme.FA)$tTable[3,2]*1.96
bsup_FA_b2 <- FA_b2 + summary(fit.lme.FA)$tTable[3,2]*1.96

IC_FA_b1 <- str_c("[", round(binf_FA_b1, 3), ";", round(bsup_FA_b1, 3), "]", 
                  sep = " ") 
IC_FA_b2 <- str_c("[", round(binf_FA_b2, 3), ";", round(bsup_FA_b2, 3), "]", 
                  sep = " ")


#------Nausees et Vomissements-------------------------------------------------
#Estimateurs
NV_b1 <- summary(fit.lme.NV)$tTable[2,1]
NV_b2 <- summary(fit.lme.NV)$tTable[3,1]
pval_NV_b1 <- summary(fit.lme.NV)$tTable[2,5]
pval_NV_b2 <- summary(fit.lme.NV)$tTable[3,5]

#Beta1 + Beta2
NV_b1_b2 <- NV_b1 + NV_b2

#IC
binf_NV_b1 <- NV_b1 - summary(fit.lme.NV)$tTable[2,2]*1.96
bsup_NV_b1 <- NV_b1 + summary(fit.lme.NV)$tTable[2,2]*1.96
binf_NV_b2 <- NV_b2 - summary(fit.lme.NV)$tTable[3,2]*1.96
bsup_NV_b2 <- NV_b2 + summary(fit.lme.NV)$tTable[3,2]*1.96

IC_NV_b1 <- str_c("[", round(binf_NV_b1, 3), ";", round(bsup_NV_b1, 3), "]", 
                  sep = " ") 
IC_NV_b2 <- str_c("[", round(binf_NV_b2, 3), ";", round(bsup_NV_b2, 3), "]", 
                  sep = " ")


#------Douleurs----------------------------------------------------------------
#Estimateurs
PA_b1 <- summary(fit.lme.PA)$tTable[2,1]
PA_b2 <- summary(fit.lme.PA)$tTable[3,1]
pval_PA_b1 <- summary(fit.lme.PA)$tTable[2,5]
pval_PA_b2 <- summary(fit.lme.PA)$tTable[3,5]

#Beta1 + Beta2
PA_b1_b2 <- PA_b1 + PA_b2

#IC
binf_PA_b1 <- PA_b1 - summary(fit.lme.PA)$tTable[2,2]*1.96
bsup_PA_b1 <- PA_b1 + summary(fit.lme.PA)$tTable[2,2]*1.96
binf_PA_b2 <- PA_b2 - summary(fit.lme.PA)$tTable[3,2]*1.96
bsup_PA_b2 <- PA_b2 + summary(fit.lme.PA)$tTable[3,2]*1.96

IC_PA_b1 <- str_c("[", round(binf_PA_b1, 3), ";", round(bsup_PA_b1, 3), "]", 
                  sep = " ") 
IC_PA_b2 <- str_c("[", round(binf_PA_b2, 3), ";", round(bsup_PA_b2, 3), "]", 
                  sep = " ")


#------Dyspnee-----------------------------------------------------------------
#Estimateurs
DY_b1 <- summary(fit.lme.DY)$tTable[2,1]
DY_b2 <- summary(fit.lme.DY)$tTable[3,1]
pval_DY_b1 <- summary(fit.lme.DY)$tTable[2,5]
pval_DY_b2 <- summary(fit.lme.DY)$tTable[3,5]

#Beta1 + Beta2
DY_b1_b2 <- DY_b1 + DY_b2

#IC
binf_DY_b1 <- DY_b1 - summary(fit.lme.DY)$tTable[2,2]*1.96
bsup_DY_b1 <- DY_b1 + summary(fit.lme.DY)$tTable[2,2]*1.96
binf_DY_b2 <- DY_b2 - summary(fit.lme.DY)$tTable[3,2]*1.96
bsup_DY_b2 <- DY_b2 + summary(fit.lme.DY)$tTable[3,2]*1.96

IC_DY_b1 <- str_c("[", round(binf_DY_b1, 3), ";", round(bsup_DY_b1, 3), "]", 
                  sep = " ") 
IC_DY_b2 <- str_c("[", round(binf_DY_b2, 3), ";", round(bsup_DY_b2, 3), "]", 
                  sep = " ")


#------Insomnies---------------------------------------------------------------
#Estimateurs
SL_b1 <- summary(fit.lme.SL)$tTable[2,1]
SL_b2 <- summary(fit.lme.SL)$tTable[3,1]
pval_SL_b1 <- summary(fit.lme.SL)$tTable[2,5] 
pval_SL_b2 <- summary(fit.lme.SL)$tTable[3,5] 

#Beta1 + Beta2
SL_b1_b2 <- SL_b1 + SL_b2

#IC
binf_SL_b1 <- SL_b1 - summary(fit.lme.SL)$tTable[2,2]*1.96 
bsup_SL_b1 <- SL_b1 + summary(fit.lme.SL)$tTable[2,2]*1.96 
binf_SL_b2 <- SL_b2 - summary(fit.lme.SL)$tTable[3,2]*1.96 
bsup_SL_b2 <- SL_b2 + summary(fit.lme.SL)$tTable[3,2]*1.96 

IC_SL_b1 <- str_c("[", round(binf_SL_b1, 3), ";", round(bsup_SL_b1, 3), "]", 
                  sep = " ") 
IC_SL_b2 <- str_c("[", round(binf_SL_b2, 3), ";", round(bsup_SL_b2, 3), "]", 
                  sep = " ")


#------Pertes d'appetit--------------------------------------------------------
#Estimateurs
AP_b1 <- summary(fit.lme.AP)$tTable[2,1] 
AP_b2 <- summary(fit.lme.AP)$tTable[3,1] 
pval_AP_b1 <- summary(fit.lme.AP)$tTable[2,5] 
pval_AP_b2 <- summary(fit.lme.AP)$tTable[3,5] 

#Beta1 + Beta2
AP_b1_b2 <- AP_b1 + AP_b2

#IC
binf_AP_b1 <- AP_b1 - summary(fit.lme.AP)$tTable[2,2]*1.96 
bsup_AP_b1 <- AP_b1 + summary(fit.lme.AP)$tTable[2,2]*1.96 
binf_AP_b2 <- AP_b2 - summary(fit.lme.AP)$tTable[3,2]*1.96 
bsup_AP_b2 <- AP_b2 + summary(fit.lme.AP)$tTable[3,2]*1.96 

IC_AP_b1 <- str_c("[", round(binf_AP_b1, 3), ";", round(bsup_AP_b1, 3), "]", 
                  sep = " ") 
IC_AP_b2 <- str_c("[", round(binf_AP_b2, 3), ";", round(bsup_AP_b2, 3), "]", 
                  sep = " ")


#------Constipations--------------------------------------------------
#Estimateurs
CO_b1 <- summary(fit.lme.CO)$tTable[2,1] 
CO_b2 <- summary(fit.lme.CO)$tTable[3,1] 
pval_CO_b1 <- summary(fit.lme.CO)$tTable[2,5] 
pval_CO_b2 <- summary(fit.lme.CO)$tTable[3,5] 

#Beta1 + Beta2
CO_b1_b2 <- CO_b1 + CO_b2

#IC
binf_CO_b1 <- CO_b1 - summary(fit.lme.CO)$tTable[2,2]*1.96 
bsup_CO_b1 <- CO_b1 + summary(fit.lme.CO)$tTable[2,2]*1.96 
binf_CO_b2 <- CO_b2 - summary(fit.lme.CO)$tTable[3,2]*1.96 
bsup_CO_b2 <- CO_b2 + summary(fit.lme.CO)$tTable[3,2]*1.96 

IC_CO_b1 <- str_c("[", round(binf_CO_b1, 3), ";", round(bsup_CO_b1, 3), "]", 
                  sep = " ") 
IC_CO_b2 <- str_c("[", round(binf_CO_b2, 3), ";", round(bsup_CO_b2, 3), "]", 
                  sep = " ")


#------Diarrhees--------------------------------------------------
#Estimateurs
DI_b1 <- summary(fit.lme.DI)$tTable[2,1] 
DI_b2 <- summary(fit.lme.DI)$tTable[3,1] 
pval_DI_b1 <- summary(fit.lme.DI)$tTable[2,5] 
pval_DI_b2 <- summary(fit.lme.DI)$tTable[3,5] 

#Beta1 + Beta2
DI_b1_b2 <- DI_b1 + DI_b2

#IC
binf_DI_b1 <- DI_b1 - summary(fit.lme.DI)$tTable[2,2]*1.96 
bsup_DI_b1 <- DI_b1 + summary(fit.lme.DI)$tTable[2,2]*1.96 
binf_DI_b2 <- DI_b2 - summary(fit.lme.DI)$tTable[3,2]*1.96 
bsup_DI_b2 <- DI_b2 + summary(fit.lme.DI)$tTable[3,2]*1.96 

IC_DI_b1 <- str_c("[", round(binf_DI_b1, 3), ";", round(bsup_DI_b1, 3), "]", 
                  sep = " ") 
IC_DI_b2 <- str_c("[", round(binf_DI_b2, 3), ";", round(bsup_DI_b2, 3), "]", 
                  sep = " ")


#------Dificultees financieres--------------------------------------------------
#Estimateurs
FI_b1 <- summary(fit.lme.FI)$tTable[2,1] 
FI_b2 <- summary(fit.lme.FI)$tTable[3,1] 
pval_FI_b1 <- summary(fit.lme.FI)$tTable[2,5] 
pval_FI_b2 <- summary(fit.lme.FI)$tTable[3,5] 

#Beta1 + Beta2
FI_b1_b2 <- FI_b1 + FI_b2

#IC
binf_FI_b1 <- FI_b1 - summary(fit.lme.FI)$tTable[2,2]*1.96 
bsup_FI_b1 <- FI_b1 + summary(fit.lme.FI)$tTable[2,2]*1.96 
binf_FI_b2 <- FI_b2 - summary(fit.lme.FI)$tTable[3,2]*1.96 
bsup_FI_b2 <- FI_b2 + summary(fit.lme.FI)$tTable[3,2]*1.96 

IC_FI_b1 <- str_c("[", round(binf_FI_b1, 3), ";", round(bsup_FI_b1, 3), "]", 
                  sep = " ") 
IC_FI_b2 <- str_c("[", round(binf_FI_b2, 3), ";", round(bsup_FI_b2, 3), "]", 
                  sep = " ")

#~~~~~~~~~~~~~Tableau Final~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TAB <- data.frame(matrix("", 17, 8))
TAB <- data.frame(c("Statut Global Sante", "ECHELLES FONCTIONNELLES", 
                    "Physique", "Personnelles", "Emotionnelles", "Cognitives",
                    "Sociales", "ECHELLES SYMPTOMATIQUES", "Fatigue", 
                    "Nausees et vomissements", "Douleurs", "Dyspnee",
                    "Insomnies", "Pertes d'appetit", "Constipations", 
                    "Diarrhees", "Difficultees Financieres"), 
                  c(round(QL2_b1, 3), "-----", round(PF_b1, 3), round(RF_b1, 3),
                    round(EF_b1, 3), round(CF_b1, 3), round(SF_b1, 3), "-----",
                    round(FA_b1, 3), round(NV_b1, 3), round(PA_b1, 3), 
                    round(DY_b1, 3), round(SL_b1, 3), round(AP_b1, 3), 
                    round(CO_b1, 3), round(DI_b1, 3), round(FI_b1, 3)), 
                  c(IC_QL2_b1, "-----", IC_PF_b1, 
                    IC_RF_b1, IC_EF_b1, IC_CF_b1,
                    IC_SF_b1, "-----", IC_FA_b1, 
                    IC_NV_b1, IC_PA_b1, IC_DY_b1,
                    IC_SL_b1, IC_AP_b1, IC_CO_b1,
                    IC_DI_b1, IC_FI_b1), 
                  c(round(pval_QL2_b1, 3), "-----", round(pval_PF_b1, 3), 
                    round(pval_RF_b1, 3), round(pval_EF_b1, 3), 
                    round(pval_CF_b1, 3), round(pval_SF_b1, 3), "-----",
                    round(pval_FA_b1, 3), round(pval_NV_b1, 3), 
                    round(pval_PA_b1, 3), round(pval_DY_b1, 3),
                    round(pval_SL_b1, 3), round(pval_AP_b1, 3), 
                    round(pval_CO_b1, 3), round(pval_DI_b1, 3),
                    round(pval_FI_b1, 3)),
                  c(round(QL2_b1_b2, 3), "-----", round(PF_b1_b2, 3), 
                    round(RF_b1_b2, 3), round(EF_b1_b2, 3), round(CF_b1_b2, 3), 
                    round(SF_b1_b2, 3), "-----", round(FA_b1_b2, 3), 
                    round(NV_b1_b2, 3), round(PA_b1_b2, 3), round(DY_b1_b2, 3), 
                    round(SL_b1_b2, 3), round(AP_b1_b2, 3), round(CO_b1_b2, 3),
                    round(DI_b1_b2, 3), round(FI_b1_b2, 3)), 
                  c(round(QL2_b2, 3), "-----", round(PF_b2, 3), round(RF_b2, 3),
                    round(EF_b2, 3), round(CF_b2, 3), round(SF_b2, 3), "-----",
                    round(FA_b2, 3), round(NV_b2, 3), round(PA_b2, 3), 
                    round(DY_b2, 3), round(SL_b2, 3), round(AP_b2, 3), 
                    round(CO_b2, 3), round(DI_b2, 3), round(FI_b2, 3)), 
                  c(IC_QL2_b2, "-----", IC_PF_b2, 
                    IC_RF_b2, IC_EF_b2, IC_CF_b2,
                    IC_SF_b2, "-----", IC_FA_b2, 
                    IC_NV_b2, IC_PA_b2, IC_DY_b2,
                    IC_SL_b2, IC_AP_b2, IC_CO_b2, 
                    IC_DI_b2, IC_FI_b2),
                  c(round(pval_QL2_b2, 3), "-----", round(pval_PF_b2, 3), 
                    round(pval_RF_b2, 3), round(pval_EF_b2, 3), 
                    round(pval_CF_b2, 3), round(pval_SF_b2, 3), "-----", 
                    round(pval_FA_b2, 3), round(pval_NV_b2, 3),
                    round(pval_PA_b2, 3), round(pval_DY_b2, 3), 
                    round(pval_SL_b2, 3), round(pval_AP_b2, 3), 
                    round(pval_CO_b2, 3), round(pval_DI_b2, 3), 
                    round(pval_FI_b2, 3)))
names(TAB) <- c("", "Beta1_chap", "95% CI", "p-value", "Beta1 + Beta2", 
                "Beta2_chap", "95% CI", "p_value")
View(TAB)
# write.table(TAB, file = "tableau_estimations_mois_REML.csv",
#             sep = ";", dec = ",", row.names = FALSE)


###############################################################################
#                        GRAPHIQUES DES PREDICTIONS                           #
###############################################################################

#------------------------------------------------------------------------------
#                       Statut de sante gobale
#------------------------------------------------------------------------------

#~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.QL2$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données PF2 manquantes:
df_QL2_lmm <- df_tpreel[which(!is.na(df_tpreel$QL2)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_QL2_lmm <-  subset(df_QL2_lmm, select =-c(visit:LagDate, PF2:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_QL2_lmm )
row.names(df_QL2_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_QL2_lmm <- cbind(df_QL2_lmm, fit.lme.QL2$fitted[,1])
names(df_QL2_lmm)[ncol(df_QL2_lmm)] <- "QL2_lmm.fitted"
# View(df_QL2_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne QL2:
df_QL2_spl <- df_tpreel[which(!is.na(df_tpreel$QL2)),]
#ne garder que les colonnes qui sont utiles pour l'échelle
df_QL2_spl <-  subset(df_QL2_spl, select =-c(visit:LagDate, PF2:FI))
df_QL2_spl <- df_QL2_spl[order(df_QL2_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_QL2_spl)
row.names(df_QL2_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.QL2.1 <- lme(QL2 ~ ns(Time,3), 
                        random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                        data=df_QL2_spl, na.action = na.omit,
                        control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.QL2.2 <- lme(QL2 ~ ns(Time,3):BRAS, 
                        random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                        data=df_QL2_spl, na.action = na.omit,
                        control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.QL2.1, fit.lme.ns.QL2.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.QL2.1.ml <- update(fit.lme.ns.QL2.1, method="ML")
fit.lme.ns.QL2.2.ml <- update(fit.lme.ns.QL2.2, method="ML")
QL2_test <- lrtest(fit.lme.ns.QL2.1.ml, fit.lme.ns.QL2.2.ml)
#p-value < chi = effet significatif du bras

QL2_pval <- QL2_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_QL2_spl <- cbind(df_QL2_spl, fit.lme.ns.QL2.2$fitted[,1])
names(df_QL2_spl)[ncol(df_QL2_spl)] <- "QL2_spl.fitted"
# View(df_QL2_spl)  -> derniere colonne ajoutee

#~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
QL2_ <- ggplot(df_QL2_lmm, aes(x=Time, y=QL2_lmm.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(40,80))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Statut Global de Santé (trajectoires linéaires)")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

#SPLINES
QL2__ <- ggplot(df_QL2_spl, aes(x=Time, y=QL2_spl.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(40,80))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color = "Bras") + 
  ggtitle("Statut Global de Santé (trajectoires splines)") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

#sorties graphiques
grid.arrange(QL2_, QL2__)

#------------------------------------------------------------------------------
#                       Fonctions physiques
#------------------------------------------------------------------------------

#~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.PF$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données PF2 manquantes:
df_PF_lmm <- df_tpreel[which(!is.na(df_tpreel$PF2)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_PF_lmm <-  subset(df_PF_lmm, select =-c(visit:LagDate, QL2, RF2:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_PF_lmm )
row.names(df_PF_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_PF_lmm <- cbind(df_PF_lmm, fit.lme.PF$fitted[,1])
names(df_PF_lmm)[ncol(df_PF_lmm)] <- "PF_lmm.fitted"
# View(df_PF_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne PF2
df_PF_spl <- df_tpreel[which(!is.na(df_tpreel$PF2)),]
#ne garder que les colonnes qui sont utiles pour l'échelle
df_PF_spl <-  subset(df_PF_spl, select =-c(visit:LagDate, QL2, RF2:FI))
df_PF_spl <- df_PF_spl[order(df_PF_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_PF_spl)
row.names(df_PF_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.PF.1 <- lme(PF2 ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_PF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.PF.2 <- lme(PF2 ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_PF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.PF.1, fit.lme.ns.PF.2)
#p-value tres petite = effet significatif du bras

#Tester les modeles
fit.lme.ns.PF.1.ml <- update(fit.lme.ns.PF.1, method="ML")
fit.lme.ns.PF.2.ml <- update(fit.lme.ns.PF.2, method="ML")
PF_test <- lrtest(fit.lme.ns.PF.1.ml, fit.lme.ns.PF.2.ml)
#p-value < chi = effet significatif du bras

PF_pval <- PF_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_PF_spl <- cbind(df_PF_spl, fit.lme.ns.PF.2$fitted[,1])
names(df_PF_spl)[ncol(df_PF_spl)] <- "PF_spl.fitted"
# View(df_PF_spl)  -> derniere colonne ajoutee

#~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRE
PF_ <- ggplot(df_PF_lmm, aes(x=Time, y=PF_lmm.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(65,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction Physique (trajectoires linéaires)")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

#SPLINES
PF__ <- ggplot(df_PF_spl, aes(x=Time, y=PF_spl.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(65,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color = "Bras") + 
  ggtitle("Fonction Physique (trajectoires splines)") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

grid.arrange(PF_, PF__)

#------------------------------------------------------------------------------
#                       Fonctions personnelles
#------------------------------------------------------------------------------

#~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.RF$fitted[,1]) # #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees RF2 manquantes:
df_RF_lmm <- df_tpreel[which(!is.na(df_tpreel$RF2)),] # retire (si il y en a)
#on garde que les colonnes utilies pour l'échelle:
df_RF_lmm <-  subset(df_RF_lmm, select =-c(visit:LagDate, QL2:PF2, EF:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_RF_lmm )
row.names(df_RF_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_RF_lmm <- cbind(df_RF_lmm, fit.lme.RF$fitted[,1])
names(df_RF_lmm)[ncol(df_RF_lmm)] <- "RF_lmm.fitted"
# View(df_RF_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne RF:
df_RF_spl <- df_tpreel[which(!is.na(df_tpreel$RF2)),]
#on garde que les colonnes utiles pour l'échelle:
df_RF_spl <-  subset(df_RF_spl, select =-c(visit:LagDate, QL2:PF2, EF:FI))
df_RF_spl <- df_RF_spl[order(df_RF_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_RF_spl)
row.names(df_RF_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.RF.1 <- lme(RF2 ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))), 
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_RF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.RF.2 <- lme(RF2 ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_RF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.RF.1, fit.lme.ns.RF.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.RF.1.ml <- update(fit.lme.ns.RF.1, method="ML")
fit.lme.ns.RF.2.ml <- update(fit.lme.ns.RF.2, method="ML")
RF_test <- lrtest(fit.lme.ns.RF.1.ml, fit.lme.ns.RF.2.ml) 
RF_pval <- RF_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_RF_spl <- cbind(df_RF_spl, fit.lme.ns.RF.2$fitted[,1])
names(df_RF_spl)[ncol(df_RF_spl)] <- "RF_spl.fitted"
# View(df_RF_spl)  -> derniere colonne ajoutee

#~~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
RF_ <- ggplot(df_RF_lmm, aes(x=Time, y=RF_lmm.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(65,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction Personnelle (trajectoires linéaires)")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

#SPLINES
RF__ <- ggplot(df_RF_spl, aes(x=Time, y=RF_spl.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(65,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color = "Bras") + 
  ggtitle("Fonction Personnelle (trajectoires splines)") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

grid.arrange(RF_, RF__)

#------------------------------------------------------------------------------
#                       Fonctions sociale
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.SF$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees SF manquantes:
df_SF_lmm <- df_tpreel[which(!is.na(df_tpreel$SF)),] # retire (si il y en a)
#on garde les colonnes utiles pour l'échelle:
df_SF_lmm <-  subset(df_SF_lmm, select =-c(visit:LagDate, QL2:CF, FA:FI))
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_SF_lmm )
row.names(df_SF_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_SF_lmm <- cbind(df_SF_lmm, fit.lme.SF$fitted[,1])
names(df_SF_lmm)[ncol(df_SF_lmm)] <- "SF_lmm.fitted"
# View(df_SF_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne SF
df_SF_spl <- df_tpreel[which(!is.na(df_tpreel$SF)),]
#on garde les colonnes utiles pour l'échelle:
df_SF_spl <-  subset(df_SF_spl, select =-c(visit:LagDate, QL2:CF, FA:FI))
df_SF_spl <- df_SF_spl[order(df_SF_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_SF_spl)
row.names(df_SF_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.SF.1 <- lme(SF ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_SF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.SF.2 <- lme(SF ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_SF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.SF.1, fit.lme.ns.SF.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.SF.1.ml <- update(fit.lme.ns.SF.1, method="ML")
fit.lme.ns.SF.2.ml <- update(fit.lme.ns.SF.2, method="ML")
SF_test <- lrtest(fit.lme.ns.SF.1.ml, fit.lme.ns.SF.2.ml) 
SF_pval <- SF_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_SF_spl <- cbind(df_SF_spl, fit.lme.ns.SF.2$fitted[,1])
names(df_SF_spl)[ncol(df_SF_spl)] <- "SF_spl.fitted"
# View(df_SF_spl)  -> derniere colonne ajoutee

#~~~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
SF_ <- ggplot(df_SF_lmm, aes(x=Time, y=SF_lmm.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(60,95))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction Sociale (trajectoires linéaires)")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

#SPLINES
SF__ <- ggplot(df_SF_spl, aes(x=Time, y=SF_spl.fitted, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(60,95))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color = "Bras") + 
  ggtitle("Fonction Sociale (trajectoires splines)") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))

grid.arrange(SF_, SF__)

#------------------------------------------------------------------------------
#                       Fonctions emotionnelles
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~LMMs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.EF$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte 
# on construit un fichier ou on supprime les lignes avec donnees EF manquantes:
df_EF_lmm <- df_tpreel[which(!is.na(df_tpreel$EF)),] # retire (si il y en a)
#on ne garde que les colonnes utiles pour l'échelle:
df_EF_lmm <-  subset(df_EF_lmm, select =-c(visit:LagDate, QL2:RF2, CF:FI))
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_EF_lmm )
row.names(df_EF_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_EF_lmm <- cbind(df_EF_lmm, fit.lme.EF$fitted[,1])
names(df_EF_lmm)[ncol(df_EF_lmm)] <- "EF_lmm.fitted"
# View(df_EF_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne EF:
df_EF_spl <- df_tpreel[which(!is.na(df_tpreel$EF)),]
#on ne garde que les colonnes utiles à l'échelle:
df_EF_spl <-  subset(df_EF_spl, select =-c(visit:LagDate, QL2:RF2, CF:FI))
df_EF_spl <- df_EF_spl[order(df_EF_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_EF_spl)
row.names(df_EF_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.EF.1 <- lme(EF ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_EF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.EF.2 <- lme(EF ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_EF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.EF.1, fit.lme.ns.EF.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.EF.1.ml <- update(fit.lme.ns.EF.1, method="ML")
fit.lme.ns.EF.2.ml <- update(fit.lme.ns.EF.2, method="ML")
EF_test <- lrtest(fit.lme.ns.EF.1.ml, fit.lme.ns.EF.2.ml) 
EF_pval <- EF_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_EF_spl <- cbind(df_EF_spl, fit.lme.ns.EF.2$fitted[,1])
names(df_EF_spl)[ncol(df_EF_spl)] <- "EF_spl.fitted"
# View(df_EF_spl)  -> derniere colonne ajoutee

#~~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_EF_lmm$Time[which(df_EF_lmm$BRAS==0)],
     df_EF_lmm$EF_lmm.fitted[which(df_EF_lmm$BRAS==0)], 
     ylim = c(60,90), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", 
     main = "Fonctions Emotionnelles (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_EF_lmm$Time[which(df_EF_lmm$BRAS==1)],
      df_EF_lmm$EF_lmm.fitted[which(df_EF_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 90, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_EF_spl$Time[which(df_EF_spl$BRAS==0)], 
     df_EF_spl$EF_spl.fitted[which(df_EF_spl$BRAS==0)], 
     col=c("blue"), ylim=c(60,90), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Fonctions Emotionnelles (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_EF_spl$Time[which(df_EF_spl$BRAS==1)], 
      df_EF_spl$EF_spl.fitted[which(df_EF_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 90, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                       Fonctions cognitives
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.CF$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte 
#on construit un fichier ou on supprime les lignes avec donnees CF manquantes:
df_CF_lmm <- df_tpreel[which(!is.na(df_tpreel$CF)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_CF_lmm <-  subset(df_CF_lmm, select =-c(visit:LagDate, QL2:EF, SF:FI))
#vérfication des dimensions des tableaux
dim(df_tpreel)
dim(df_CF_lmm )
row.names(df_CF_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_CF_lmm <- cbind(df_CF_lmm, fit.lme.CF$fitted[,1])
names(df_CF_lmm)[ncol(df_CF_lmm)] <- "CF_lmm.fitted"
# View(df_CF_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne CF:
df_CF_spl <- df_tpreel[which(!is.na(df_tpreel$CF)),]
#on ne garde que les colonnes utiles à l'échelle:
df_CF_spl <-  subset(df_CF_spl, select =-c(visit:LagDate, QL2:EF, SF:FI))
df_CF_spl <- df_CF_spl[order(df_CF_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_CF_spl)
row.names(df_CF_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.CF.1 <- lme(CF ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_CF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.CF.2 <- lme(CF ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_CF_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.CF.1, fit.lme.ns.CF.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.CF.1.ml <- update(fit.lme.ns.CF.1, method="ML")
fit.lme.ns.CF.2.ml <- update(fit.lme.ns.CF.2, method="ML")
CF_test <- lrtest(fit.lme.ns.CF.1.ml, fit.lme.ns.CF.2.ml) 
CF_pval <- CF_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_CF_spl <- cbind(df_CF_spl, fit.lme.ns.CF.2$fitted[,1])
names(df_CF_spl)[ncol(df_CF_spl)] <- "CF_spl.fitted"
# View(df_CF_spl)  -> derniere colonne ajoutee

#~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_CF_lmm$Time[which(df_CF_lmm$BRAS==0)],
     df_CF_lmm$CF_lmm.fitted[which(df_CF_lmm$BRAS==0)], 
     ylim = c(65,90), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", 
     main = "Fonctions Cognitives (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_CF_lmm$Time[which(df_CF_lmm$BRAS==1)],
      df_CF_lmm$CF_lmm.fitted[which(df_CF_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 70, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_CF_spl$Time[which(df_CF_spl$BRAS==0)], 
     df_CF_spl$CF_spl.fitted[which(df_CF_spl$BRAS==0)], 
     col=c("blue"), ylim=c(65,90), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Fonctions Cognitives (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_CF_spl$Time[which(df_CF_spl$BRAS==1)], 
      df_CF_spl$CF_spl.fitted[which(df_CF_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 70, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)


#------------------------------------------------------------------------------
#                             Fatigue
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.FA$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte 
#on construit un fichier ou on supprime les lignes avec donnees FA manquantes:
df_FA_lmm <- df_tpreel[which(!is.na(df_tpreel$FA)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_FA_lmm <-  subset(df_FA_lmm, select =-c(visit:LagDate, QL2:SF, NV:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_FA_lmm )
row.names(df_FA_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_FA_lmm <- cbind(df_FA_lmm, fit.lme.FA$fitted[,1])
names(df_FA_lmm)[ncol(df_FA_lmm)] <- "FA_lmm.fitted"
# View(df_FA_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne FA:
df_FA_spl <- df_tpreel[which(!is.na(df_tpreel$FA)),]
#on ne garde que les colonnes utiles à l'échelle:
df_FA_spl <-  subset(df_FA_spl, select =-c(visit:LagDate, QL2:SF, NV:FI))
df_FA_spl <- df_FA_spl[order(df_FA_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_FA_spl)
row.names(df_FA_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.FA.1 <- lme(FA ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_FA_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.FA.2 <- lme(FA ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_FA_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.FA.1, fit.lme.ns.FA.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.FA.1.ml <- update(fit.lme.ns.FA.1, method="ML")
fit.lme.ns.FA.2.ml <- update(fit.lme.ns.FA.2, method="ML")
FA_test <- lrtest(fit.lme.ns.FA.1.ml, fit.lme.ns.FA.2.ml) 
FA_pval <- FA_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_FA_spl <- cbind(df_FA_spl, fit.lme.ns.FA.2$fitted[,1])
names(df_FA_spl)[ncol(df_FA_spl)] <- "FA_spl.fitted"
# View(df_FA_spl)  -> derniere colonne ajoutee

#~~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_FA_lmm$Time[which(df_FA_lmm$BRAS==0)],
     df_FA_lmm$FA_lmm.fitted[which(df_FA_lmm$BRAS==0)], 
     ylim = c(20,50), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Fatigue (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_FA_lmm$Time[which(df_FA_lmm$BRAS==1)],
      df_FA_lmm$FA_lmm.fitted[which(df_FA_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 27, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_FA_spl$Time[which(df_FA_spl$BRAS==0)], 
     df_FA_spl$FA_spl.fitted[which(df_FA_spl$BRAS==0)], 
     col=c("blue"), ylim=c(20,50), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Fatigue (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_FA_spl$Time[which(df_FA_spl$BRAS==1)], 
      df_FA_spl$FA_spl.fitted[which(df_FA_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 27, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                    Nausees et Vomissements
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.NV$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte 
#on construit un fichier ou on supprime les lignes avec donnees NV manquantes:
df_NV_lmm <- df_tpreel[which(!is.na(df_tpreel$NV)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_NV_lmm <-  subset(df_NV_lmm, select =-c(visit:LagDate, QL2:FA, PA:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_NV_lmm )
row.names(df_NV_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_NV_lmm <- cbind(df_NV_lmm, fit.lme.NV$fitted[,1])
names(df_NV_lmm)[ncol(df_NV_lmm)] <- "NV_lmm.fitted"
# View(df_NV_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne NV:
df_NV_spl <- df_tpreel[which(!is.na(df_tpreel$NV)),]
#on ne garde que les colonnes utiles à l'échelle:
df_NV_spl <-  subset(df_NV_spl, select =-c(visit:LagDate, QL2:FA, PA:FI))
df_NV_spl <- df_NV_spl[order(df_NV_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_NV_spl)
row.names(df_NV_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.NV.1 <- lme(NV ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_NV_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.NV.2 <- lme(NV ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_NV_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.NV.1, fit.lme.ns.NV.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.NV.1.ml <- update(fit.lme.ns.NV.1, method="ML")
fit.lme.ns.NV.2.ml <- update(fit.lme.ns.NV.2, method="ML")
NV_test <- lrtest(fit.lme.ns.NV.1.ml, fit.lme.ns.NV.2.ml) 
NV_pval <- NV_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_NV_spl <- cbind(df_NV_spl, fit.lme.ns.NV.2$fitted[,1])
names(df_NV_spl)[ncol(df_NV_spl)] <- "NV_spl.fitted"
# View(df_NV_spl)  -> derniere colonne ajoutee

#~~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_NV_lmm$Time[which(df_NV_lmm$BRAS==0)],
     df_NV_lmm$NV_lmm.fitted[which(df_NV_lmm$BRAS==0)], 
     ylim = c(0,10), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Nausees et Vomissements (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_NV_lmm$Time[which(df_NV_lmm$BRAS==1)],
      df_NV_lmm$NV_lmm.fitted[which(df_NV_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 2, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_NV_spl$Time[which(df_NV_spl$BRAS==0)], 
     df_NV_spl$NV_spl.fitted[which(df_NV_spl$BRAS==0)], 
     col=c("blue"), ylim=c(0,10), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Nausees et Vomissements (trajectoires splines)", 
     cex.main = 1, type = "l")
lines(df_NV_spl$Time[which(df_NV_spl$BRAS==1)], 
      df_NV_spl$NV_spl.fitted[which(df_NV_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 2, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                                Douleurs
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.PA$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte 
#on construit un fichier ou on supprime les lignes avec donnees PA manquantes :
df_PA_lmm <- df_tpreel[which(!is.na(df_tpreel$PA)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_PA_lmm <-  subset(df_PA_lmm, select =-c(visit:LagDate, QL2:NV, DY:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_PA_lmm )
row.names(df_PA_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_PA_lmm <- cbind(df_PA_lmm, fit.lme.PA$fitted[,1])
names(df_PA_lmm)[ncol(df_PA_lmm)] <- "PA_lmm.fitted"
# View(df_PA_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne PA:
df_PA_spl <- df_tpreel[which(!is.na(df_tpreel$PA)),]
#on en garde que les colonnes utiles à l'échelle:
df_PA_spl <-  subset(df_PA_spl, select =-c(visit:LagDate, QL2:NV, DY:FI))
df_PA_spl <- df_PA_spl[order(df_PA_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_PA_spl)
row.names(df_PA_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.PA.1 <- lme(PA ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_PA_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.PA.2 <- lme(PA ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_PA_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.PA.1, fit.lme.ns.PA.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.PA.1.ml <- update(fit.lme.ns.PA.1, method="ML")
fit.lme.ns.PA.2.ml <- update(fit.lme.ns.PA.2, method="ML")
PA_test <- lrtest(fit.lme.ns.PA.1.ml, fit.lme.ns.PA.2.ml) 
PA_pval <- PA_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_PA_spl <- cbind(df_PA_spl, fit.lme.ns.PA.2$fitted[,1])
names(df_PA_spl)[ncol(df_PA_spl)] <- "PA_spl.fitted"
# View(df_PA_spl)  -> derniere colonne ajoutee

#~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_PA_lmm$Time[which(df_PA_lmm$BRAS==0)],
     df_PA_lmm$PA_lmm.fitted[which(df_PA_lmm$BRAS==0)], 
     ylim = c(15,50), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Douleur (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_PA_lmm$Time[which(df_PA_lmm$BRAS==1)],
      df_PA_lmm$PA_lmm.fitted[which(df_PA_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 50, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_PA_spl$Time[which(df_PA_spl$BRAS==0)], 
     df_PA_spl$PA_spl.fitted[which(df_PA_spl$BRAS==0)], 
     col=c("blue"), ylim=c(15,50), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Douleur (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_PA_spl$Time[which(df_PA_spl$BRAS==1)], 
      df_PA_spl$PA_spl.fitted[which(df_PA_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 50, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                       Dyspnee
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.DY$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees DY manquantes:
df_DY_lmm <- df_tpreel[which(!is.na(df_tpreel$DY)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_DY_lmm <-  subset(df_DY_lmm, select =-c(visit:LagDate, QL2:PA, SL:FI))
#vérfication des dimensions des tableaux:
dim(df_tpreel)
dim(df_DY_lmm )
row.names(df_DY_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_DY_lmm <- cbind(df_DY_lmm, fit.lme.DY$fitted[,1])
names(df_DY_lmm)[ncol(df_DY_lmm)] <- "DY_lmm.fitted"
# View(df_DY_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne DY:
df_DY_spl <- df_tpreel[which(!is.na(df_tpreel$DY)),]
#on ne garde les colonnes utiles à l'échelle:
df_DY_spl <-  subset(df_DY_spl, select =-c(visit:LagDate, QL2:PA, SL:FI))
df_DY_spl <- df_DY_spl[order(df_DY_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_DY_spl)
row.names(df_DY_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.DY.1 <- lme(DY ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_DY_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.DY.2 <- lme(DY ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_DY_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.DY.1, fit.lme.ns.DY.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.DY.1.ml <- update(fit.lme.ns.DY.1, method="ML")
fit.lme.ns.DY.2.ml <- update(fit.lme.ns.DY.2, method="ML")
DY_test <- lrtest(fit.lme.ns.DY.1.ml, fit.lme.ns.DY.2.ml) 
DY_pval <- DY_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_DY_spl <- cbind(df_DY_spl, fit.lme.ns.DY.2$fitted[,1])
names(df_DY_spl)[ncol(df_DY_spl)] <- "DY_spl.fitted"
# View(df_DY_spl)  -> derniere colonne ajoutee

#~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_DY_lmm$Time[which(df_DY_lmm$BRAS==0)],
     df_DY_lmm$DY_lmm.fitted[which(df_DY_lmm$BRAS==0)], 
     ylim = c(0,30), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Dyspnee (trajctoires linéaires)", 
     cex.main = 1, type="l")
lines(df_DY_lmm$Time[which(df_DY_lmm$BRAS==1)],
      df_DY_lmm$DY_lmm.fitted[which(df_DY_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 7, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_DY_spl$Time[which(df_DY_spl$BRAS==0)], 
     df_DY_spl$DY_spl.fitted[which(df_DY_spl$BRAS==0)], 
     col=c("blue"), ylim=c(0,30), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Dyspnee (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_DY_spl$Time[which(df_DY_spl$BRAS==1)], 
      df_DY_spl$DY_spl.fitted[which(df_DY_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 7, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                             Insomnies
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.SL$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees SL manquantes:
df_SL_lmm <- df_tpreel[which(!is.na(df_tpreel$SL)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_SL_lmm <-  subset(df_SL_lmm, select =-c(visit:LagDate, QL2:DY, AP:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_SL_lmm )
row.names(df_SL_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_SL_lmm <- cbind(df_SL_lmm, fit.lme.SL$fitted[,1])
names(df_SL_lmm)[ncol(df_SL_lmm)] <- "SL_lmm.fitted"
# View(df_SL_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne SL:
df_SL_spl <- df_tpreel[which(!is.na(df_tpreel$SL)),]
#on ne garde que les colonnes utiles à l'échelle:
df_SL_spl <-  subset(df_SL_spl, select =-c(visit:LagDate, QL2:DY, AP:FI))
df_SL_spl <- df_SL_spl[order(df_SL_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_SL_spl)
row.names(df_SL_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.SL.1 <- lme(SL ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_SL_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.SL.2 <- lme(SL ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_SL_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.SL.1, fit.lme.ns.SL.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.SL.1.ml <- update(fit.lme.ns.SL.1, method="ML")
fit.lme.ns.SL.2.ml <- update(fit.lme.ns.SL.2, method="ML")
SL_test <- lrtest(fit.lme.ns.SL.1.ml, fit.lme.ns.SL.2.ml) 
SL_pval <- SL_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_SL_spl <- cbind(df_SL_spl, fit.lme.ns.SL.2$fitted[,1])
names(df_SL_spl)[ncol(df_SL_spl)] <- "SL_spl.fitted"
# View(df_SL_spl)  -> derniere colonne ajoutee

#~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_SL_lmm$Time[which(df_SL_lmm$BRAS==0)],
     df_SL_lmm$SL_lmm.fitted[which(df_SL_lmm$BRAS==0)], 
     ylim = c(30,50), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Insomnie (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_SL_lmm$Time[which(df_SL_lmm$BRAS==1)],
      df_SL_lmm$SL_lmm.fitted[which(df_SL_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 50, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_SL_spl$Time[which(df_SL_spl$BRAS==0)], 
     df_SL_spl$SL_spl.fitted[which(df_SL_spl$BRAS==0)], 
     col=c("blue"), ylim=c(30,50), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Insomnie (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_SL_spl$Time[which(df_SL_spl$BRAS==1)], 
      df_SL_spl$SL_spl.fitted[which(df_SL_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 50, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                       Perte d'appetit
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.AP$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees AP manquantes:
df_AP_lmm <- df_tpreel[which(!is.na(df_tpreel$AP)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_AP_lmm <-  subset(df_AP_lmm, select =-c(visit:LagDate, QL2:SL, CO:FI))
#vérfication des dimensions des tableaux:
dim(df_tpreel)
dim(df_AP_lmm )
row.names(df_AP_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_AP_lmm <- cbind(df_AP_lmm, fit.lme.AP$fitted[,1])
names(df_AP_lmm)[ncol(df_AP_lmm)] <- "AP_lmm.fitted"
# View(df_AP_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne AP:
df_AP_spl <- df_tpreel[which(!is.na(df_tpreel$AP)),]
#on ne garde que les colonnes utiles à l'échelle:
df_AP_spl <-  subset(df_AP_spl, select =-c(visit:LagDate, QL2:sL, CO:FI))
df_AP_spl <- df_AP_spl[order(df_AP_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_AP_spl)
row.names(df_AP_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.AP.1 <- lme(AP ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_AP_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.AP.2 <- lme(AP ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_AP_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.AP.1, fit.lme.ns.AP.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.AP.1.ml <- update(fit.lme.ns.AP.1, method="ML")
fit.lme.ns.AP.2.ml <- update(fit.lme.ns.AP.2, method="ML")
AP_test <- lrtest(fit.lme.ns.AP.1.ml, fit.lme.ns.AP.2.ml) 
AP_pval <- AP_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_AP_spl <- cbind(df_AP_spl, fit.lme.ns.AP.2$fitted[,1])
names(df_AP_spl)[ncol(df_AP_spl)] <- "AP_spl.fitted"
# View(df_AP_spl)  -> derniere colonne ajoutee

#~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_AP_lmm$Time[which(df_AP_lmm$BRAS==0)],
     df_AP_lmm$AP_lmm.fitted[which(df_AP_lmm$BRAS==0)], 
     ylim = c(0,20), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Perte d'Appetit (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_AP_lmm$Time[which(df_AP_lmm$BRAS==1)],
      df_AP_lmm$AP_lmm.fitted[which(df_AP_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 5, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_AP_spl$Time[which(df_AP_spl$BRAS==0)], 
     df_AP_spl$AP_spl.fitted[which(df_AP_spl$BRAS==0)], 
     col=c("blue"), ylim=c(0,20), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Perte d'Appetit (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_AP_spl$Time[which(df_AP_spl$BRAS==1)], 
      df_AP_spl$AP_spl.fitted[which(df_AP_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 5, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                       Constipations
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.CO$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees CO manquantes:
df_CO_lmm <- df_tpreel[which(!is.na(df_tpreel$CO)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_CO_lmm <-  subset(df_CO_lmm, select =-c(visit:LagDate, QL2:AP, DI:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_CO_lmm )
row.names(df_CO_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_CO_lmm <- cbind(df_CO_lmm, fit.lme.CO$fitted[,1])
names(df_CO_lmm)[ncol(df_CO_lmm)] <- "CO_lmm.fitted"
# View(df_CO_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne CO:
df_CO_spl <- df_tpreel[which(!is.na(df_tpreel$CO)),]
#on ne garde que le colonnes utiles à l'échelles
df_CO_spl <-  subset(df_CO_spl, select =-c(visit:LagDate, QL2:AP, DI:FI))
df_CO_spl <- df_CO_spl[order(df_CO_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_CO_spl)
row.names(df_CO_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.CO.1 <- lme(CO ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_CO_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.CO.2 <- lme(CO ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_CO_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.CO.1, fit.lme.ns.CO.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.CO.1.ml <- update(fit.lme.ns.CO.1, method="ML")
fit.lme.ns.CO.2.ml <- update(fit.lme.ns.CO.2, method="ML")
CO_test <- lrtest(fit.lme.ns.CO.1.ml, fit.lme.ns.CO.2.ml) 
CO_pval <- CO_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_CO_spl <- cbind(df_CO_spl, fit.lme.ns.CO.2$fitted[,1])
names(df_CO_spl)[ncol(df_CO_spl)] <- "CO_spl.fitted"
# View(df_CO_spl)  -> derniere colonne ajoutee

#~~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_CO_lmm$Time[which(df_CO_lmm$BRAS==0)],
     df_CO_lmm$CO_lmm.fitted[which(df_CO_lmm$BRAS==0)], 
     ylim = c(10,30), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Constipation (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_CO_lmm$Time[which(df_CO_lmm$BRAS==1)],
      df_CO_lmm$CO_lmm.fitted[which(df_CO_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 30, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_CO_spl$Time[which(df_CO_spl$BRAS==0)], 
     df_CO_spl$CO_spl.fitted[which(df_CO_spl$BRAS==0)], 
     col=c("blue"), ylim=c(10,30), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Consitpation (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_CO_spl$Time[which(df_CO_spl$BRAS==1)], 
      df_CO_spl$CO_spl.fitted[which(df_CO_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 30, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                       Diarrhees
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.DI$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees DI manquantes:
df_DI_lmm <- df_tpreel[which(!is.na(df_tpreel$DI)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_DI_lmm <-  subset(df_DI_lmm, select =-c(visit:LagDate, QL2:CO, FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_DI_lmm )
row.names(df_DI_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_DI_lmm <- cbind(df_DI_lmm, fit.lme.DI$fitted[,1])
names(df_DI_lmm)[ncol(df_DI_lmm)] <- "DI_lmm.fitted"
# View(df_DI_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne DI:
df_DI_spl <- df_tpreel[which(!is.na(df_tpreel$DI)),]
#on ne garde que les colonnes utiles à l'échelle:
df_DI_spl <-  subset(df_DI_spl, select =-c(visit:LagDate, QL2:CO, FI))
df_DI_spl <- df_DI_spl[order(df_DI_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux
dim(df_tpreel)
dim(df_DI_spl)
row.names(df_DI_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.DI.1 <- lme(DI ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_DI_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.DI.2 <- lme(DI ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_DI_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.DI.1, fit.lme.ns.DI.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.DI.1.ml <- update(fit.lme.ns.DI.1, method="ML")
fit.lme.ns.DI.2.ml <- update(fit.lme.ns.DI.2, method="ML")
DI_test <- lrtest(fit.lme.ns.DI.1.ml, fit.lme.ns.DI.2.ml) 
DI_pval <- DI_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_DI_spl <- cbind(df_DI_spl, fit.lme.ns.DI.2$fitted[,1])
names(df_DI_spl)[ncol(df_DI_spl)] <- "DI_spl.fitted"
# View(df_DI_spl)  -> derniere colonne ajoutee

#~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_DI_lmm$Time[which(df_DI_lmm$BRAS==0)],
     df_DI_lmm$DI_lmm.fitted[which(df_DI_lmm$BRAS==0)], 
     ylim = c(0,15), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Diarrhee (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_DI_lmm$Time[which(df_DI_lmm$BRAS==1)],
      df_DI_lmm$DI_lmm.fitted[which(df_DI_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 4, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_DI_spl$Time[which(df_DI_spl$BRAS==0)], 
     df_DI_spl$DI_spl.fitted[which(df_DI_spl$BRAS==0)], 
     col=c("blue"), ylim=c(0,15), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Diarrhee (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_DI_spl$Time[which(df_DI_spl$BRAS==1)], 
      df_DI_spl$DI_spl.fitted[which(df_DI_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 4, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#------------------------------------------------------------------------------
#                      Difficultees financieres
#------------------------------------------------------------------------------

#~~~~~~~~~~~~~Linéaires~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(fit.lme.FI$fitted[,1]) #vérification du nombre de lignes : données 
                              #manquantes non prises en compte
#on construit un fichier ou on supprime les lignes avec donnees FI manquantes:
df_FI_lmm <- df_tpreel[which(!is.na(df_tpreel$FI)),] # retire (si il y en a)
#on ne garde que les colonnes utiles à l'échelle:
df_FI_lmm <-  subset(df_FI_lmm, select =-c(visit:LagDate, QL2:DI))
#vérifivation des dimensions des tableaux:
dim(df_tpreel)
dim(df_FI_lmm )
row.names(df_FI_lmm) <- NULL

#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_FI_lmm <- cbind(df_FI_lmm, fit.lme.FI$fitted[,1])
names(df_FI_lmm)[ncol(df_FI_lmm)] <- "FI_lmm.fitted"
# View(df_FI_lmm)  -> derniere colonne ajoutee

#~~~~~~~~~~~~~Splines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-- Modele --------------------------------------------------------------------
#pour ne pas avoir les NA dans la colonne FI:
df_FI_spl <- df_tpreel[which(!is.na(df_tpreel$FI)),]
#on ne garde que les colonnes utiles à l'échelle:
df_FI_spl <-  subset(df_FI_spl, select =-c(visit:LagDate, QL2:DI))
df_FI_spl <- df_FI_spl[order(df_FI_spl$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_FI_spl)
row.names(df_FI_spl) <- NULL

#-- Comparaison modeles splines -----------------------------------------------
fit.lme.ns.FI.1 <- lme(FI ~ ns(Time,3), 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_FI_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #sans effet bras

fit.lme.ns.FI.2 <- lme(FI ~ ns(Time,3):BRAS, 
                       random = list(npat = pdDiag(form = ~ ns(Time, 3))),
                 #matrice de variance covariance des effets aleatoires diagonale
                       data=df_FI_spl, na.action = na.omit,
                       control = lmeControl(opt="optim")) #avec effet bras

lrtest(fit.lme.ns.FI.1, fit.lme.ns.FI.2)
#p-value tres petite = effet significatif de l'effet bras

#Tester les modeles
fit.lme.ns.FI.1.ml <- update(fit.lme.ns.FI.1, method="ML")
fit.lme.ns.FI.2.ml <- update(fit.lme.ns.FI.2, method="ML")
FI_test <- lrtest(fit.lme.ns.FI.1.ml, fit.lme.ns.FI.2.ml) 
FI_pval <- FI_test$`Pr(>Chisq)`[2]

#-- Valeurs predites ----------------------------------------------------------
df_FI_spl <- cbind(df_FI_spl, fit.lme.ns.FI.2$fitted[,1])
names(df_FI_spl)[ncol(df_FI_spl)] <- "FI_spl.fitted"
# View(df_FI_spl)  -> derniere colonne ajoutee

#~~~~~~~~Graphiques~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(2,1))
#LINEAIRES
plot(df_FI_lmm$Time[which(df_FI_lmm$BRAS==0)],
     df_FI_lmm$FI_lmm.fitted[which(df_FI_lmm$BRAS==0)], 
     ylim = c(5,30), col="blue", ylab="Score moyen prédit (%)", 
     xlab="Temps (mois)", cex.lab=0.8, 
     main = "Difficultees financieres (trajectoires linéaires)", 
     cex.main = 1, type="l")
lines(df_FI_lmm$Time[which(df_FI_lmm$BRAS==1)],
      df_FI_lmm$FI_lmm.fitted[which(df_FI_lmm$BRAS==1)], 
      col="red", type="l")
legend(0, 30, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)

#SPLINES
plot(df_FI_spl$Time[which(df_FI_spl$BRAS==0)], 
     df_FI_spl$FI_spl.fitted[which(df_FI_spl$BRAS==0)], 
     col=c("blue"), ylim=c(5,30), cex.lab=0.8, 
     xlab = "Temps (mois)", ylab = "Score moyen prédit (%)", 
     main = "Difficultees financieres (trajectoires splines)", cex.main = 1, 
     type = "l")
lines(df_FI_spl$Time[which(df_FI_spl$BRAS==1)], 
      df_FI_spl$FI_spl.fitted[which(df_FI_spl$BRAS==1)], 
      col="red", type = "l")
legend(0, 30, lty=c(1,1), col=c("red","blue"), 
       legend=c("APAD (exp)", "Usual Care (ctrl)"), cex = 0.6)


#------------------------------------------------------------------------------
#                       TABLEAU FINAL DES SPLINES
#------------------------------------------------------------------------------

TAB_spl <- data.frame(matrix("", 17, 2))
TAB_spl <- data.frame(c("Statut Global Sante", "ECHELLES FONCTIONNELLES", 
                    "Physique", "Personnelles", "Emotionnelles", "Cognitives",
                    "Sociales", "ECHELLES SYMPTOMATIQUES", "Fatigue", 
                    "Nausees et vomissements", "Douleurs", "Dyspnee",
                    "Insomnies", "Pertes d'appetit", "Constipations", 
                    "Diarrhees", "Difficultees Financieres"), 
                  c(round(QL2_pval, 5), "-----", round(PF_pval, 5),
                    round(RF_pval, 5), round(EF_pval, 5), round(CF_pval, 5), 
                    round(SF_pval, 5), "-----", round(FA_pval, 5), 
                    round(NV_pval, 5), round(PA_pval, 5), round(DY_pval, 5), 
                    round(SL_pval, 5), round(AP_pval, 5), round(CO_pval, 5), 
                    round(DI_pval, 5), round(FI_pval, 5)))
names(TAB_spl) <- c("", "p-value")
View(TAB_spl)
# write.table(TAB_spl, file = "tableau_estimations_mois_SPLINES.csv",
#             sep = ";", dec = ",", row.names = FALSE)


###############################################################################
#         Test d'hypotheses sur les 4 échelles choisies                       #
###############################################################################

#------------------------------------------------------------------------------
#            Statut Global de santé
#------------------------------------------------------------------------------

#~~~Distribution residus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
#pour ne pas avoir les NA dans la colonne QL2:
df_QL2_res <- df_tpreel[which(!is.na(df_tpreel$QL2)),]
#on ne garde que les colonnes utiles pour l'échelle:
df_QL2_res <-  subset(df_QL2_res, select =-c(visit:LagDate, PF2:FI))
df_QL2_res <- df_QL2_res[order(df_QL2_res$Time), ] # ordonner par les temps
#vérifier les dimensions des tableaux:
dim(df_tpreel)
dim(df_QL2_res)
row.names(df_QL2_res) <- NULL
residus_QL2 <- residuals(fit.lme.QL2) #extraction des résidus du modèle
df_QL2_res$residus_QL2 <- residus_QL2 #ajout des résidus au tableau

#---Graphique------------------------------------------------------------------
p_QL2_res <- ggplot(df_QL2_res, aes(x=residus_QL2)) + 
  geom_histogram(aes(y=..density..))
p_QL2_res <- p_QL2_res + 
  stat_function(fun = dnorm, aes(colour = "loi normale"),
                args = list(mean = mean(df_QL2_res$residus_QL2, na.rm = TRUE), 
                            sd = sd(df_QL2_res$residus_QL2, na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_QL2_res <- p_QL2_res + theme(plot.title = element_text(size = 12, 
                                                         face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus du \n statut global de santé: linéaire")

#SPLINES
#pour ne pas avoir les NA dans la colonne QL2:
df_QL2_res_spl <- df_tpreel[which(!is.na(df_tpreel$QL2)),]
#on ne garde que les colonnes à l'échelle:
df_QL2_res_spl <-  subset(df_QL2_res_spl, select =-c(visit:LagDate, PF2:FI))
#ordonner par les temps:
df_QL2_res_spl <- df_QL2_res_spl[order(df_QL2_res_spl$Time), ]
#vérifier les dimensions des tableaux:
dim(df_tpreel)
dim(df_QL2_res_spl)
row.names(df_QL2_res_spl) <- NULL
residus_QL2_spl <- residuals(fit.lme.ns.QL2.2) #extraction des résidus du modèle
df_QL2_res_spl$residus_QL2_spl <- residus_QL2_spl #ajout au tableau

#---Graphique------------------------------------------------------------------
p_QL2_res_spl <- ggplot(df_QL2_res_spl, aes(x=residus_QL2_spl)) + 
  geom_histogram(aes(y=..density..))
p_QL2_res_spl <- p_QL2_res_spl + 
  stat_function(fun = dnorm, aes(colour = "loi normale"),
                        args = list(mean = mean(df_QL2_res_spl$residus_QL2_spl, 
                                                                 na.rm = TRUE),
                                        sd = sd(df_QL2_res_spl$residus_QL2_spl, 
                                                               na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_QL2_res_spl <- p_QL2_res_spl + 
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus du \n statut global de santé : splines")

grid.arrange(p_QL2_res, p_QL2_res_spl)

#~~~QQ-PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,2))
#LINEAIRES
QL2_qqplot <- qqPlot(df_QL2_res$residus_QL2, envelope = .95, id=FALSE, 
                     xlab = "Quantiles de la loi normale",
                     ylab = "Résidus standardisés", 
                     main = "QQ-plot statut global \n de santé: linéaire")
#SPLINES
QL2_qqplot_spl <- qqPlot(df_QL2_res_spl$residus_QL2_spl, envelope = .95, 
                         id=FALSE, xlab = "Quantiles de la loi normale",
                         ylab = "Résidus standardisés", 
                         main = "QQ-plot statut global \n de santé : spline")

#~~~Graphe des résidus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
QL2_res <- plot(fit.lme.QL2, 
                main="Résidus du statut global de santé :\n linéaire")
#SPLINES
QL2_res_spl <- plot(fit.lme.ns.QL2.2, 
                    main="Résidus du statut global de santé : \n spline")

grid.arrange(QL2_res, QL2_res_spl)


#------------------------------------------------------------------------------
#            Fonction physique
#------------------------------------------------------------------------------

#~~~Distribution residus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
#pour ne pas avoir les NA dans la colonne PF:
df_PF_res <- df_tpreel[which(!is.na(df_tpreel$PF2)),]
#on ne garde que les colonnes utiles à l'échelle:
df_PF_res <-  subset(df_PF_res, select =-c(visit:LagDate, QL2, RF2:FI))
df_PF_res <- df_PF_res[order(df_PF_res$Time), ] # ordonner par les temps
#vérifier les dimensions des tableaux
dim(df_tpreel)
dim(df_PF_res)
row.names(df_PF_res) <- NULL
residus_PF <- residuals(fit.lme.PF) #extraction des résidus du modèle
df_PF_res$residus_PF <- residus_PF #ajout des résidus au tableau

#---Graphique------------------------------------------------------------------
p_PF_res <- ggplot(df_PF_res, aes(x=residus_PF)) + 
  geom_histogram(aes(y=..density..))
p_PF_res <- p_PF_res + 
  stat_function(fun = dnorm, aes(colour = "loi normale"), 
                args = list(mean = mean(df_PF_res$residus_PF, na.rm = TRUE), 
                            sd = sd(df_PF_res$residus_PF, na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_PF_res <- p_PF_res + theme(plot.title = element_text(size = 12, 
                                                       face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus de la \n fonction physique: linéaire")

#SPLINES
#pour ne pas avoir les NA dans la colonne PF:
df_PF_res_spl <- df_tpreel[which(!is.na(df_tpreel$PF2)),]
#on ne garde que les colonnes utiles à l'échelle:
df_PF_res_spl <-  subset(df_PF_res_spl, select =-c(visit:LagDate, QL2, RF2:FI))
#ordonner par les temps:
df_PF_res_spl <- df_PF_res_spl[order(df_PF_res_spl$Time), ]
#vérifier les dimensions des tableaux:
dim(df_tpreel)
dim(df_PF_res_spl)
row.names(df_PF_res_spl) <- NULL
residus_PF_spl <- residuals(fit.lme.ns.PF.2) #extraction des résidus du modèle
df_PF_res_spl$residus_PF_spl <- residus_PF_spl #ajout au tableau

#---Graphique------------------------------------------------------------------
p_PF_res_spl <- ggplot(df_PF_res_spl, aes(x=residus_PF_spl)) + 
  geom_histogram(aes(y=..density..))
p_PF_res_spl <- p_PF_res_spl + 
  stat_function(fun = dnorm, aes(colour = "loi normale"), 
                args = list(mean = mean(df_PF_res_spl$residus_PF_spl, 
                                        na.rm = TRUE), 
                            sd = sd(df_PF_res_spl$residus_PF_spl, 
                                    na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_PF_res_spl <- p_PF_res_spl + theme(plot.title = element_text(size = 12, 
                                                               face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus de la \n fonction physique : splines")
#print(p_PF_res_spl)

grid.arrange(p_PF_res, p_PF_res_spl)

#~~~QQ-PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,2))
#LINEAIRES
PF_qqplot <- qqPlot(df_PF_res$residus_PF, envelope = .95, id=FALSE, 
                    xlab = "Quantiles de la loi normale",
                    ylab = "Résidus standardisés", 
                    main = "QQ-plot fonction \n physique : linéaire")

#SPLINES
PF_qqplot_spl <- qqPlot(df_PF_res_spl$residus_PF_spl, envelope = .95, id=FALSE, 
                        xlab = "Quantiles de la loi normale",
                        ylab = "Résidus standardisés", 
                        main = "QQ-plot fonction \n physique : spline")

#~~~Graphe des résidus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
PF_res <- plot(fit.lme.PF, 
               main="Résidus de la fonction physique :\n linéaire")
#SPLINES
PF_res_spl <- plot(fit.lme.ns.PF.2, 
                   main="Résidus de la fonction physique : \n spline")

grid.arrange(PF_res, PF_res_spl)

#------------------------------------------------------------------------------
#            Fonction personelle
#------------------------------------------------------------------------------

#~~~Distribution residus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
#pour ne pas avoir les NA dans la colonne RF:
df_RF_res <- df_tpreel[which(!is.na(df_tpreel$RF2)),]
#on ne garde que les colonnes utiles à l'échelle:
df_RF_res <-  subset(df_RF_res, select =-c(visit:LagDate, QL2:PF2, EF:FI))
df_RF_res <- df_RF_res[order(df_RF_res$Time), ] # ordonner par les temps
#vérifier les dimensions des tableaux:
dim(df_tpreel)
dim(df_RF_res)
row.names(df_RF_res) <- NULL
residus_RF <- residuals(fit.lme.RF) #extraction des résidus du modèle
df_RF_res$residus_RF <- residus_RF #ajout des résidus au tableau

#---Graphique------------------------------------------------------------------
p_RF_res <- ggplot(df_RF_res, aes(x=residus_RF)) + 
  geom_histogram(aes(y=..density..))
p_RF_res <- p_RF_res + 
  stat_function(fun = dnorm, aes(colour = "loi normale"), 
                args = list(mean = mean(df_RF_res$residus_RF, na.rm = TRUE), 
                            sd = sd(df_RF_res$residus_RF, na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_RF_res <- p_RF_res + theme(plot.title = element_text(size = 12, 
                                                       face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus de la \n fonction personnelle: linéaire")

#SPLINES
#pour ne pas avoir les NA dans la colonne RF:
df_RF_res_spl <- df_tpreel[which(!is.na(df_tpreel$RF2)),]
#on ne garde que les colonnes utiles à l'échelle:
df_RF_res_spl <-  subset(df_RF_res_spl, 
                         select =-c(visit:LagDate, QL2:PF2, EF:FI))
#ordonner par les temps:
df_RF_res_spl <- df_RF_res_spl[order(df_RF_res_spl$Time), ]
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_RF_res_spl)
row.names(df_RF_res_spl) <- NULL
residus_RF_spl <- residuals(fit.lme.ns.RF.2) #extraction des résidus du modèle
df_RF_res_spl$residus_RF_spl <- residus_RF_spl #ajout au tableau

#---Graphique------------------------------------------------------------------
p_RF_res_spl <- ggplot(df_RF_res_spl, aes(x=residus_RF_spl)) + 
  geom_histogram(aes(y=..density..))
p_RF_res_spl <- p_RF_res_spl + 
  stat_function(fun = dnorm, aes(colour = "loi normale"), 
                args = list(mean = mean(df_RF_res_spl$residus_RF_spl, 
                                        na.rm = TRUE), 
                            sd = sd(df_RF_res_spl$residus_RF_spl, 
                                    na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_RF_res_spl <- p_RF_res_spl + theme(plot.title = element_text(size = 12, 
                                                               face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus de la \n fonction personnelle : splines")

grid.arrange(p_RF_res, p_RF_res_spl)


#~~~QQ-PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,2))
#LINEAIRES
RF_qqplot <- qqPlot(df_RF_res$residus_RF, envelope = .95, id=FALSE, 
                    xlab = "Quantiles de la loi normale",
                    ylab = "Résidus standardisés", 
                    main = "QQ-plot fonction \n personnelle : linéaire")
#SPLINES
RF_qqplot_spl <- qqPlot(df_RF_res_spl$residus_RF_spl, envelope = .95, id=FALSE, 
                        xlab = "Quantiles de la loi normale",
                        ylab = "Résidus standardisés", 
                        main = "QQ-plot fonction \n personelle : spline")

#~~~Graphe des résidus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
RF_res <- plot(fit.lme.RF, 
               main="Résidus de la fonction personnelle :\n linéaire")
#SPLINES
RF_res_spl <- plot(fit.lme.ns.RF.2, 
                   main="Résidus de la fonction personnelle : \n spline")

grid.arrange(RF_res, RF_res_spl)


#------------------------------------------------------------------------------
#            Fonction personelle
#------------------------------------------------------------------------------

#~~~Distribution residus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINEAIRES
#pour ne pas avoir les NA dans la colonne SF:
df_SF_res <- df_tpreel[which(!is.na(df_tpreel$SF)),]
#on ne garde que les colonnes utiles à l'échelle:
df_SF_res <-  subset(df_SF_res, select =-c(visit:LagDate, QL2:CF, FA:FI))
df_SF_res <- df_SF_res[order(df_SF_res$Time), ] # ordonner par les temps
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_SF_res)
row.names(df_SF_res) <- NULL
residus_SF <- residuals(fit.lme.SF) #extraction des résidus du modèle
df_SF_res$residus_SF <- residus_SF #ajout des résidus dans le tableau

#---Graphique------------------------------------------------------------------
p_SF_res <- ggplot(df_SF_res, aes(x=residus_SF)) + 
  geom_histogram(aes(y=..density..))
p_SF_res <- p_SF_res + 
  stat_function(fun = dnorm, aes(colour = "loi normale"), 
                args = list(mean = mean(df_SF_res$residus_SF, na.rm = TRUE), 
                            sd = sd(df_SF_res$residus_SF, na.rm = TRUE)))+
  scale_colour_manual("Densité", values = c("red"))
p_SF_res <- p_SF_res + theme(plot.title = element_text(size = 12, 
                                                       face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité") + 
  ggtitle("Histogramme des résidus de la \n fonction sociale: linéaire")

#SPLINES
#pour ne pas avoir les NA dans la colonne SF:
df_SF_res_spl <- df_tpreel[which(!is.na(df_tpreel$SF)),]
#on ne garde que les colonnes utiles à l'échelle:
df_SF_res_spl <-  subset(df_SF_res_spl, 
                         select =-c(visit:LagDate, QL2:CF, FA:FI))
#ordonner par les temps
df_SF_res_spl <- df_SF_res_spl[order(df_SF_res_spl$Time), ]
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_SF_res_spl)
row.names(df_SF_res_spl) <- NULL
residus_SF_spl <- residuals(fit.lme.ns.SF.2) #extraction des résidus du modèle
df_SF_res_spl$residus_SF_spl <- residus_SF_spl #ajout des résidus au tableau

#---Graphique------------------------------------------------------------------
p_SF_res_spl <- ggplot(df_SF_res_spl, aes(x=residus_SF_spl)) + 
  geom_histogram(aes(y=..density..))
p_SF_res_spl <- p_SF_res_spl + 
  stat_function(fun = dnorm, aes(colour = "loi normale"), 
                args = list(mean = mean(df_SF_res_spl$residus_SF_spl, 
                                        na.rm = TRUE), 
                            sd = sd(df_SF_res_spl$residus_SF_spl, 
                                    na.rm = TRUE))) +
  scale_colour_manual("Densité", values = c("red"))
p_SF_res_spl <- p_SF_res_spl + theme(plot.title = element_text(size = 12, 
                                                               face = "bold")) +
  theme(legend.position="right") +
  labs(x = "Distribution du résidu", y = "Densité")+ 
  ggtitle("Histogramme des résidus de la \n fonction sociale : splines")

grid.arrange(p_SF_res, p_SF_res_spl, ncol=2)

#~~~QQ-PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,2))
#LINEAIRES
RF_qqplot <- qqPlot(df_SF_res$residus_SF, envelope = .95, id=FALSE, 
                    xlab = "Quantiles de la loi normale",
                    ylab = "Résidus standardisés", 
                    main = "QQ-plot fonction \n sociale : linéaire")
#SPLINES
SF_qqplot_spl <- qqPlot(df_SF_res_spl$residus_SF_spl, envelope = .95, id=FALSE, 
                        xlab = "Quantiles de la loi normale",
                        ylab = "Résidus standardisés", 
                        main = "QQ-plot fonction \n sociale : spline")

#~~~Graphe des résidus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LINAIRES
SF_res <- plot(fit.lme.SF, 
               main="Résidus de la fonction sociale :\n linéaire")
#SPLINES
SF_res_spl <- plot(fit.lme.ns.SF.2, 
                   main="Résidus de la fonction sociale : \n spline")

grid.arrange(SF_res, SF_res_spl)

###############################################################################
#                         SHAPIRO TEST                                        #
###############################################################################

#-----------QL2----------------------------------------------------------------
shapiro_QL2_lmm <- shapiro.test(residuals(fit.lme.QL2))
shapiro_QL2_spl <- shapiro.test(residuals(fit.lme.ns.QL2.2))

#-----------PF-----------------------------------------------------------------
shapiro_PF_lmm <- shapiro.test(residuals(fit.lme.PF))
shapiro_PF_spl <- shapiro.test(residuals(fit.lme.ns.PF.2))

#-----------RF-----------------------------------------------------------------
shapiro_RF_lmm <- shapiro.test(residuals(fit.lme.RF))
shapiro_RF_spl <- shapiro.test(residuals(fit.lme.ns.RF.2))

#-----------SF-----------------------------------------------------------------
shapiro_SF_lmm <- shapiro.test(residuals(fit.lme.SF))
shapiro_SF_spl <- shapiro.test(residuals(fit.lme.ns.SF.2))

#------------TABLEAU-----------------------------------------------------------
TAB_shapiro <- data.frame(matrix("", 4, 3))
TAB_shapiro <- data.frame(c("Statut Global Sante", "Fonction Physique", 
                            "Fonction Personnelle", "Fonction Sociale"), 
                          c(shapiro_QL2_lmm$p.value, shapiro_PF_lmm$p.value,
                            shapiro_RF_lmm$p.value, shapiro_SF_lmm$p.value),
                          c(shapiro_QL2_spl$p.value, shapiro_PF_spl$p.value,
                            shapiro_RF_spl$p.value, shapiro_SF_spl$p.value)
)
names(TAB_shapiro) <- c("p-value Shapiro test", "Linéaire", "Spline")
View(TAB_shapiro)
