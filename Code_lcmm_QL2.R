###############################################################################
#                             PACKAGES UTILES                                 #
###############################################################################

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
library(car)

###############################################################################
#                           LCMM                                              #
###############################################################################

#------------------------------------------------------------------------------
#                 DIMENSION QL2
#------------------------------------------------------------------------------

#---Trajectoires linéaires + lien linéaire-------------------------------------
lien <- "linear" #choix de la fonction de lien pour la transformation
QL2_lin_lin <- lcmm(QL2 ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                    data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_QL2_ll <- subset(df_tpreel, select =-c(visit:LagDate, PF2:FI))
#on retire les lignes avec des NA:
df_QL2_ll <- df_QL2_ll[which(!is.na(df_QL2_ll$QL2)),]
#ajout du score transformé:
df_QL2_ll <- cbind(df_QL2_ll, QL2_lin_lin$pred$pred_m)
names(df_QL2_ll)[names(df_QL2_ll) == "QL2_lin_lin$pred$pred_m"] <- "score_transf_QL2_lin"
df_QL2_ll <- df_QL2_ll[order(df_QL2_ll$Time), ] #ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
QL2_ll_p <- ggplot(df_QL2_ll, 
                   aes(x=Time, y=score_transf_QL2_lin, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue"))+
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires linéaires et lien linéaire (LCMM) :
          Statut global de santé")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(QL2_ll_p)

#LMM
length(fit.lme.QL2$fitted[,1]) #vérification du nombre de lignes : données 
#manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données QL2 manquantes:
df_QL2_lmm <- df_tpreel[which(!is.na(df_tpreel$QL2)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_QL2_lmm <-  subset(df_QL2_lmm, select =-c(visit:LagDate, PF2:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_QL2_lmm)
row.names(df_QL2_lmm) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_QL2_lmm <- cbind(df_QL2_lmm, fit.lme.QL2$fitted[,1])
names(df_QL2_lmm)[ncol(df_QL2_lmm)] <- "QL2_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees QL2 manquantes :
df_QL2_lcmm <- df_tpreel[which(!is.na(df_tpreel$QL2)),] # retire (si il y en a)
df_QL2_lcmm <-  subset(df_QL2_lcmm, select =-c(visit:LagDate, PF2:FI))
df_QL2_lcmm <- df_QL2_lcmm[order(df_QL2_lcmm$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_QL2_lcmm)
row.names(df_QL2_lcmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(QL2_lin_lin, newdata = df_QL2_lcmm)
df_QL2_lcmm <- cbind(df_QL2_lcmm, pred_marg_transform$pred)
names(df_QL2_lcmm)[ncol(df_QL2_lcmm)] <- "QL2_lcmm.fitted"

pQL2_lcmm <- ggplot(data = df_QL2_lmm, aes(x=Time, QL2_lmm.fitted, 
                                           group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_QL2_lcmm,aes(x=Time, y=QL2_lcmm.fitted, 
                                   group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(40,80))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Statut global de santé :
  Trajectoires linéaires (LMM) : pleins
  Trajectoires linéaires et lien linéaires (LCMM) : pointillés")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pQL2_lcmm)

### Test des hypothèses #######################################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.QL2)
resid_transform <- residuals(QL2_lin_lin)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LMM): \n Statut global de santé")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LCMM): \n Statut global de santé")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot statut global \n de santé: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE,
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot statut global \n de santé: LCMM")


#Résidu prédit vs modèle prédit
QL2_lmm_res <- plot(model_initial_QL2_spl,
                    main="Graphique des résidus (LCMM): \n Statut global de santé")
QL2_lcmm_res <- plot(model_transform_QL2_spl,
                     main="Graphique des résidus (LMM): \n Statut global de santé")
grid.arrange(QL2_lmm_res, QL2_lcmm_res)



#--Trajectoires lineaires + lien spline----------------------------------------
lien <- "splines"
QL2_lin_spl <- lcmm(QL2 ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                    data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_QL2_ls <- subset(df_tpreel, select =-c(visit:LagDate, PF2:FI))
#on retire les lignes avec des NA:
df_QL2_ls <- df_QL2_ls[which(!is.na(df_QL2_ls$QL2)),]
#ajout du score transformé:
df_QL2_ls <- cbind(df_QL2_ls, QL2_lin_spl$pred$pred_m)
names(df_QL2_ls)[names(df_QL2_ls) == "QL2_lin_spl$pred$pred_m"] <- "score_transf_QL2_lin_spl"
df_QL2_ls <- df_QL2_ls[order(df_QL2_ls$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
QL2_ls_p <- ggplot(df_QL2_ls, aes(x=Time, y=score_transf_QL2_lin_spl, 
                                  group = BRAS)) + geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires linéaires et lien spline (LCMM) :
          Statut global de santé")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(QL2_ls_p)

#LMM
length(fit.lme.QL2$fitted[,1]) #vérification du nombre de lignes : données 
#manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données QL2 manquantes:
df_QL2_lmm <- df_tpreel[which(!is.na(df_tpreel$QL2)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_QL2_lmm <-  subset(df_QL2_lmm, select =-c(visit:LagDate, PF2:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_QL2_lmm )
row.names(df_QL2_lmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
df_QL2_lmm <- cbind(df_QL2_lmm, fit.lme.QL2$fitted[,1])
names(df_QL2_lmm)[ncol(df_QL2_lmm)] <- "QL2_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees QL2 manquantes :
df_QL2_lcmm_ls <- df_tpreel[which(!is.na(df_tpreel$QL2)),] # retire (si il y en a)
df_QL2_lcmm_ls <-  subset(df_QL2_lcmm_ls, select =-c(visit:LagDate, PF2:FI))
df_QL2_lcmm_ls <- df_QL2_lcmm_ls[order(df_QL2_lcmm_ls$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_QL2_lcmm_ls)
row.names(df_QL2_lcmm_ls) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(QL2_lin_spl, newdata = df_QL2_lcmm_ls)
df_QL2_lcmm_ls <- cbind(df_QL2_lcmm_ls, pred_marg_transform$pred)
names(df_QL2_lcmm_ls)[ncol(df_QL2_lcmm_ls)] <- "QL2_lcmm_ls.fitted"

pQL2_lcmm_ls <- ggplot(data = df_QL2_lmm, aes(x=Time, QL2_lmm.fitted, 
                                              group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_QL2_lcmm_ls,aes(x=Time, y=QL2_lcmm_ls.fitted, 
                                      group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(60,80))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Statut global de santé :
            Trajectoires linéaires (LMM) : pleins
          Trajectoires linéaires et lien spline (LCMM) : pointillés")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pQL2_lcmm_ls)

### Test des hypothèses ###########################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.QL2)
resid_transform <- residuals(QL2_lin_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LMM): \n Statut global de santé")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LCMM): \n Statut global de santé")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot statut global \n de santé: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot statut global \n de santé: LCMM")

#Résidu prédit vs modèle prédit
QL2_lmm_res <- plot(fit.lme.QL2,
                    main="Graphique des résidus (LCMM): \n Statut global de santé")
QL2_lcmm_res <- plot(QL2_lin_spl,
                     main="Graphique des résidus (LMM): \n Statut global de santé")
grid.arrange(QL2_lmm_res, QL2_lcmm_res)


#---Trajectoires splines + lien spline-----------------------------------------
lien <- "splines"
#retirer les NA s'il y en a:
df_QL2_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$QL2)),]
QL2_spl_spl <- lcmm(fixed = QL2 ~ ns(Time,3):BRAS, random = ~ 1 + Time,
                    subject = "npat",data=df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_QL2_ss <- subset(df_tpreel, select =-c(visit:LagDate, PF2:FI))
#on retire les lignes avec des NA:
df_QL2_ss <- df_QL2_ss[which(!is.na(df_QL2_ss$QL2)),]
#ajout du score transformé:
df_QL2_ss <- cbind(df_QL2_ss, QL2_spl_spl$pred$pred_m)
names(df_QL2_ss)[names(df_QL2_ss) == "QL2_spl_spl$pred$pred_m"] <- "score_transf_QL2_spl"
df_QL2_ss <- df_QL2_ss[order(df_QL2_ss$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
QL2_ss_p <- ggplot(df_QL2_ss, aes(x=Time, y=score_transf_QL2_spl, 
                                  group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires splines et lien spline (LCMM) :
          Statut global de santé")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(QL2_ss_p)

#LMM
#pour ne pas avoir les NA dans la colonne QL2:
df_QL2_spl <- df_tpreel[which(!is.na(df_tpreel$QL2)),]
#ne garder que les colonnes qui sont utiles pour l'échelle
df_QL2_spl <-  subset(df_QL2_spl, select =-c(visit:LagDate, PF2:FI))
df_QL2_spl <- df_QL2_spl[order(df_QL2_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_QL2_spl)
row.names(df_QL2_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
df_QL2_spl <- cbind(df_QL2_spl, fit.lme.ns.QL2.2$fitted[,1])
names(df_QL2_spl)[ncol(df_QL2_spl)] <- "QL2_spl.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees QL2 manquantes :
df_QL2_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$QL2)),] # retire (si il y en a)
df_QL2_lcmm_spl <-  subset(df_QL2_lcmm_spl, select =-c(visit:LagDate, PF2:FI))
df_QL2_lcmm_spl <- df_QL2_lcmm_spl[order(df_QL2_lcmm_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_QL2_lcmm_spl)
row.names(df_QL2_lcmm_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(QL2_spl_spl, newdata = df_QL2_lcmm_spl)
df_QL2_lcmm_spl <- cbind(df_QL2_lcmm_spl, pred_marg_transform$pred)
names(df_QL2_lcmm_spl)[ncol(df_QL2_lcmm_spl)] <- "QL2_lcmm_spl.fitted"

pQL2_lcmm_spl <- ggplot(data = df_QL2_spl, aes(x=Time, y=QL2_spl.fitted, 
                                               group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_QL2_lcmm_spl,aes(x=Time, y=QL2_lcmm_spl.fitted, 
                                       group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(40,80))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Statut global de santé :
  Trajectoires splines (LMM) : pleins
  Trajectoires splines et lien spline (LCMM) : pointillés") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pQL2_lcmm_spl)

### Test des hypothèses ###########################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.ns.QL2.2)
resid_transform <- residuals(QL2_spl_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LMM): \n Statut global de santé")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LCMM): \n Statut global de santé")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot statut global \n de santé: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot statut global \n de santé: LCMM")


#Résidu prédit vs modèle prédit
QL2_lmm_res <- plot(fit.lme.ns.QL2.2,
                    main="Graphique des résidus (LCMM): \n Statut global de santé")
QL2_lcmm_res <- plot(QL2_spl_spl,
                     main="Graphique des résidus (LMM): \n Statut global de santé")
grid.arrange(QL2_lmm_res, QL2_lcmm_res)


###############################################################################
#                 PLOT CURVILINEARITE
###############################################################################

col <- rainbow(3)
plot(QL2_lin_lin, which = "linkfunction", bty = "l", ylab = "Score prédit (%)", 
     xlab = "Processus latent sous-jacent", lwd = 2, col = col[1], 
     main = "Fonction de lien estimée \n pour le statut global de santé")
plot(QL2_lin_spl, which = "linkfunction", add = TRUE, col = col[2], lwd = 2)
plot(QL2_spl_spl, which = "linkfunction", add = TRUE, col = col[3], lwd = 2)
legend(x = "topleft",
       legend = c("trajectoire linéaire et\n fonction de lien linéaire", 
                  "trajectoire linéaire et\n fonction de lien spline", 
                  "trajectoire spline et\n fonction de lien spline"), lty = 1, 
       cex=0.8, col = col, bty = "n", lwd = 2)