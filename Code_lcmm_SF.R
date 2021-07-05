###############################################################################
#                           LCMM                                              #
###############################################################################

#------------------------------------------------------------------------------
#                 DIMENSION SF
#------------------------------------------------------------------------------

#---Trajectoires linéaires + lien linéaire-------------------------------------
lien <- "linear" #choix de la fonction de lien pour la transformation
SF_lin_lin <- lcmm(SF ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                   data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_SF_ll <- subset(df_tpreel, select =-c(visit:LagDate, QL2:CF, FA:FI))
#on retire les lignes avec des NA:
df_SF_ll <- df_SF_ll[which(!is.na(df_SF_ll$SF)),]
#ajout du score transformé:
df_SF_ll <- cbind(df_SF_ll, SF_lin_lin$pred$pred_m)
names(df_SF_ll)[names(df_SF_ll) == "SF_lin_lin$pred$pred_m"] <- "score_transf_SF_lin"
df_SF_ll <- df_SF_ll[order(df_SF_ll$Time), ] #ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
SF_ll_p <- ggplot(df_SF_ll, 
                  aes(x=Time, y=score_transf_SF_lin, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue"))+
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires linéaires et lien linéaire (LCMM) :
          Fonction sociale")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(SF_ll_p)

#LMM
length(fit.lme.SF$fitted[,1]) #vérification du nombre de lignes : données 
#manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données SF manquantes:
df_SF_lmm <- df_tpreel[which(!is.na(df_tpreel$SF)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_SF_lmm <-  subset(df_SF_lmm, select =-c(visit:LagDate, SF:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_SF_lmm )
row.names(df_SF_lmm) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_SF_lmm <- cbind(df_SF_lmm, fit.lme.SF$fitted[,1])
names(df_SF_lmm)[ncol(df_SF_lmm)] <- "SF_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees SF manquantes :
df_SF_lcmm <- df_tpreel[which(!is.na(df_tpreel$SF)),] # retire (si il y en a)
df_SF_lcmm <-  subset(df_SF_lcmm, select =-c(visit:LagDate, QL2:CF, FA:FI))
df_SF_lcmm <- df_SF_lcmm[order(df_SF_lcmm$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_SF_lcmm)
row.names(df_SF_lcmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(SF_lin_lin, newdata = df_SF_lcmm)
df_SF_lcmm <- cbind(df_SF_lcmm, pred_marg_transform$pred)
names(df_SF_lcmm)[ncol(df_SF_lcmm)] <- "SF_lcmm.fitted"

pSF_lcmm <- ggplot(data = df_SF_lmm, aes(x=Time, SF_lmm.fitted, 
                                            group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_SF_lcmm,aes(x=Time, y=SF_lcmm.fitted, 
                                     group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(45,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction sociale :
  Trajectoires linéaires (LMM) : pleins
  Trajectoires linéaires et lien linéaires (LCMM) : pointillés")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pSF_lcmm)

### Test des hypothèses #######################################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.SF)
resid_transform <- residuals(SF_lin_lin)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LMM): \n Fonction sociale")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LCMM): \n Fonction sociale")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction sociale: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE,
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction sociale: LCMM")


#Résidu prédit vs modèle prédit
plot(fit.lme.SF, main="Graphique des résidus (LMM): \n Fonction sociale")
plot(SF_lin_lin, main="Graphique des résidus (LCMM): \n Fonction sociale")


#--Trajectoires lineaires + lien spline----------------------------------------
lien <- "splines"
SF_lin_spl <- lcmm(SF ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                   data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_SF_ls <- subset(df_tpreel, select =-c(visit:LagDate, QL2:CF, FA:FI))
#on retire les lignes avec des NA:
df_SF_ls <- df_SF_ls[which(!is.na(df_SF_ls$SF)),]
#ajout du score transformé:
df_SF_ls <- cbind(df_SF_ls, SF_lin_spl$pred$pred_m)
names(df_SF_ls)[names(df_SF_ls) == "SF_lin_spl$pred$pred_m"] <- "score_transf_SF_lin_spl"
df_SF_ls <- df_SF_ls[order(df_SF_ls$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
SF_ls_p <- ggplot(df_SF_ls, aes(x=Time, y=score_transf_SF_lin_spl, 
                                group = BRAS)) + geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires linéaires et lien spline (LCMM) :
          Fonction sociale")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(SF_ls_p)

#LMM
length(fit.lme.SF$fitted[,1]) #vérification du nombre de lignes : données 
#manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données SF manquantes:
df_SF_lmm <- df_tpreel[which(!is.na(df_tpreel$SF)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_SF_lmm <-  subset(df_SF_lmm, select =-c(visit:LagDate, SF:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_SF_lmm )
row.names(df_SF_lmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
df_SF_lmm <- cbind(df_SF_lmm, fit.lme.SF$fitted[,1])
names(df_SF_lmm)[ncol(df_SF_lmm)] <- "SF_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees SF manquantes :
df_SF_lcmm_ls <- df_tpreel[which(!is.na(df_tpreel$SF)),] # retire (si il y en a)
df_SF_lcmm_ls <-  subset(df_SF_lcmm_ls, select =-c(visit:LagDate, QL2:CF, FA:FI))
df_SF_lcmm_ls <- df_SF_lcmm_ls[order(df_SF_lcmm_ls$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_SF_lcmm_ls)
row.names(df_SF_lcmm_ls) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(SF_lin_spl, newdata = df_SF_lcmm_ls)
df_SF_lcmm_ls <- cbind(df_SF_lcmm_ls, pred_marg_transform$pred)
names(df_SF_lcmm_ls)[ncol(df_SF_lcmm_ls)] <- "SF_lcmm_ls.fitted"

pSF_lcmm_ls <- ggplot(data = df_SF_lmm, aes(x=Time, SF_lmm.fitted, 
                                             group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_SF_lcmm_ls,aes(x=Time, y=SF_lcmm_ls.fitted, 
                                      group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(70,95))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction sociale :
            Trajectoires linéaires (LMM) : pleins
          Trajectoires linéaires et lien spline (LCMM) : pointillés")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pSF_lcmm_ls)

### Test des hypothèses ###########################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.SF)
resid_transform <- residuals(SF_lin_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LMM): \n Fonction sociale")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LCMM): \n Fonction sociale")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction sociale: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction sociale: LCMM")



#Résidu prédit vs modèle prédit
plot(fit.lme.SF, main="Graphique des résidus (LMM): \n Fonction sociale")
plot(SF_lin_spl, main="Graphique des résidus (LCMM): \n Fonction sociale")


#---Trajectoires splines + lien spline-----------------------------------------
lien <- "splines"
#retirer les NA s'il y en a:
df_SF_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$SF)),]
SF_spl_spl <- lcmm(fixed = SF ~ ns(Time,3):BRAS, random = ~ 1 + Time,
                   subject = "npat",data=df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_SF_ss <- subset(df_tpreel, select =-c(visit:LagDate, QL2:CF, FA:FI))
#on retire les lignes avec des NA:
df_SF_ss <- df_SF_ss[which(!is.na(df_SF_ss$SF)),]
#ajout du score transformé:
df_SF_ss <- cbind(df_SF_ss, SF_spl_spl$pred$pred_m)
names(df_SF_ss)[names(df_SF_ss) == "SF_spl_spl$pred$pred_m"] <- "score_transf_SF_spl"
df_SF_ss <- df_SF_ss[order(df_SF_ss$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
SF_ss_p <- ggplot(df_SF_ss, aes(x=Time, y=score_transf_SF_spl, 
                                group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires splines et lien spline (LCMM) :
          Fonction sociale")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(SF_ss_p)

#LMM
#pour ne pas avoir les NA dans la colonne SF:
df_SF_spl <- df_tpreel[which(!is.na(df_tpreel$SF)),]
#ne garder que les colonnes qui sont utiles pour l'échelle
df_SF_spl <-  subset(df_SF_spl, select =-c(visit:LagDate, QL2:CF, FA:FI))
df_SF_spl <- df_SF_spl[order(df_SF_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_SF_spl)
row.names(df_SF_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
df_SF_spl <- cbind(df_SF_spl, fit.lme.ns.SF.2$fitted[,1])
names(df_SF_spl)[ncol(df_SF_spl)] <- "SF_spl.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees SF manquantes :
df_SF_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$SF)),] # retire (si il y en a)
df_SF_lcmm_spl <-  subset(df_SF_lcmm_spl, select =-c(visit:LagDate, QL2:CF, FA:FI))
df_SF_lcmm_spl <- df_SF_lcmm_spl[order(df_SF_lcmm_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_SF_lcmm_spl)
row.names(df_SF_lcmm_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(SF_spl_spl, newdata = df_SF_lcmm_spl)
df_SF_lcmm_spl <- cbind(df_SF_lcmm_spl, pred_marg_transform$pred)
names(df_SF_lcmm_spl)[ncol(df_SF_lcmm_spl)] <- "SF_lcmm_spl.fitted"

pSF_lcmm_spl <- ggplot(data = df_SF_spl, aes(x=Time, y=SF_spl.fitted, 
                                             group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_SF_lcmm_spl,aes(x=Time, y=SF_lcmm_spl.fitted, 
                                      group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(45,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction sociale :
  Trajectoires splines (LMM) : pleins
  Trajectoires splines et lien spline (LCMM) : pointillés") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pSF_lcmm_spl)

### Test des hypothèses ###########################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.ns.SF.2)
resid_transform <- residuals(SF_spl_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LMM): \n Fonction sociale")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LCMM): \n Fonction sociale")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction sociale: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction sociale: LCMM")


#Résidu prédit vs modèle prédit
plot(fit.lme.ns.SF.2, main="Graphique des résidus (LMM): \n Fonction sociale")
plot(SF_spl_spl, main="Graphique des résidus (LCMM): \n Fonction sociale")


###############################################################################
#                 PLOT CURVILINEARITE
###############################################################################

col <- rainbow(3)
plot(SF_lin_lin, which = "linkfunction", bty = "l", ylab = "Score prédit (%)", 
     xlab = "Processus latent sous-jacent", lwd = 2, col = col[1],
     main = "Fonction de lien estimée \n pour la fonction sociale")
plot(SF_lin_spl, which = "linkfunction", add = TRUE, col = col[2], lwd = 2)
plot(SF_spl_spl, which = "linkfunction", add = TRUE, col = col[3], lwd = 2)
legend(x = "topleft",
       legend = c("trajectoire linéaire et\n fonction de lien linéaire", 
                  "trajectoire linéaire et\n fonction de lien spline", 
                  "trajectoire spline et\n fonction de lien spline"), lty = 1, 
       cex=0.8, col = col, bty = "n", lwd = 2)

