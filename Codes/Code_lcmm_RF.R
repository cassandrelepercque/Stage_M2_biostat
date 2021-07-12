###############################################################################
#                           LCMM                                              #
###############################################################################

#------------------------------------------------------------------------------
#                 DIMENSION RF
#------------------------------------------------------------------------------

#---Trajectoires linéaires + lien linéaire-------------------------------------
lien <- "linear" #choix de la fonction de lien pour la transformation
RF_lin_lin <- lcmm(RF2 ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                   data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_RF_ll <- subset(df_tpreel, select =-c(visit:LagDate, QL2:PF2, EF:FI))
#on retire les lignes avec des NA:
df_RF_ll <- df_RF_ll[which(!is.na(df_RF_ll$RF2)),]
#ajout du score transformé:
df_RF_ll <- cbind(df_RF_ll, RF_lin_lin$pred$pred_m)
names(df_RF_ll)[names(df_RF_ll) == "RF_lin_lin$pred$pred_m"] <- "score_transf_RF_lin"
df_RF_ll <- df_RF_ll[order(df_RF_ll$Time), ] #ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
RF_ll_p <- ggplot(df_RF_ll, 
                  aes(x=Time, y=score_transf_RF_lin, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue"))+
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires linéaires et lien linéaire (LCMM) :
          Fonction personnelle")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(RF_ll_p)

#LMM
length(fit.lme.RF$fitted[,1]) #vérification du nombre de lignes : données 
#manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données RF manquantes:
df_RF_lmm <- df_tpreel[which(!is.na(df_tpreel$RF2)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_RF_lmm <-  subset(df_RF_lmm, select =-c(visit:LagDate, QL2:PF2, EF:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_RF_lmm)
row.names(df_RF_lmm) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_RF_lmm <- cbind(df_RF_lmm, fit.lme.RF$fitted[,1])
names(df_RF_lmm)[ncol(df_RF_lmm)] <- "RF_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees RF manquantes :
df_RF_lcmm <- df_tpreel[which(!is.na(df_tpreel$RF2)),] # retire (si il y en a)
df_RF_lcmm <-  subset(df_RF_lcmm, select =-c(visit:LagDate, QL2:PF2, EF:FI))
df_RF_lcmm <- df_RF_lcmm[order(df_RF_lcmm$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_RF_lcmm)
row.names(df_RF_lcmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(RF_lin_lin, newdata = df_RF_lcmm)
df_RF_lcmm <- cbind(df_RF_lcmm, pred_marg_transform$pred)
names(df_RF_lcmm)[ncol(df_RF_lcmm)] <- "RF_lcmm.fitted"

pRF_lcmm <- ggplot(data = df_RF_lmm, aes(x=Time, RF_lmm.fitted, 
                                         group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_RF_lcmm,aes(x=Time, y=RF_lcmm.fitted, 
                                  group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(75, 100))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction personnelle :
  Trajectoires linéaires (LMM) : pleins
  Trajectoires linéaires et lien linéaires (LCMM) : pointillés")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pRF_lcmm)

### Test des hypothèses #######################################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.RF)
resid_transform <- residuals(RF_lin_lin)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LMM): \n Fonction personnelle")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LCMM): \n Fonction personnelle")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction personnelle: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE,
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction personnelle: LCMM")


#Résidu prédit vs modèle prédit
plot(fit.lme.RF, main="Graphique des résidus (LMM): \n Fonction personnelle")
plot(RF_lin_lin, main="Graphique des résidus (LCMM):\n Fonction personnelle")


#--Trajectoires lineaires + lien spline----------------------------------------
lien <- "splines"
RF_lin_spl <- lcmm(RF2 ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                   data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_RF_ls <- subset(df_tpreel, select =-c(visit:LagDate, QL2:PF2, EF:FI))
#on retire les lignes avec des NA:
df_RF_ls <- df_RF_ls[which(!is.na(df_RF_ls$RF2)),]
#ajout du score transformé:
df_RF_ls <- cbind(df_RF_ls, RF_lin_spl$pred$pred_m)
names(df_RF_ls)[names(df_RF_ls) == "RF_lin_spl$pred$pred_m"] <- "score_transf_RF_lin_spl"
df_RF_ls <- df_RF_ls[order(df_RF_ls$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
RF_ls_p <- ggplot(df_RF_ls, aes(x=Time, y=score_transf_RF_lin_spl, 
                                group = BRAS)) + geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires linéaires et lien spline (LCMM) :
          Fonction personnelle")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(RF_ls_p)

#LMM
length(fit.lme.RF$fitted[,1]) #vérification du nombre de lignes : données 
#manquantes non prises en compte
#on construit un fichier où on supprime les lignes avec données RF manquantes:
df_RF_lmm <- df_tpreel[which(!is.na(df_tpreel$RF2)),] # retire (si il y en a)
#on retire les colonnes inutiles à l'échelle:
df_RF_lmm <-  subset(df_RF_lmm, select =-c(visit:LagDate, QL2:PF2, EF:FI))
#vérification des dimensions des tableaux:
dim(df_tpreel)
dim(df_RF_lmm )
row.names(df_RF_lmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
df_RF_lmm <- cbind(df_RF_lmm, fit.lme.RF$fitted[,1])
names(df_RF_lmm)[ncol(df_RF_lmm)] <- "RF_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees RF manquantes :
df_RF_lcmm_ls <- df_tpreel[which(!is.na(df_tpreel$RF)),] # retire (si il y en a)
df_RF_lcmm_ls <-  subset(df_RF_lcmm_ls, select =-c(visit:LagDate, QL2:PF2, EF:FI))
df_RF_lcmm_ls <- df_RF_lcmm_ls[order(df_RF_lcmm_ls$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_RF_lcmm_ls)
row.names(df_RF_lcmm_ls) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(RF_lin_spl, newdata = df_RF_lcmm_ls)
df_RF_lcmm_ls <- cbind(df_RF_lcmm_ls, pred_marg_transform$pred)
names(df_RF_lcmm_ls)[ncol(df_RF_lcmm_ls)] <- "RF_lcmm_ls.fitted"

pRF_lcmm_ls <- ggplot(data = df_RF_lmm, aes(x=Time, RF_lmm.fitted, 
                                            group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_RF_lcmm_ls,aes(x=Time, y=RF_lcmm_ls.fitted, 
                                     group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(80,97))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction personnelle :
            Trajectoires linéaires (LMM) : pleins
          Trajectoires linéaires et lien spline (LCMM) : pointillés")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pRF_lcmm_ls)

### Test des hypothèses ###########################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.RF)
resid_transform <- residuals(RF_lin_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des résidus (LMM): \n Fonction personnelle")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LCMM): \n Fonction personnelle")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction personnelle: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction personnelle: LCMM")


#Résidu prédit vs modèle prédit
plot(fit.lme.RF, main="Graphique des résidus (LMM): \n Fonction personnelle")
plot(RF_lin_spl, main="Graphique des résidus (LCMM): \n Fonction personnelle")


#---Trajectoires splines + lien spline-----------------------------------------
lien <- "splines"
#retirer les NA s'il y en a:
df_RF_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$RF2)),]
RF_spl_spl <- lcmm(fixed = RF2 ~ ns(Time,3):BRAS, random = ~ 1 + Time,
                   subject = "npat",data=df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous intéressent
df_RF_ss <- subset(df_tpreel, select =-c(visit:LagDate, QL2:PF2, EF:FI))
#on retire les lignes avec des NA:
df_RF_ss <- df_RF_ss[which(!is.na(df_RF_ss$RF2)),]
#ajout du score transformé:
df_RF_ss <- cbind(df_RF_ss, RF_spl_spl$pred$pred_m)
names(df_RF_ss)[names(df_RF_ss) == "RF_spl_spl$pred$pred_m"] <- "score_transf_RF_spl"
df_RF_ss <- df_RF_ss[order(df_RF_ss$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'échelle du lcmm
RF_ss_p <- ggplot(df_RF_ss, aes(x=Time, y=score_transf_RF_spl, 
                                group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transformé", limits = c(-2,10))+
  scale_color_manual(labels = c("Contrôle", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires splines et lien spline (LCMM) :
          Fonction personnelle")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(RF_ss_p)

#LMM
#pour ne pas avoir les NA dans la colonne RF:
df_RF_spl <- df_tpreel[which(!is.na(df_tpreel$RF2)),]
#ne garder que les colonnes qui sont utiles pour l'échelle
df_RF_spl <-  subset(df_RF_spl, select =-c(visit:LagDate, QL2:PF2, EF:FI))
df_RF_spl <- df_RF_spl[order(df_RF_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_RF_spl)
row.names(df_RF_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
df_RF_spl <- cbind(df_RF_spl, fit.lme.ns.RF.2$fitted[,1])
names(df_RF_spl)[ncol(df_RF_spl)] <- "RF_spl.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees RF manquantes :
df_RF_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$RF2)),] # retire (si il y en a)
df_RF_lcmm_spl <-  subset(df_RF_lcmm_spl, select =-c(visit:LagDate, QL2:PF2, EF:FI))
df_RF_lcmm_spl <- df_RF_lcmm_spl[order(df_RF_lcmm_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_RF_lcmm_spl)
row.names(df_RF_lcmm_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(RF_spl_spl, newdata = df_RF_lcmm_spl)
df_RF_lcmm_spl <- cbind(df_RF_lcmm_spl, pred_marg_transform$pred)
names(df_RF_lcmm_spl)[ncol(df_RF_lcmm_spl)] <- "RF_lcmm_spl.fitted"

pRF_lcmm_spl <- ggplot(data = df_RF_spl, aes(x=Time, y=RF_spl.fitted, 
                                             group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_RF_lcmm_spl,aes(x=Time, y=RF_lcmm_spl.fitted, 
                                      group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen prédit (%)", limits = c(65,100))+
  scale_color_manual(labels = c("Contrôle", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction personnelle :
  Trajectoires splines (LMM) : pleins
  Trajectoires splines et lien spline (LCMM) : pointillés") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pRF_lcmm_spl)

### Test des hypothèses ###########################################

#Récupération des résidus
resid_initial <- residuals(fit.lme.ns.RF.2)
resid_transform <- residuals(RF_spl_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LMM): \n Fonction personnelle")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des résidus (LCMM): \n Fonction personnelle")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction personnelle: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "Résidus standardisés", 
       main = "QQ-plot Fonction personnelle: LCMM")


#Résidu prédit vs modèle prédit
plot(fit.lme.ns.RF.2, 
     main="Graphique des résidus (LMM): \n Fonction personnelle")
plot(RF_spl_spl, main="Graphique des résidus (LCMM): \n Fonction personnelle")


###############################################################################
#                 PLOT CURVILINEARITE
###############################################################################

col <- rainbow(3)
plot(RF_lin_lin, which = "linkfunction", bty = "l", ylab = "Score prédit (%)", 
     xlab = "Processus latent sous-jacent", lwd = 2, col = col[1],
     main = "Fonction de lien estimée \n pour la fonction personnelle")
plot(RF_lin_spl, which = "linkfunction", add = TRUE, col = col[2], lwd = 2)
plot(RF_spl_spl, which = "linkfunction", add = TRUE, col = col[3], lwd = 2)
legend(x = "topleft",
       legend = c("trajectoire linéaire et\n fonction de lien linéaire", 
                  "trajectoire linéaire et\n fonction de lien spline", 
                  "trajectoire spline et\n fonction de lien spline"), lty = 1, 
       cex=0.8, col = col, bty = "n", lwd = 2)
