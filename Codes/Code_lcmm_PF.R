###############################################################################
#                           LCMM                                              #
###############################################################################

#------------------------------------------------------------------------------
#                 DIMENSION PF
#------------------------------------------------------------------------------

#---Trajectoires lin�aires + lien lin�aire-------------------------------------
lien <- "linear" #choix de la fonction de lien pour la transformation
PF_lin_lin <- lcmm(PF2 ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                    data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous int�ressent
df_PF_ll <- subset(df_tpreel, select =-c(visit:LagDate, QL2, RF2:FI))
#on retire les lignes avec des NA:
df_PF_ll <- df_PF_ll[which(!is.na(df_PF_ll$PF2)),]
#ajout du score transform�:
df_PF_ll <- cbind(df_PF_ll, PF_lin_lin$pred$pred_m)
names(df_PF_ll)[names(df_PF_ll) == "PF_lin_lin$pred$pred_m"] <- "score_transf_PF_lin"
df_PF_ll <- df_PF_ll[order(df_PF_ll$Time), ] #ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'�chelle du lcmm
PF_ll_p <- ggplot(df_PF_ll, 
                   aes(x=Time, y=score_transf_PF_lin, group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transform�", limits = c(-2,10))+
  scale_color_manual(labels = c("Contr�le", "APAD"), values = c("red", "blue"))+
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires lin�aires et lien lin�aire (LCMM) :
          Fonction physique")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(PF_ll_p)

#LMM
length(fit.lme.PF$fitted[,1]) #v�rification du nombre de lignes : donn�es 
#manquantes non prises en compte
#on construit un fichier o� on supprime les lignes avec donn�es PF manquantes:
df_PF_lmm <- df_tpreel[which(!is.na(df_tpreel$PF2)),] # retire (si il y en a)
#on retire les colonnes inutiles � l'�chelle:
df_PF_lmm <-  subset(df_PF_lmm, select =-c(visit:LagDate, QL2, RF2:FI))
#v�rification des dimensions des tableaux:
dim(df_tpreel)
dim(df_PF_lmm)
row.names(df_PF_lmm) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
df_PF_lmm <- cbind(df_PF_lmm, fit.lme.PF$fitted[,1])
names(df_PF_lmm)[ncol(df_PF_lmm)] <- "PF_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees PF manquantes :
df_PF_lcmm <- df_tpreel[which(!is.na(df_tpreel$PF2)),] # retire (si il y en a)
df_PF_lcmm <-  subset(df_PF_lcmm, select =-c(visit:LagDate, QL2, RF2:FI))
df_PF_lcmm <- df_PF_lcmm[order(df_PF_lcmm$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_PF_lcmm)
row.names(df_PF_lcmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(PF_lin_lin, newdata = df_PF_lcmm)
df_PF_lcmm <- cbind(df_PF_lcmm, pred_marg_transform$pred)
names(df_PF_lcmm)[ncol(df_PF_lcmm)] <- "PF_lcmm.fitted"

pPF_lcmm <- ggplot(data = df_PF_lmm, aes(x=Time, PF_lmm.fitted, 
                                           group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_PF_lcmm,aes(x=Time, y=PF_lcmm.fitted, 
                                   group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen pr�dit (%)", limits = c(75, 100))+
  scale_color_manual(labels = c("Contr�le", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction physique :
  Trajectoires lin�aires (LMM) : pleins
  Trajectoires lin�aires et lien lin�aires (LCMM) : pointill�s")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pPF_lcmm)

### Test des hypoth�ses #######################################################

#R�cup�ration des r�sidus
resid_initial <- residuals(fit.lme.PF)
resid_transform <- residuals(PF_lin_lin)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des r�sidus (LMM): \n Fonction physique")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des r�sidus (LCMM): \n Fonction physique")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "R�sidus standardis�s", 
       main = "QQ-plot Fonction physique: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE,
       xlab = "Quantiles de la loi normale",
       ylab = "R�sidus standardis�s", 
       main = "QQ-plot Fonction physique: LCMM")


#R�sidu pr�dit vs mod�le pr�dit
plot(fit.lme.PF,
     main="Graphique des r�sidus (LMM): \n Fonction physique")
plot(PF_lin_lin,
     main="Graphique des r�sidus (LCMM): \n Fonction physique")


#--Trajectoires lineaires + lien spline----------------------------------------
lien <- "splines"
PF_lin_spl <- lcmm(PF2 ~ Time + Time:BRAS, random = ~ Time, subject = "npat",
                    data = df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous int�ressent
df_PF_ls <- subset(df_tpreel, select =-c(visit:LagDate, QL2, RF2:FI))
#on retire les lignes avec des NA:
df_PF_ls <- df_PF_ls[which(!is.na(df_PF_ls$PF2)),]
#ajout du score transform�:
df_PF_ls <- cbind(df_PF_ls, PF_lin_spl$pred$pred_m)
names(df_PF_ls)[names(df_PF_ls) == "PF_lin_spl$pred$pred_m"] <- "score_transf_PF_lin_spl"
df_PF_ls <- df_PF_ls[order(df_PF_ls$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'�chelle du lcmm
PF_ls_p <- ggplot(df_PF_ls, aes(x=Time, y=score_transf_PF_lin_spl, 
                                  group = BRAS)) + geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transform�", limits = c(-2,10))+
  scale_color_manual(labels = c("Contr�le", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires lin�aires et lien spline (LCMM) :
          Fonction physique")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(PF_ls_p)

#LMM
length(fit.lme.PF$fitted[,1]) #v�rification du nombre de lignes : donn�es 
#manquantes non prises en compte
#on construit un fichier o� on supprime les lignes avec donn�es PF manquantes:
df_PF_lmm <- df_tpreel[which(!is.na(df_tpreel$PF2)),] # retire (si il y en a)
#on retire les colonnes inutiles � l'�chelle:
df_PF_lmm <-  subset(df_PF_lmm, select =-c(visit:LagDate, QL2, RF2:FI))
#v�rification des dimensions des tableaux:
dim(df_tpreel)
dim(df_PF_lmm )
row.names(df_PF_lmm) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
df_PF_lmm <- cbind(df_PF_lmm, fit.lme.PF$fitted[,1])
names(df_PF_lmm)[ncol(df_PF_lmm)] <- "PF_lmm.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees PF manquantes :
df_PF_lcmm_ls <- df_tpreel[which(!is.na(df_tpreel$PF)),] # retire (si il y en a)
df_PF_lcmm_ls <-  subset(df_PF_lcmm_ls, select =-c(visit:LagDate, QL2, RF2:FI))
df_PF_lcmm_ls <- df_PF_lcmm_ls[order(df_PF_lcmm_ls$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_PF_lcmm_ls)
row.names(df_PF_lcmm_ls) <- NULL
#-- Valeurs predites ----
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(PF_lin_spl, newdata = df_PF_lcmm_ls)
df_PF_lcmm_ls <- cbind(df_PF_lcmm_ls, pred_marg_transform$pred)
names(df_PF_lcmm_ls)[ncol(df_PF_lcmm_ls)] <- "PF_lcmm_ls.fitted"

pPF_lcmm_ls <- ggplot(data = df_PF_lmm, aes(x=Time, PF_lmm.fitted, 
                                              group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_PF_lcmm_ls,aes(x=Time, y=PF_lcmm_ls.fitted, 
                                      group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen pr�dit (%)", limits = c(75,90))+
  scale_color_manual(labels = c("Contr�le", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction physique :
          Trajectoires lin�aires (LMM) : pleins
          Trajectoires lin�aires et lien spline (LCMM) : pointill�s")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pPF_lcmm_ls)

### Test des hypoth�ses ###########################################

#R�cup�ration des r�sidus
resid_initial <- residuals(fit.lme.PF)
resid_transform <- residuals(PF_lin_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40, 
     main="Histogramme des r�sidus (LMM): \n Fonction physique")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des r�sidus (LCMM): \n Fonction physique")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "R�sidus standardis�s", 
       main = "QQ-plot Fonction physique: LMM")

qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "R�sidus standardis�s", 
       main = "QQ-plot Fonction physique: LCMM")

#R�sidu pr�dit vs mod�le pr�dit
plot(fit.lme.PF,
                    main="Graphique des r�sidus (LMM): \n Fonction physique")
plot(PF_lin_spl,
                     main="Graphique des r�sidus (LCMM): \n Fonction physique")


#---Trajectoires splines + lien spline-----------------------------------------
lien <- "splines"
#retirer les NA s'il y en a:
df_PF_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$PF2)),]
PF_spl_spl <- lcmm(fixed = PF2 ~ ns(Time,3):BRAS, random = ~ 1 + Time,
                    subject = "npat",data=df_tpreel, link = lien)
#on ne garde que les colonnes du tableau initial qui nous int�ressent
df_PF_ss <- subset(df_tpreel, select =-c(visit:LagDate, QL2, RF2:FI))
#on retire les lignes avec des NA:
df_PF_ss <- df_PF_ss[which(!is.na(df_PF_ss$PF2)),]
#ajout du score transform�:
df_PF_ss <- cbind(df_PF_ss, PF_spl_spl$pred$pred_m)
names(df_PF_ss)[names(df_PF_ss) == "PF_spl_spl$pred$pred_m"] <- "score_transf_PF_spl"
df_PF_ss <- df_PF_ss[order(df_PF_ss$Time), ] # ordonner par les temps

### Partie Graphique ##########################################################
#Graphe dans l'�chelle du lcmm
PF_ss_p <- ggplot(df_PF_ss, aes(x=Time, y=score_transf_PF_spl, 
                                  group = BRAS)) + 
  geom_line(aes(color=BRAS)) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score transform�", limits = c(-2,10))+
  scale_color_manual(labels = c("Contr�le", "APAD"), values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Trajectoires splines et lien spline (LCMM) :
          Fonction physique")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(PF_ss_p)

#LMM
#pour ne pas avoir les NA dans la colonne PF:
df_PF_spl <- df_tpreel[which(!is.na(df_tpreel$PF2)),]
#ne garder que les colonnes qui sont utiles pour l'�chelle
df_PF_spl <-  subset(df_PF_spl, select =-c(visit:LagDate, QL2, RF2:FI))
df_PF_spl <- df_PF_spl[order(df_PF_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_PF_spl)
row.names(df_PF_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
df_PF_spl <- cbind(df_PF_spl, fit.lme.ns.PF.2$fitted[,1])
names(df_PF_spl)[ncol(df_PF_spl)] <- "PF_spl.fitted"

#LCMM
# on construit un fichier ou on supprime les lignes avec donnees PF manquantes :
df_PF_lcmm_spl <- df_tpreel[which(!is.na(df_tpreel$PF2)),] # retire (si il y en a)
df_PF_lcmm_spl <-  subset(df_PF_lcmm_spl, select =-c(visit:LagDate, QL2, RF2:FI))
df_PF_lcmm_spl <- df_PF_lcmm_spl[order(df_PF_lcmm_spl$Time), ] # ordonner par les temps
dim(df_tpreel)
dim(df_PF_lcmm_spl)
row.names(df_PF_lcmm_spl) <- NULL
#-- Valeurs predites ----------------------------------------------------------
# on ajoute les valeurs predites au tableau de donnees :
pred_marg_transform <- predictY(PF_spl_spl, newdata = df_PF_lcmm_spl)
df_PF_lcmm_spl <- cbind(df_PF_lcmm_spl, pred_marg_transform$pred)
names(df_PF_lcmm_spl)[ncol(df_PF_lcmm_spl)] <- "PF_lcmm_spl.fitted"

pPF_lcmm_spl <- ggplot(data = df_PF_spl, aes(x=Time, y=PF_spl.fitted, 
                                               group = BRAS)) +
  geom_line(aes(color=BRAS))+
  geom_line(data = df_PF_lcmm_spl,aes(x=Time, y=PF_lcmm_spl.fitted, 
                                       group = BRAS, color=BRAS),
            linetype = "longdash", size = 1) +
  scale_x_continuous(name = "Temps (mois)", limits = c(0,30))+
  scale_y_continuous(name = "Score moyen pr�dit (%)", limits = c(65,100))+
  scale_color_manual(labels = c("Contr�le", "APAD"),
                     values = c("red", "blue")) +
  theme(legend.position="right") + labs(color="Bras")+ 
  ggtitle("Fonction physique :
  Trajectoires splines (LMM) : pleins
  Trajectoires splines et lien spline (LCMM) : pointill�s") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8))
plot(pPF_lcmm_spl)

### Test des hypoth�ses ###########################################

#R�cup�ration des r�sidus
resid_initial <- residuals(fit.lme.ns.PF.2)
resid_transform <- residuals(PF_spl_spl)

#Shapiro-Wilk Test
shapiro.test(resid_initial)
shapiro.test(resid_transform)

#Histogramme
par(mfrow=c(1,2))
hist(resid_initial, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des r�sidus (LMM): \n Fonction physique")
Range = seq(min(resid_initial), max(resid_initial), 
            length = length(resid_initial))
Norm = dnorm(Range, mean = mean(resid_initial), sd = sd(resid_initial))
lines(Range, Norm, col = "red", lwd = 2)

hist(resid_transform, prob=TRUE, col="darkgray", breaks=40,
     main="Histogramme des r�sidus (LCMM): \n Fonction physique")
Range = seq(min(resid_transform), max(resid_transform), 
            length = length(resid_transform))
Norm = dnorm(Range, mean = mean(resid_transform), sd = sd(resid_transform))
lines(Range, Norm, col = "red", lwd = 2)

#QQ-Plot
qqPlot(resid_initial, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "R�sidus standardis�s", 
       main = "QQ-plot Fonction physique: LMM")
qqPlot(resid_transform, envelope = .95, id=FALSE, 
       xlab = "Quantiles de la loi normale",
       ylab = "R�sidus standardis�s", 
       main = "QQ-plot Fonction physique: LCMM")


#R�sidu pr�dit vs mod�le pr�dit
plot(fit.lme.ns.PF.2,
     main="Graphique des r�sidus (LMM): \n Fonction physique")
plot(PF_spl_spl, main="Graphique des r�sidus (LCMM): \n Fonction physique")


###############################################################################
#                 PLOT CURVILINEARITE
###############################################################################

col <- rainbow(3)
plot(PF_lin_lin, which = "linkfunction", bty = "l", ylab = "Score pr�dit (%)", 
     xlab = "Processus latent sous-jacent", lwd = 2, col = col[1],
     main = "Fonction de lien estim�e \n pour la fonction physique")
plot(PF_lin_spl, which = "linkfunction", add = TRUE, col = col[2], lwd = 2)
plot(PF_spl_spl, which = "linkfunction", add = TRUE, col = col[3], lwd = 2)
legend(x = "topleft",
       legend = c("trajectoire lin�aire et\n fonction de lien lin�aire", 
                  "trajectoire lin�aire et\n fonction de lien spline", 
                  "trajectoire spline et\n fonction de lien spline"), lty = 1, 
       cex=0.8, col = col, bty = "n", lwd = 2)
