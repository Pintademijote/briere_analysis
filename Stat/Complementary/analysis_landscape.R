setwd("C:/Users/pglem/Documents/Master/Stage M1/Rendu/Data/Data")
#setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/Data")

require("enaR")
require("igraph")
require("NetIndices")
require("nlme")
require("mgcv")
library("ggplot2")
library("MASS")
library("bestglm")

# Import trophic link database
tab_trophic=read.delim("tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
# Import taxon occurrence data
tab_releve=read.delim("tab_releve.txt")

sites=levels(tab_releve$Id_Sites)
for (i in sites) {
  x=tab_releve[tab_releve$Id_Sites==i,]
  if(all(c("Aeshnidae","Anisoptere","Libellulidae") %in% x$Taxon)){
    tab_releve[tab_releve$Id_Sites==i & tab_releve$Taxon=="Anisoptere",7]="Libellulidae"
  }else{
    if(all(c("Aeshnidae","Anisoptere") %in% x$Taxon)){
      tab_releve[tab_releve$Id_Sites==i & tab_releve$Taxon=="Aeshnidae",7]="Anisoptere"
    }else{
      if(all(c("Libellulidae","Anisoptere") %in% x$Taxon)){
        tab_releve[tab_releve$Id_Sites==i & tab_releve$Taxon=="Libellulidae",7]="Anisoptere"
      }
    }
    
    
  }
}

levels(tab_releve$Taxon)[which(!(levels(tab_releve$Taxon) %in% colnames(Matrix_Reg_Manual)))]

tab_releve=tab_releve[tab_releve$Taxon!="Couleuvre a collier",]

# Import environmental variables
envi_data=read.table("envi_data.txt",header=T)
# Import distance from marsh
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
# Import manual matrix
Matrix_Reg_Manual=read.delim2("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
data=read.delim("data.txt")

# Import function to create matrix
#source("S:/2018/pglemasle/Data/Script/matrix_function.R")
source("C:/Users/pglem/Documents/Master/Stage M1/Rendu/Data/Script_Function/matrix_function.R")
##########Marais/Mare

releve_marais=tab_releve[tab_releve$Milieux=="Marais",]
releve_marais=droplevels(releve_marais)
marais=levels(releve_marais$Id_Sites)

##########

data_format(tab_releve,envi_data,dist)

#create_matrix(tab_trophic)

#create_matrix_marre(tab_releve,Matrix_tab_trophic,envi_data,Count = F, year=F,name_tab = "data")

create_matrix_marre(tab_releve,Matrix_Reg_Manual,envi_data,Count = T, year=F,name_tab = "data")


###########

glmN=datacount[,c("n","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmN)[1]="y"

my_QP <- function(link="log") {
  QP <- quasipoisson(link=link)
  QP$aic <- function (y, n, mu, wt, dev) {
    nobs <- length(y)
    disp <- dev/n
    dev0 <- -2*sum((-mu + y*log(mu) -lfactorial(y))*wt/disp)
    dev0 + 2  ## extra penalty parameter
  }
  QP
}


choiceN=bestglm(Xy = glmN,
        family = my_QP,
        IC = "AIC",                 # Information criteria for
        method = "exhaustive")

choiceN$BestModels


summary(choiceN$BestModel)

write.table(summary(choiceN$BestModel)$coefficients,"clipboard",sep="\t")

###

glmC=datacount[,c("C","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmC)[1]="y"

glmC=glmC[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

choiceC=bestglm(Xy = glmC,
                family=gaussian,
                IC = "AIC",                 # Information criteria for
                method = "exhaustive")

choiceC$BestModels

summary(choiceC$BestModel)

write.table(summary(choiceC$BestModel)$coefficients,"clipboard",sep="\t")
###

glmTL=datacount[,c("Mean_trophic_level","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmTL)[1]="y"

glmTL=glmTL[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

glmTL=bestglm(Xy = glmTL,
                family=gaussian,
                IC = "AIC",                 # Information criteria for
                method = "exhaustive")

glmTL$BestModels

summary(glmTL$BestModel)


write.table(summary(glmTL$BestModel)$coefficients,"clipboard",sep="\t")


###

glmOI=datacount[,c("Mean_Omnivory_Index","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmOI)[1]="y"

glmOI=glmOI[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

glmOI=bestglm(Xy = glmOI,
              family=gaussian,
              IC = "AIC",                 # Information criteria for
              method = "exhaustive")

glmOI$BestModels

summary(glmOI$BestModel)

write.table(summary(glmOI$BestModel)$coefficients,"clipboard",sep="\t")

###


glmpB=datacount[,c("Percentage_basal","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmpB)[1]="y"

glmpB=glmpB[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

glmpB=bestglm(Xy = glmpB,
              family=gaussian,
              IC = "AIC",                 # Information criteria for
              method = "exhaustive")

glmpB$BestModels

summary(glmpB$BestModel)


write.table(summary(glmpB$BestModel)$coefficients,"clipboard",sep="\t")

###

glmpT=datacount[,c("Percentage_top","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmpT)[1]="y"

glmpT=glmpT[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

glmpT=bestglm(Xy = glmpT,
              family=gaussian,
              IC = "AIC",                 # Information criteria for
              method = "exhaustive")

glmpT$BestModels

summary(glmpT$BestModel)

write.table(summary(glmpT$BestModel)$coefficients,"clipboard",sep="\t")

###



glmpI=datacount[,c("Percentage_intermediate","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmpI)[1]="y"

glmpI=glmpI[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

glmpI=bestglm(Xy = glmpI,
              family=gaussian,
              IC = "AIC",                 # Information criteria for
              method = "exhaustive")

glmpI$BestModels

summary(glmpI$BestModel)

write.table(summary(glmpI$BestModel)$coefficients,"clipboard",sep="\t")


####

glmmodu=datacount[,c("Modularity","Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies")]

colnames(glmmodu)[1]="y"

glmmodu=glmmodu[,c("Bois","Cultures","Prairies","Surf_Urbanisees","Longueur_haies","y")]

glmmodu=bestglm(Xy = glmmodu,
              family=gaussian,
              IC = "AIC",                 # Information criteria for
              method = "exhaustive")

glmmodu$BestModels

summary(glmmodu$BestModel)

write.table(summary(glmmodu$BestModel)$coefficients,"clipboard",sep="\t")



