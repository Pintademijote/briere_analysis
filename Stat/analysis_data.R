setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/PG_Sept2018/Data/Data")

require("enaR")
require("igraph")
require("NetIndices")
require("nlme")
require("mgcv")
library("ggplot2")

# Import trophic link database
tab_trophic=read.delim("tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
# Import taxon occurrence data
tab_releve=read.delim("tab_releve.txt")
# Import manual matrix
Matrix_Reg_Manual=read.delim2("Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)

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

#data=read.delim("data.txt")

# Import function to create matrix
#source("S:/2018/pglemasle/Data/Script/matrix_function.R")
source("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/PG_Sept2018/Data/Script_Function/matrix_function.R")
##########Marais/Mare

releve_marais=tab_releve[tab_releve$Milieux=="Marais",]
releve_marais=droplevels(releve_marais)
marais=levels(releve_marais$Id_Sites)

##########

data_format(tab_releve,envi_data,dist)

#create_matrix(tab_trophic)

#create_matrix_marre(tab_releve,Matrix_tab_trophic,envi_data,Count = F, year=F,name_tab = "data")

create_matrix_marre(tab_releve,Matrix_Reg_Manual,envi_data,Count = T, year=F,name_tab = "data")


#data=data[!is.na(data$Hydroperiode),]
write.table(datacount,"data_observed.txt", sep="\t")
data=datacount

#########

div=c(1:length(marais))
for (i in 1:length(marais)) {
  
  div[i]=length(levels(droplevels(releve_marais[releve_marais$Id_Sites==marais[i],7])))
}

mean(div)

mean(data$n)

meanresults=matrix(NA,nrow=length(div)+length(data$n),ncol=2)
meanresults=as.data.frame(meanresults)
meanresults[1:155,1]=data$n
meanresults[1:155,2]="Ponds"
meanresults[156:173,1]=div
meanresults[156:173,2]="Marsh"
colnames(meanresults)=c("Richesse_spécifique","Lieu")

wilcox.test(meanresults[1:155,1],meanresults[156:173,1])

########################## RELATION AMONG DESCRIPTORS
# Environmental
pairs(data[,c(13:15,17:20,25)])
# Topological
pairs(data[,c(2:7)])

cor_envi=cor(data[,c(13:15,17:20,25)])
cor_envi>=0.4 | cor_envi<=-0.4

cor_topo=cor(data[,c(2:7)])
cor_topo>=0.4 | cor_topo<=-0.4

#########Connectance

m0 <- lm(C ~  scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
         , data = data)

anova(m0)

summary(m0)

write.table(summary(m0)$coefficients,"clipboard",sep="\t")

#########SPecies Richness

m1=glm(n ~ scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data, family = quasipoisson(link = "log"))

anova(m1, test="Chisq")

summary(m1)

write.table(summary(m1)$coefficients,"clipboard",sep="\t")

#########Modularity
m2=lm(Modularity ~  scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data)

anova(m2, test="Chisq")

summary(m2)

write.table(summary(m2)$coefficients,"clipboard",sep="\t")

m2=gam(n ~ s(scale(log(Surface)))+ Distance
       +  scale(log(Profondeur_maximum))
       + Vegetation_aquatique 
       + Hydroperiode + Anthropique
       , data = data)


summary(m2)

plot(m2)

#########Mean Trophic Level

m3=lm(Mean_trophic_level ~  scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data)

anova(m3, test="Chisq")

summary(m3)

write.table(summary(m3)$coefficients,"clipboard",sep="\t")

m3=gam(Mean_trophic_level ~ s(scale(log(Surface)))+ Distance
       +  scale(log(Profondeur_maximum))
       + Vegetation_aquatique 
       + Hydroperiode + Anthropique
       , data = data)


summary(m3)

plot(m3)



#########Mean_Omnivory_Index

m4=lm(Mean_Omnivory_Index ~  scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data)

anova(m4, test="Chisq")

  summary(m4)

write.table(summary(m4)$coefficients,"clipboard",sep="\t")

m4=gam(Mean_Omnivory_Index ~ s(scale(log(Surface)))+ Distance
    +  scale(log(Profondeur_maximum))
    + Vegetation_aquatique 
    + Hydroperiode + Anthropique
    , data = data)

summary(m4)

plot(m4)




