#setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
setwd("C:/Users/pglem/Documents/Master/Stage M1/Rendu/Data/Data")

require("enaR")
require("igraph")
require("NetIndices")
require("nlme")
require("mgcv")
library("ggplot2")
library(ade4)
library(ade4TkGUI)

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


###############

data_topo=datacount[,c(2:18,26:30,33:36)] #Keep only non correlated variables
row.names(data_topo)=datacount$Site

data_topo$Vegetation_aquatique=as.character(data_topo$Vegetation_aquatique)

data_topo[data_topo[,18]=="peu",18]=1
data_topo[data_topo[,18]=="intermediaire",18]=2
data_topo[data_topo[,18]=="beaucoup",18]=3

data_topo$Vegetation_aquatique=as.numeric(data_topo$Vegetation_aquatique)
data_topo$Hydroperiode=as.numeric(data_topo$Hydroperiode)
data_topo$Anthropique=as.numeric(data_topo$Anthropique)
nahydro=!(is.na(data_topo$Hydroperiode))
data_topo=data_topo[nahydro,]

data_topo=data_topo[,-c(2,3,6,8,9,10,14:17,19,25)]

data_spe=matrix(as.numeric(as.matrix(datacount[,37:94])), nrow = 155, ncol =58 ) #Keep only species count
row.names(data_spe)=datacount$Site
colnames(data_spe)=colnames(datacount[,37:94])
data_spe=data_spe[nahydro,]


test=cor(data_topo)>=0.8 | cor(data_topo)<=-0.8 #Verified correlation

cor(data_spe)>=0.8 | cor(data_spe)<=-0.8




data_spe=data_spe[,!((colSums(data_spe)/sum(data_spe)*100)<1)] #Remove rare species representing less than 1% of total count

##########

AFC=dudi.coa(df = data_spe, scannf = FALSE, nf = 2) #COA on species count

barplot(AFC$eig*100, main="Pourcentage d'inertie expliqué", ylim=c(0,35), col = )




s.label(AFC$li, label = row.names(data_spe) )
s.label(AFC$co)

contribution<-inertia.dudi(AFC, col.inertia=TRUE)
seuil=100/length(as.data.frame(data_spe))

choice=contribution$col.abs
choice>seuil


############

ACP=dudi.pca(df = data_topo, row.w = AFC$lw, scannf = FALSE, nf = 2) #PCA on 

barplot(ACP$eig/sum(ACP$eig)*100, main="Pourcentage d'inertie expliquÃ©", ylim=c(0,80))



s.corcircle(ACP$co)

s.label(ACP$li )



contribution<-inertia.dudi(ACP, col.inertia=TRUE)
seuil=100/length(as.data.frame(data_topo))

choice=contribution$col.abs
choice>seuil


##############


COI=coinertia(dudiX = ACP, dudiY = AFC, scannf = FALSE, nf = 2)

randtest(COI)

plot(randtest(COI))


iner=inertia.dudi(COI,col.inertia=T,row.inertia=T)

iner

summary(COI)

plot(COI, col=c(1:8))

ade4TkGUI()
