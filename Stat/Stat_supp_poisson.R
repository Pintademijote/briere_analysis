setwd("C:/Users/pglem/Documents/Master/Stage M1/Rendu/Data/Data")
#setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/Data")

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


#data=data[!is.na(data$Hydroperiode),]
write.table(datacount,"data_observed.txt", sep="\t")
data=datacount
#####

colnames(Matrix_Reg_Manual)

pois=c("Anguille","Breme spp.","Brochet","Carassin","Carpe","Epinoche","Epinochette",
  "Gambusie","Gardon","Perche-soleil","Perche franche","Poisson-chat")

mean(colSums(Matrix_Reg_Manual[,pois]))
sd(colSums(Matrix_Reg_Manual[,pois]))

mean(colSums(Matrix_Reg_Manual[,!(colnames(Matrix_Reg_Manual) %in% pois)]))
sd(colSums(Matrix_Reg_Manual[,!(colnames(Matrix_Reg_Manual) %in% pois)]))

wilcox.test(colSums(Matrix_Reg_Manual[,pois]),colSums(Matrix_Reg_Manual[,!(colnames(Matrix_Reg_Manual) %in% pois)]))

length(colSums(Matrix_Reg_Manual[,pois]))
length(colSums(Matrix_Reg_Manual[,!(colnames(Matrix_Reg_Manual) %in% pois)]))

nlinks=matrix(ncol=2,nrow=length(colSums(Matrix_Reg_Manual[,pois]))+length(colSums(Matrix_Reg_Manual[,!(colnames(Matrix_Reg_Manual) %in% pois)])))
colnames(nlinks)=c("Liens","Cat")
nlinks=as.data.frame(nlinks)
nlinks[1:length(colSums(Matrix_Reg_Manual[,pois])),1]=colSums(Matrix_Reg_Manual[,pois])
nlinks[1:length(colSums(Matrix_Reg_Manual[,pois])),2]="Poissons"
nlinks[(length(colSums(Matrix_Reg_Manual[,pois]))+1):length(nlinks[,1]),1]=colSums(Matrix_Reg_Manual[,!(colnames(Matrix_Reg_Manual) %in% pois)])
nlinks[(length(colSums(Matrix_Reg_Manual[,pois]))+1):length(nlinks[,1]),2]='Autres'


BOXPLOT=ggplot(nlinks, aes(x=Cat,y=Liens))+
  geom_boxplot(aes(fill=Cat), fatten = 5)+
  xlab("") + ylab("Nombre de \n liens de prédations") +
  expand_limits(y=70)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=30, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=26),
        axis.text=element_text(size=20, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+ 
  scale_fill_manual(values = c("gray80", "gray95"))+ 
  guides(fill=FALSE)

df2 <- data.frame(a = c(1, 1,2, 2), b = c(60, 62, 62, 60))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=63), label = "***",size=8)


mean(colSums(Matrix_Reg_Manual[pois,pois]))
mean(colSums(Matrix_Reg_Manual[!(colnames(Matrix_Reg_Manual) %in% pois),pois]))



data$Poissons=as.factor(data$Poissons)
levels(data$Poissons)=c("Absent","Présent")

wilcox.test(data$n~data$Poissons)
mean(data[data$Poissons=="Présent","n"])
mean(data[data$Poissons=="Absent","n"])


wilcox.test(log(data$Surface)~data$Poissons)
mean(log(data[data$Poissons=="Présent","Surface"]))
mean(log(data[data$Poissons=="Absent","Surface"]))
length(data[data$Poissons=="Présent","Surface"])
length(data[data$Poissons=="Absent","Surface"])

wilcox.test(log(data$Surface)~data$Anthropique)
mean(log(data[data$Anthropique=="Agricole","Surface"]))
mean(log(data[data$Anthropique=="Recreationnel","Surface"]))


length(data[data$Anthropique=="Recreationnel","Poissons"])
length(data[data$Anthropique=="Recreationnel" & data$Poissons=="Presence","Poissons"])
27/59*100

BOXPLOT=ggplot(data, aes(x=Poissons,y=log(Surface)))+
  geom_boxplot(aes(fill=Poissons), fatten = 5)+
  xlab("") + ylab("Surface") +
  expand_limits(y=9)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=30, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=26),
        axis.text=element_text(size=20, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+ 
  scale_fill_manual(values = c("gray80", "gray95"))+ 
  guides(fill=FALSE)

df2 <- data.frame(a = c(1, 1,2, 2), b = c(9.5, 9.7, 9.7, 9.5))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=9.8), label = "***",size=8)





BOXPLOT=ggplot(data, aes(x=Poissons,y=log(Surface)))+
  geom_boxplot(aes(fill=Poissons), fatten = 5)+
  xlab("Poissons") + ylab("log(Surface)") +
  expand_limits(y=9)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=30, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=26),
        axis.text=element_text(size=20, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+ 
  scale_fill_manual(values = c("gray80", "gray95"))+ 
  guides(fill=FALSE)

df2 <- data.frame(a = c(1, 1,2, 2), b = c(9.5, 10, 10, 9.5))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=10.1), label = "***",size=8)







BOXPLOT=ggplot(data, aes(x=Anthropique,y=log(Surface)))+
  geom_boxplot(aes(fill=Anthropique), fatten = 5)+
  xlab("") + ylab("Richesse spécifique") +
  expand_limits(y=9)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=30, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=26),
        axis.text=element_text(size=20, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+ 
  scale_fill_manual(values = c("gray80", "gray95"))+ 
  guides(fill=FALSE)

df2 <- data.frame(a = c(1, 1,2, 2), b = c(9.5, 10, 10, 9.5))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=10.1), label = "***",size=8)



BOXPLOT=ggplot(data, aes(x=Vegetation_aquatique,y=n))+
  geom_boxplot(aes(fill=Vegetation_aquatique), fatten = 5)+
  xlab("") + ylab("Richesse spécifique") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=30, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=26),
        axis.text=element_text(size=20, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+ 
  scale_fill_manual(values = c("gray40","gray80", "gray95"))+ 
  guides(fill=FALSE)

df2 <- data.frame(a = c(1, 1,2, 2), b = c(9.5, 10, 10, 9.5))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=10.1), label = "***",size=8)

BOXPLOT=ggplot(data, aes(x=Poissons,y=n))+
  geom_boxplot(aes(fill=Poissons), fatten = 5)+
  xlab("Poissons") + ylab("Richesse spécifique") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=30, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=26),
        axis.text=element_text(size=20, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+ 
  scale_fill_manual(values = c("gray40","gray80", "gray95"))+ 
  guides(fill=FALSE)

df2 <- data.frame(a = c(1, 1,2, 2), b = c(35, 36, 36, 35))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=37), label = "***",size=8)

wilcox.test(data$n~data$Poissons)
