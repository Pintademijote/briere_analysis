require("enaR")
require("igraph")
require("NetIndices")
require("nlme")
require("mgcv")
library("ggplot2")
library(cheddar)

# Import trophic link database
tab_trophic=read.delim("./Data/tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
# Import taxon occurrence data
tab_releve=read.delim("./Data/tab_releve.txt")
# Import manual matrix
Matrix_Reg_Manual=read.delim2("./Data/Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)

sites=unique(tab_releve$Id_Sites)
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

unique(tab_releve$Taxon)[which(!(unique(tab_releve$Taxon) %in% colnames(Matrix_Reg_Manual)))]

tab_releve=tab_releve[tab_releve$Taxon!="Couleuvre a collier",]

# Import environmental variables
envi_data=read.table("./Data/envi_data.txt",header=T)
# Import distance from marsh
dist=read.delim2("./Data/Distance.txt")
colnames(dist)[2]="Distance"

#data=read.delim("data.txt")

# Import function to create matrix
#source("S:/2018/pglemasle/Data/Script/matrix_function.R")
source("./Script_Function/matrix_function.R")
##########Marais/Mare

releve_marais=tab_releve[tab_releve$Milieux=="Marais",]
releve_marais=droplevels(releve_marais)
marais=unique(releve_marais$Id_Sites)

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
  
  div[i]=length(unique(releve_marais[releve_marais$Id_Sites==marais[i],7]))
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

hist(meanresults[1:155,1], breaks=100)



BOXPLOT=ggplot(meanresults, aes(x=Lieu,y=Richesse_spécifique))+
  geom_boxplot(aes(fill=Lieu), fatten = 5)+
  xlab("") + ylab("Richesse_spécifique") +
  expand_limits(y=45)+
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

df2 <- data.frame(a = c(1, 1,2, 2), b = c(40, 41, 41, 40))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=42), label = "***",size=8)
