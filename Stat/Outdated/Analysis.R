  setwd("S:/2018/pglemasle/Data/Script/")
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
  # Import environmental variables
  envi_data=read.table("envi_data.txt",header=T)
  # Import distance from marsh
  dist=read.delim2("Distance.txt")
  colnames(dist)[2]="Distance"
  # Import manual matrix
  Matrix_Reg_Manual=read.delim2("Matrix_Manual_Chekck.txt", row.names = 1)
  colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
  Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
  data=read.delim("data.txt")
  
  # Import function to create matrix
  #source("S:/2018/pglemasle/Data/Script/matrix_function.R")
  source("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/matrix_function.R")
  ##########Marais/Mare
  
  releve_marais=tab_releve[tab_releve$Milieux=="Marais",]
  releve_marais=droplevels(releve_marais)
  marais=levels(releve_marais$Id_Sites)
  
  ##########
  
  data_format(tab_releve,envi_data,dist)
  
  create_matrix(tab_trophic)
  
  #create_matrix_marre(tab_releve,Matrix_tab_trophic,envi_data,Count = F, year=F,name_tab = "data")
  
  create_matrix_marre(tab_releve,Matrix_Reg_Manual,envi_data,Count = F, year=F,name_tab = "data")

  
  data=data[!is.na(data$Hydroperiode),]
  write.table(data,"data_observed.txt", sep="\t")

#########

div=c(1:length(marais))
for (i in 1:length(marais)) {
  
  div[i]=length(levels(droplevels(releve_marais[releve_marais$Id_Sites==marais[i],7])))
}

mean(div)

mean(data$n)

meanresults=matrix(NA,nrow=length(div)+length(data$n),ncol=2)
meanresults=as.data.frame(meanresults)
meanresults[1:140,1]=data$n
meanresults[1:140,2]="Ponds"
meanresults[141:158,1]=div
meanresults[141:158,2]="Marsh"
colnames(meanresults)=c("Richesse_spécifique","Lieu")

hist(meanresults[1:155,1], breaks=100)

wilcox.test(meanresults[1:140,1],meanresults[141:158,1])

########################## RELATION AMONG DESCRIPTORS
# Environmental
pairs(data[,c(13:15,17:20,25)])
# Topological
pairs(data[,c(2:7)])

cor_envi=cor(data[,c(13:15,17:20,25)])
cor_envi>=0.4 | cor_envi<=-0.4

cor_topo=cor(data[,c(2:7)])
cor_topo>=0.4 | cor_topo<=-0.4

########################## CONNECTANCE ANALYSIS
plot(data$C ~ data$Surface)
plot(data$C ~ data$Distance)

m0 <- lm(C ~ scale(log(Surface)) +  scale(log(Profondeur_maximum))  +
           Vegetation_aquatique+ Distance + Hydroperiode + Anthropique, data = data)
predict_C=predict(m0)
summary(m0)

write.table(m0$coefficients,"clipboard",sep="\t")


########################## Specific richness ANALYSIS
plot(data$n ~ data$Surface)
plot(data$n ~ data$Distance)
hist(data$n)

m0_Obs <- glm(n ~ scale(log(Surface)) + scale(I(log(Surface)^2)) + Distance+
           Vegetation_aquatique + Hydroperiode + Anthropique, data = data, family = quasipoisson(link = "log"))

predict_n=predict(m)

m=summary(m)
write.table(m$coefficients,"clipboard",sep="\t")

m <- gam(n ~ s(scale(log(Surface)))+ Distance+
           Vegetation_aquatique + Hydroperiode + Anthropique, data = data, family = quasipoisson(link = "log"))
plot(m)
m=summary(m)
write.table(m$p.table,"clipboard",sep="\t")
write.table(m$s.table,"clipboard",sep="\t")

########################## MEAN TROPHIC LEVEL ANALYSIS
hist(data$Mean_trophic_level, breaks = 50)
plot(log(data$Mean_trophic_level) ~ log(data$Surface))
m2_Obs <- lm(log(Mean_trophic_level) ~ scale(log(Surface))+ Distance +Vegetation_aquatique + Hydroperiode + Anthropique, data = data)
summary(m2)

predict_TL=predict(m2)

m2=summary(m2)
write.table(m2$coefficients,"clipboard",sep="\t")

m2 <- gam(log(Mean_trophic_level) ~ s(scale(log(Surface))) + Distance+Vegetation_aquatique + Hydroperiode + Anthropique, data = data)
plot(m2)
summary(m2)
########################## Modularity ANALYSIS
hist(data$Modularity, breaks = 50)
plot(data$Modularity ~ log(data$Surface))
m3_Obs <- lm(Modularity ~ scale(log(Surface)) + scale(I(log(Surface)^2)) + + Distance+Vegetation_aquatique+Hydroperiode + Anthropique, data = data)
summary(m3_Obs)

predict_modularity=predict(m3)

m3=summary(m3)
write.table(m3$coefficients,"clipboard",sep="\t")

m3 <- gam(Modularity ~ s(scale(log(Surface))) + Distance+Vegetation_aquatique + Hydroperiode + Anthropique, data = data)
summary(m3)
plot(m3)
m3=summary(m3)
write.table(m3$p.table,"clipboard",sep="\t")
write.table(m3$s.table,"clipboard",sep="\t")


#########################





