setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
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

create_matrix_marre(tab_releve,Matrix_Reg_Manual,envi_data,Count = T, year=F,name_tab = "data")


#####

p1=ggplot(data = datacount, aes(x=Distance,y=Percentage_top))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p2=ggplot(data = datacount, aes(x=Distance,y=Percentage_basal))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p3=ggplot(data = datacount, aes(x=Distance,y=Percentage_intermediate))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p4=ggplot(data = datacount, aes(x=log(Surface),y=Percentage_top))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p5=ggplot(data = datacount, aes(x=log(Surface),y=Percentage_basal))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p6=ggplot(data = datacount, aes(x=log(Surface),y=Percentage_intermediate))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p7=ggplot(data = datacount, aes(x=log(Surface),y=Percentage_intermediate))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p8=ggplot(data = datacount, aes(x=log(Surface),y=indegree))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p9=ggplot(data = datacount, aes(x=log(Surface),y=outdegree))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,NULL,ncol=2,nrow=5,common.legend = T)
