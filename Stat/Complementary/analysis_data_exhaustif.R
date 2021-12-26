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
tab_releve=tab_releve[tab_releve$Qualite.de.l.echantillonnage=="Assez exhaustif",]
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

data_observed=datacount


data_ALL=data_observed


data_ALL$Simu="Observed"



#########Connectance

m0 <- lm(C ~ scale(log(Surface)) +  scale(log(Profondeur_maximum)) + Distance + Vegetation_aquatique 
         + Hydroperiode + Anthropique
         , data = data_ALL)

anova(m0)

summary(m0)

write.table(summary(m0)$coefficients,"clipboard",sep="\t")

#########SPecies Richness

m1=glm(n ~ scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
       +  scale(log(Profondeur_maximum))
       + Vegetation_aquatique 
       + Hydroperiode + Anthropique
       , data = data_ALL, family = quasipoisson(link = "log"))

anova(m1, test="Chisq")

summary(m1)

write.table(summary(m1)$coefficients,"clipboard",sep="\t")

#########Modularity
m2=lm(Modularity ~ scale(log(Surface))+scale(I(log(Surface)^2)) + Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data_ALL)

anova(m2, test="Chisq")

summary(m2)

write.table(summary(m2)$coefficients,"clipboard",sep="\t")

m2=gam(n ~ s(scale(log(Surface)))+ Distance
       +  scale(log(Profondeur_maximum))
       + Vegetation_aquatique 
       + Hydroperiode + Anthropique
       , data = data_ALL)


summary(m2)

plot(m2)

#########Mean Trophic Level

m3=lm(Mean_trophic_level ~ scale(log(Surface))+scale(I(log(Surface)^2))+ Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data_ALL)

anova(m3, test="Chisq")

summary(m3)

write.table(summary(m3)$coefficients,"clipboard",sep="\t")

m3=gam(Mean_trophic_level ~ s(scale(log(Surface)))+ Distance
       +  scale(log(Profondeur_maximum))
       + Vegetation_aquatique 
       + Hydroperiode + Anthropique
       , data = data_ALL)


summary(m3)

plot(m3)



#########Mean_Omnivory_Index

m4=lm(Mean_Omnivory_Index ~ scale(log(Surface))+scale(I(log(Surface)^2))+ Distance
      +  scale(log(Profondeur_maximum))
      + Vegetation_aquatique 
      + Hydroperiode + Anthropique
      , data = data_ALL)

anova(m4, test="Chisq")

summary(m4)

write.table(summary(m4)$coefficients,"clipboard",sep="\t")

m4=gam(Mean_Omnivory_Index ~ s(scale(log(Surface)))+ Distance
       +  scale(log(Profondeur_maximum))
       + Vegetation_aquatique 
       + Hydroperiode + Anthropique
       , data = data_ALL)

summary(m4)

plot(m4)








layout(matrix(c(1,3,5,2,4,6,7,7,7),3,3, byrow = T), heights = c(0.45,0.45,0.1))

sur=expression(paste(log(surface),"(",m^2,")"))

p1=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=Distance, y=Modularity, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) +
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.1), labels= c(" 0"," 0.1", "0.2","  0.3", "0.4"), limits=c(0,0.4))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "",y="\n Modularity")




p3=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=Distance, y=n, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) +
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 60, by = 20), labels= c("    0","   20", "   40","   60"), limits=c(0,65))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="\n Species richness")



p5=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=Distance, y=Mean_trophic_level, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) +
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(1, 3.8, by = 0.5), labels=c(" 1","  1.5", " 2", " 2.5", " 3", " 3.5"))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="Mean \n trophic level")





p7=ggplot(data_ALL, aes(x=Distance, y=Mean_Omnivory_Index, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="Mean \n Omnivory Index")

p9=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=Distance, y=C, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) +
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="Mean \n trophic level")




#####Surface

p4=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=log(Surface), y=n, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_y_continuous(breaks = seq(0, 60, by = 20), labels= c("   0","   20", "   40","   60"), limits=c(0,65))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="")







p6=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=log(Surface), y=Mean_trophic_level, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_y_continuous(breaks = seq(1, 3.8, by = 0.5), labels=c(" 1","  1.5", " 2", " 2.5", " 3", " 3.5"))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="")




p8=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=log(Surface), y=Mean_Omnivory_Index, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="")



p2=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=log(Surface), y=Modularity, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.1), labels= c("  0","   0.1", "   0.2","   0.3", "0.4"), limits=c(0,0.4))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "",y="")





p10=ggplot(data_ALL[data_ALL$Simu=="Observed",], aes(x=log(Surface), y=C, color=Simu)) +
  geom_point(aes(shape = Simu), size = 3, alpha=0.8) +
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=26),
        axis.title.y = element_text(color="black", size=28),
        axis.text=element_text(size=21, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=25),
        legend.title=element_blank(),
        legend.position="none",
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="")



########*


ggarrange(p3,p4,p5,p6,p7,p8,p9,p10,p1,p2,ncol=2,nrow=5,labels=c("A","B"),font.label=c(size=32), vjust=1
          ,legend=c("bottom", "right"),common.legend = T)
