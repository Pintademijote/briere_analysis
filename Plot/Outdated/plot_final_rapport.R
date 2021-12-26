library(ggplot2)

setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
Matrix_Reg_Manual=read.delim2("Matrix_Manual_Chekck.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
envi_data=read.table("envi_data.txt",header=T)
tab_releve=read.delim("tab_releve.txt")
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
simulation_TTIB=read.delim("simulation_TTIB.txt")
simulation_TTIB=simulation_TTIB[!is.na(simulation_TTIB$Hydroperiode),]
simulation_TIB=read.delim("simulation_TIB.txt")
simulation_TIB=simulation_TIB[!is.na(simulation_TIB$Hydroperiode),]
data_observed=read.delim("data_observed.txt")
sim_TTIB=simulation_TTIB[,c(1,2,4,6,8,10,12,23,25,31)]
colnames(sim_TTIB)=colnames(data_observed[,c(1:7,18,20,26)])
sim_TIB=simulation_TIB[,c(1,2,4,6,8,10,12,23,25,31)]
colnames(sim_TIB)=colnames(data_observed[,c(1:7,18,20,26)])

data_ALL=rbind(data_observed[,c(1:7,18,20,26)],sim_TTIB,sim_TIB)

data_ALL$Simu=rep(NA,420)
data_ALL[1:140,11]="Observed"
data_ALL[141:280,11]="TTIB"
data_ALL[281:420,11]="TIB"


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

p1=p1+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=Modularity,shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=Modularity,linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
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
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "Distance (m)",y="\n Modularity")


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

p3=p3+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=n, shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=n, linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 60, by = 20), labels= c("    0","   20", "   40","   60"), limits=c(0,65))+
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

p5=p5+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=Mean_trophic_level, shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=Mean_trophic_level, linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(1, 3.8, by = 0.5), labels=c(" 1","  1.5", " 2", " 2.5", " 3", " 3.5"))+
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

p9=p9+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=C, shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=Distance, y=C, linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
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
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="\n Connectance")


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



p4=p4+
  geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=n,shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=n,linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
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

p6=p6+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=Mean_trophic_level,shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=Mean_trophic_level,linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
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

p8=p8+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=Mean_Omnivory_Index,shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=Mean_Omnivory_Index,linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
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

p2=p2+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=Modularity,shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=Modularity,linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
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
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "log(size area)",y="")



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

p10=p10+geom_point(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=C, shape = Simu), size = 3, alpha=0.8) +
  geom_smooth(data=data_ALL[data_ALL$Simu==c("TIB", "TTIB"),],aes(x=log(Surface), y=C, linetype = Simu), method = "lm",formula = y ~ x, se=F,size=2)+
  scale_color_manual(values=c("black", "grey70", "grey40"))+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
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
        axis.ticks.length=unit(.25, "cm"),
        plot.margin=unit(c(0.5,0,0,0), "cm"))+
  labs(x = "", y="")


########*


ggarrange(p3,p4,p5,p6,p7,p8,p9,p10,p1,p2,ncol=2,nrow=5,labels=c("A","B"),font.label=c(size=32), vjust=1
          ,legend=c("bottom", "right"),common.legend = T)
