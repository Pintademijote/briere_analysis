library(ggplot2)


sim_TTIB=simulation_TTIB[,c(1,2,4,6,8,10,12,23,25,31)]
colnames(sim_TTIB)=colnames(data_observed[,c(1:7,18,20,26)])
sim_TIB=simulation_TIB[,c(1,2,4,6,8,10,12,23,25,31)]
colnames(sim_TIB)=colnames(data_observed[,c(1:7,18,20,26)])

data_ALL=rbind(data_observed[,c(1:7,18,20,26)],sim_TTIB,sim_TIB)

data_ALL$Simu=rep(NA,420)
data_ALL[1:140,11]="Obs"
data_ALL[141:280,11]="TTIB"
data_ALL[281:420,11]="TIB"


layout(matrix(c(1,3,5,2,4,6,7,7,7),3,3, byrow = T), heights = c(0.45,0.45,0.1))

sur=expression(paste(log(surface),"(",m^2,")"))


p1=ggplot(data_ALL, aes(x=Distance, y=Modularity, color=Simu)) +
  geom_point(aes(shape = Simu), size = 2.3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "gam", se=F)+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_color_manual(values=c("black", "grey40", "grey70"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=26, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=18, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        legend.position="none")+
  labs(x = "")


p1
p2=ggplot(data_ALL, aes(x=log(Surface), y=Modularity, color=Simu)) +
  geom_point(aes(shape = Simu), size = 2.3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "gam", se=F)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_color_manual(values=c("black", "grey40", "grey70"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=26, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=18, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        legend.position="none")+
  labs(x = "")


p3=ggplot(data_ALL, aes(x=Distance, y=n, color=Simu)) +
  geom_point(aes(shape = Simu), size = 2.3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "gam", se=F)+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_color_manual(values=c("black", "grey40", "grey70"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=26, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=18, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        legend.position="none")+
  labs(x = "Distance", y="Species richness")


p4=ggplot(data_ALL, aes(x=log(Surface), y=n, color=Simu)) +
  geom_point(aes(shape = Simu), size = 2.3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "gam", se=F)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_color_manual(values=c("black", "grey40", "grey70"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=26, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=18, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        legend.position="none")+
  labs(x = sur, y="Species richness")

p5=ggplot(data_ALL, aes(x=Distance, y=Mean_trophic_level, color=Simu)) +
  geom_point(aes(shape = Simu), size = 2.3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "gam", se=F)+
  scale_x_continuous(breaks = seq(0, 2500, by = 500), labels= seq(0, 2500, by = 500))+
  scale_color_manual(values=c("black", "grey40", "grey70"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=26, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=18, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        legend.position="none")+
  labs(x = "", y="Mean trophic level")


p6=ggplot(data_ALL, aes(x=log(Surface), y=Mean_trophic_level, color=Simu)) +
  geom_point(aes(shape = Simu), size = 2.3, alpha=0.8) + geom_smooth(aes(linetype = Simu), method = "gam", se=F)+
  scale_x_continuous(breaks = seq(0, 9, by = 1), labels= seq(0, 9, by = 1))+
  scale_color_manual(values=c("black", "grey40", "grey70"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(color="black", size=26, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=18, colour="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "Black"),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        legend.position="none")+
  labs(x = "", y="Mean trophic level")


library(ggpubr)


ggarrange(p1,p3,p5,p2,p4,p6, ncol=3,nrow=2,labels=c("A","B","C","D","E","F"),
          legend="bottom", common.legend = T)
