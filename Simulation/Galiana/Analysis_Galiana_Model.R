setwd("C:/Users/pglem/Desktop/Code_Galiana_EtAl.2018")

TIB=read.csv("output_tib_test.csv")
TTIB=read.csv("output_ttib_test.csv")
TTIB_corrected=read.csv("output_ttib_test_corrected.csv")
  

library(ggplot2)
library(ggpubr)

results=rbind(TIB,TTIB,TTIB_corrected)

results$Simu='TTIB_cor'

results[1:128117,22]="TIB"
results[128118:243000,22]="TTIB"

results=results[results$colonization_rate==0.2,]
View(results[results$Mean_Trophic_Level>5,])

p1=ggplot(data=results, aes(x=ratio,y=S, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p2=ggplot(data=results, aes(x=ratio,y=modularity_t, color=Simu))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p3=ggplot(data=results, aes(x=ratio,y=Mean_Trophic_Level, color=Simu))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p4=ggplot(data=results, aes(x=ratio,y=Mean_Omnivory, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)



p5=ggplot(data=results, aes(x=ratio,y=top, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p6=ggplot(data=results, aes(x=ratio,y=intermediate, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p7=ggplot(data=results, aes(x=ratio,y=basal, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p8=ggplot(data=results, aes(x=ratio,y=indegree, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p9=ggplot(data=results, aes(x=ratio,y=outdegree, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p10=ggplot(data=results, aes(x=ratio,y=C, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)



ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol=2,nrow=5,common.legend = T)

#####

#results=results[results$colonization_rate==0.1,]

p1=ggplot(data=results, aes(x=extinction_rate,y=S, color=Simu))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p2=ggplot(data=results, aes(x=extinction_rate,y=modularity_t, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p3=ggplot(data=results, aes(x=extinction_rate,y=Mean_Trophic_Level, color=Simu))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p4=ggplot(data=results, aes(x=extinction_rate,y=Mean_Omnivory, color=Simu))+
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)



p5=ggplot(data=results, aes(x=extinction_rate,y=top, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p6=ggplot(data=results, aes(x=extinction_rate,y=intermediate, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p7=ggplot(data=results, aes(x=extinction_rate,y=basal, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p8=ggplot(data=results, aes(x=extinction_rate,y=indegree, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p9=ggplot(data=results, aes(x=extinction_rate,y=outdegree, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p10=ggplot(data=results, aes(x=extinction_rate,y=C, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol=2,nrow=5,common.legend = T)

#####

#results=results[results$extinction_rate==0.1,]

p1=ggplot(data=results, aes(x=colonization_rate,y=S, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p2=ggplot(data=results, aes(x=colonization_rate,y=modularity_t, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p3=ggplot(data=results, aes(x=colonization_rate,y=Mean_Trophic_Level, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p4=ggplot(data=results, aes(x=colonization_rate,y=Mean_Omnivory, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)




p5=ggplot(data=results, aes(x=colonization_rate,y=top, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p6=ggplot(data=results, aes(x=colonization_rate,y=intermediate, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p7=ggplot(data=results, aes(x=colonization_rate,y=basal, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p8=ggplot(data=results, aes(x=colonization_rate,y=indegree, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p9=ggplot(data=results, aes(x=colonization_rate,y=outdegree, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)

p10=ggplot(data=results, aes(x=colonization_rate,y=C, color=Simu))+
  geom_smooth(method = "gam", formula = y ~ s(x), size = 1,
              se=T, level = 0.95)



ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol=2,nrow=5,common.legend = T)

