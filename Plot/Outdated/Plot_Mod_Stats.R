  library("ggplot2")
  require("nlme")
  require("mgcv")
  
  setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
  source("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/matrix_function.R")
  
  data_observed=read.delim("data_observed.txt")
  simulation_TTIB=read.delim("simulation_TTIB.txt")
  simulation_TIB=read.delim("simulation_TIB.txt")
  div_marais=read.delim("div_marais.txt")
  
  simulation_TIB=simulation_TIB[simulation_TIB$Site %in% data_observed$Site,]
  simulation_TTIB=simulation_TTIB[simulation_TTIB$Site %in% data_observed$Site,]
  
  layout(matrix(c(1,3,5,2,4,6,7,7,7),3,3, byrow = T), heights = c(0.45,0.45,0.1))
  
  sur=expression(paste(log(surface),"(",m^2,")"))
  
  #Mean trophic level
  par(mar=c(3,4,1,0.1))
  
  plot(data_observed$Mean_trophic_level ~ log(data_observed$Surface), ylim=c(1,4), xlim=c(2,9), xlab="",ylab="")
  abline(m2_Obs_Surf$coefficients[1],m2_Obs_Surf$coefficients[2])
  par(new=T)
  plot(simulation_TTIB$Mean_trophic_level_mean~log(simulation_TTIB$Surface), ylim=c(1,4), xlim=c(2,9),col="blue", xlab="",ylab="")
  abline(m2_TTIB_Surf$coefficients[1],m2_TTIB_Surf$coefficients[2], col=4)
  par(new=T)
  plot(simulation_TIB$Mean_trophic_level_mean~log(simulation_TIB$Surface), ylim=c(1,4), xlim=c(2,9),col="red", xlab="",ylab="")
  abline(m2_TIB_Surf$coefficients[1],m2_TIB_Surf$coefficients[2], col=2)
  par(new=F)
  
  mtext("Mean trophic level",side = 2, line =2.5)
  
  plot(data_observed$Mean_trophic_level ~ data_observed$Distance, ylim=c(1,4), xlim=c(0,3000), xlab="",ylab="")
  abline(m2_Obs_Dist$coefficients[1],m2_Obs_Dist$coefficients[4])
  par(new=T)
  plot(simulation_TTIB$Mean_trophic_level_mean ~ simulation_TTIB$Distance, ylim=c(1,4), xlim=c(0,3000), col="blue", xlab="",ylab="")
  abline(m2_TTIB_Dist$coefficients[1],m2_TTIB_Dist$coefficients[3],col=4)
  par(new=T)
  plot(simulation_TIB$Mean_trophic_level_mean ~ simulation_TIB$Distance, ylim=c(1,4), xlim=c(0,3000), col="red", xlab="",ylab="")
  abline(m2_TIB_Dist$coefficients[1],m2_TIB_Dist$coefficients[3], col=2)
  
  mtext("Mean trophic level",side = 2, line =2.5)
  #Modularity
  
  
  plot(data_observed$Modularity ~ log(data_observed$Surface), ylim=c(0,0.4), xlim=c(2,9), xlab="",ylab="")
  abline(m3_Obs_Surf$coefficients[1],m3_Obs_Surf$coefficients[2])
  par(new=T)
  plot(simulation_TTIB$Modularity_mean~log(simulation_TTIB$Surface),ylim=c(0,0.4), xlim=c(2,9),col="blue", xlab="",ylab="")
  abline(m3_TTIB_Surf$coefficients[1],m3_TTIB_Surf$coefficients[2], col=4)
  par(new=T)
  plot(simulation_TIB$Modularity_mean~log(simulation_TIB$Surface), ylim=c(0,0.4), xlim=c(2,9),col="red", xlab="",ylab="")
  abline(m3_TIB_Surf$coefficients[1],m3_TIB_Surf$coefficients[2], col=2)
  par(new=F)
  
  
  mtext("Modularity",side = 2, line =2.5)
  mtext(sur,side = 1, line =2.5)
  
  plot(data_observed$Modularity ~ data_observed$Distance, ylim=c(0,0.4), xlim=c(0,3000), xlab="",ylab="")
  abline(m3_Obs_Dist$coefficients[1],m3_Obs_Dist$coefficients[4])
  par(new=T)
  plot(simulation_TTIB$Modularity_mean ~ simulation_TTIB$Distance, ylim=c(0,0.4), xlim=c(0,3000), col="blue", xlab="",ylab="")
  abline(m3_TTIB_Dist$coefficients[1],m3_TTIB_Dist$coefficients[3], col=4)
  par(new=T)
  plot(simulation_TIB$Modularity_mean ~ simulation_TIB$Distance, ylim=c(0,0.4), xlim=c(0,3000), col="red", xlab="",ylab="")
  abline(m3_TIB_Dist$coefficients[1],m3_TIB_Dist$coefficients[3], col=2)
  
  mtext("Modularity",side = 2, line =2.5)
  mtext("Distance",side = 1, line =2.5)
  
  #N
  
  plot(data_observed$n ~ log(data_observed$Surface), ylim=c(0,60), xlim=c(2,9), xlab="",ylab="")
  abline(m0_Obs_Surf$coefficients[1],m0_Obs_Surf$coefficients[2])
  par(new=T)
  plot(simulation_TTIB$n_mean~log(simulation_TTIB$Surface), ylim=c(0,60), xlim=c(2,9),col="blue", xlab="",ylab="")
  abline(m1_TTIB_Surf$coefficients[1],m1_TTIB_Surf$coefficients[2],col=4)
  par(new=T)
  plot(simulation_TIB$n_mean~log(simulation_TIB$Surface), ylim=c(0,60), xlim=c(2,9),col="red", xlab="",ylab="")
  abline(m1_TIB_Surf$coefficients[1],m1_TIB_Surf$coefficients[2],col=2)
  par(new=F)
  
  mtext("Species richness",side = 2, line =2.5)
  
  plot(data_observed$n ~ data_observed$Distance, ylim=c(0,65), xlim=c(0,3000), xlab="",ylab="")
  abline(m0_Obs_Dist$coefficients[1],m0_Obs_Dist$coefficients[4])
  par(new=T)
  plot(simulation_TTIB$n_mean~simulation_TTIB$Distance, ylim=c(0,65), xlim=c(0,3000),col="blue", xlab="",ylab="")
  abline(m1_TTIB_Dist$coefficients[1],m1_TTIB_Dist$coefficients[3],col=4)
  par(new=T)
  plot(simulation_TIB$n_mean~simulation_TIB$Distance, ylim=c(0,65), xlim=c(0,3000),col="red", xlab="",ylab="")
  abline(m1_TIB_Dist$coefficients[1],m1_TIB_Dist$coefficients[3],col=2)
  par(new=F)
  
  mtext("Species richness",side = 2, line =2.5)
  
  ####################
  
  par(mar=c(0,1,1,0))
  plot(1,type = "n",axes = F,xlab = "",ylab = "")
  
  legend(x="top", inset = 0, legend = c("Observed", "TIB", "TTIB"), col=c("Black","Red","Blue")
         ,lwd=5,cex=1,horiz=T)
