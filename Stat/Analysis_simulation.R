require("enaR")
require("igraph")
require("NetIndices")


setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/PG_Sept2018/Data/Data")
# Matrice régionale d'adjascence
Matrix_Reg_Manual=read.delim2("Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
# var envi
envi_data=read.table("envi_data.txt",header=T)
# Taxa
tab_releve=read.delim("tab_releve.txt")
# distances
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
# Sorties des simulations
simulation_TTIB=read.delim("simulation_TTIB_Final.txt")
#simulation_TTIB=simulation_TTIB[!is.na(simulation_TTIB$Hydroperiode),]
simulation_TIB=read.delim("simulation_TIB_Final.txt")
#simulation_TIB=simulation_TIB[!is.na(simulation_TIB$Hydroperiode),]
# Output from function creata_matrix_marre (save times):
data=read.delim("data_observed.txt")
####
source("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/PG_Sept2018/Data/Script_Function/matrix_function.R")
####
data_format(tab_releve,envi_data,dist)


#####



TIB_M=lm(modularity_t ~ log(Surface)+ Distance, data = simulation_TIB)
summary(TIB_M)
TIB_N=lm(S ~ log(Surface)+ Distance, data = simulation_TIB)
summary(TIB_N)
TIB_C=lm(C ~ log(Surface)+ Distance, data = simulation_TIB)
summary(TIB_C)
TIB_TL=lm(Mean_Trophic_Level ~ log(Surface)+ Distance, data = simulation_TIB)
summary(TIB_TL)
TIB_OI=lm(Mean_Omnivory ~ log(Surface)+ Distance, data = simulation_TIB)
summary(TIB_OI)


TTIB_M=lm(modularity_t ~ log(Surface)+ Distance, data = simulation_TTIB)
summary(TTIB_M)
TTIB_N=lm(S ~ log(Surface)+ Distance, data = simulation_TTIB)
summary(TTIB_N)
TTIB_C=lm(C ~ log(Surface)+ Distance, data = simulation_TTIB)
summary(TTIB_C)
TTIB_TL=lm(Mean_Trophic_Level ~ log(Surface)+ Distance, data = simulation_TTIB)
summary(TTIB_TL)
TTIB_OI=lm(Mean_Omnivory ~ log(Surface)+ Distance, data = simulation_TTIB)
summary(TTIB_OI)
