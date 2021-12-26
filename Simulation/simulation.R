require("enaR")
require("igraph")
require("NetIndices")


setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
Matrix_Reg_Manual=read.delim2("Matrix_Manual_Chekck.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
envi_data=read.table("envi_data.txt",header=T)
tab_releve=read.delim("tab_releve.txt")
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
####
source("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/matrix_function.R")
####
data_format(tab_releve,envi_data,dist)

####
normalized  <- function(x){(x-min(x))/(max(x)-min(x))}

simu_TTIB=function(x,y,t=200,name_tab="data"){
  
  start.time = Sys.time()
  descriptors=matrix(NA,nrow = length(y$Mares), ncol = 8)
  colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level","Mean_Omnivory_Index", "Richesse_spe")
  name_site=y$Mares
  DisP=normalized(y$Distance)
  SurfP=normalized(log(y$Surface))
  P=length(x[,1])
  Base=c("Detritus","Macrophyte","Phytoplankton","Zooplankton")
  pb <- txtProgressBar(min = 0, max = length(DisP), style = 3)
  
  for (mares in 1:length(DisP)) {
    
    Sys.sleep(0.1)		
    setTxtProgressBar(pb, mares)
    temp=matrix(0,length(x[,1]),t)
    row.names(temp)=row.names(x)
    
    temp[Base,1]=1
    
    for (i in 2:t) {
      Colonization=0
      Extinction=0
      C=(DisP[mares]*(P-sum(temp[,i-1])))/100
      E=(SurfP[mares]*sum(temp[,i-1]))/100
      
      if(i>1){
        temp[,i]=temp[,i-1]
        present=names(temp[temp[,i-1]==1,i-1])
        possible=colnames(x[present,colSums(x[present,])!=0])
        possible=possible[!(possible %in% present)]
        Colonization=sample(c(0,1),1,prob = c(1-C,C))
        Extinction=sample(c(0,1),1,prob = c(1-E,E))
        
        if(Colonization==1){
          add=sample(possible,1)
          temp[add,i]=1
        }
        if(Extinction==1 & length(present[!(present %in% Base)])!=0){
          rmv=sample(present[!(present %in% Base)],1)
          temp[rmv,i]=0
          
        }
        
      }
      
    }
    present=names(temp[temp[,i]==1,i])
    
    temp_matrix=x[present,present]
    
    g=graph_from_adjacency_matrix(temp_matrix, mode="undirected")
    descriptors[mares,1]=name_site[mares]
    descriptors[mares,2:4]=enaStructure(temp_matrix)$ns[,c(1,2,3)]
    descriptors[mares,5]=modularity(g,membership(cluster_walktrap(g)))
    descriptors[mares,6]=mean(TrophInd(temp_matrix)[,1])
    descriptors[mares,7]=mean(TrophInd(temp_matrix)[,2])
    descriptors[mares,8]=length(temp_matrix[,1])
    
  }
  descriptors[,1]=as.factor(descriptors[,1])
  assign(paste(name_tab),
         merge(descriptors,y,by.x = c("Site"),by.y = c("Mares")), env = .GlobalEnv) 
close(pb)
close(pb) # Close the progress bar
end.time <- Sys.time() #Establish the time when ending the function
time.taken <- end.time - start.time # Calculate the time the function worked
time.taken 
}

simu_TTIB(Matrix_Reg_Manual,envi_data,t=1000)

require("nlme")
require("mgcv")

########################## RELATION AMONG DESCRIPTORS
# Environmental
pairs(data[,c(18:20,26)])
# Topological
pairs(data[,c(2:7)])

cor_envi=cor(data[,c(18:20,26)])
cor_envi>=0.4 | cor_envi<=-0.4

cor_topo=cor(data[,c(2:7)])
cor_topo>=0.4 | cor_topo<=-0.4

########################## CONNECTANCE ANALYSIS
plot(data$C ~ data$Surface)
plot(data$C ~ data$Distance)

m0 <- lm(C ~ scale(log(Surface)) +  scale(log(Profondeur_maximum))  +
           Vegetation_aquatique+ Distance + Hydroperiode + Anthropique, data = data)
summary(m0)



########################## Specific richness ANALYSIS
plot(data$n ~ data$Surface)
plot(data$n ~ data$Distance)
hist(data$n, breaks =50)

m <- glm(n ~ scale(log(Surface)) + scale(I(log(Surface)^2)) + Distance+
           Vegetation_aquatique + Hydroperiode + Anthropique, data = data, family = quasipoisson(link = "log"))

summary(m)

m <- gam(n ~ s(scale(log(Surface)))+ Distance+
           Vegetation_aquatique + Hydroperiode + Anthropique, data = data, family = quasipoisson(link = "log"))
plot(m)
summary(m)

########################## MEAN TROPHIC LEVEL ANALYSIS
hist(data$Mean_trophic_level, breaks = 50)
plot(log(data$Mean_trophic_level) ~ log(data$Surface))
m2 <- lm(log(Mean_trophic_level) ~ scale(log(Surface)) + scale(I(log(Surface)^2)) + Distance +Vegetation_aquatique + Hydroperiode + Anthropique, data = data)
summary(m2)

m2 <- gam(log(Mean_trophic_level) ~ s(scale(log(Surface))) + Distance+Vegetation_aquatique + Hydroperiode + Anthropique, data = data)
plot(m2)
summary(m2)
########################## Modularity ANALYSIS
hist(data$Modularity, breaks = 50)
plot(data$Modularity ~ log(data$Surface))
m3 <- lm(Modularity ~ scale(log(Surface)) + scale(I(log(Surface)^2)) + + Distance+Vegetation_aquatique+Hydroperiode + Anthropique, data = data)
summary(m3)


#########################





