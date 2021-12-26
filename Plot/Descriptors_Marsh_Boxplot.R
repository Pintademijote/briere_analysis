########Import needed data
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
# Import manual matrix
Matrix_Reg_Manual=read.delim2("Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)

sites=levels(tab_releve$Id_Sites)
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

levels(tab_releve$Taxon)[which(!(levels(tab_releve$Taxon) %in% colnames(Matrix_Reg_Manual)))]

tab_releve=tab_releve[tab_releve$Taxon!="Couleuvre a collier",]

# Import environmental variables
envi_data=read.table("envi_data.txt",header=T)
# Import distance from marsh
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"

#data=read.delim("data.txt")

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

data=datacount
##########Marais/M





x=releve_marais
y=Matrix_Reg_Manual
year=F
Count=F
descriptor=T
Matrix=F
name_tab="data_marais"

doublonPhyto=x[!duplicated(x$Id_Sites),]
doublonPhyto[,c(5,6,8,9)]=NA				# Remove useless data but keep matrix size		
doublonPhyto[,7]="Phytoplankton"
doublonZoo=x[!duplicated(x$Id_Sites),]
doublonZoo[,c(5,6,8,9)]=NA
doublonZoo[,7]="Zooplankton"
doublonDetri=x[!duplicated(x$Id_Sites),]
doublonDetri[,c(5,6,8,9)]=NA
doublonDetri[,7]="Detritus"
doublonMacrop=x[!duplicated(x$Id_Sites),]
doublonMacrop[,c(5,6,8,9)]=NA
doublonMacrop[,7]="Macrophyte"

x=rbind(x,doublonPhyto,doublonZoo,doublonDetri,doublonMacrop)

  
  start.time = Sys.time() #Establish the time when running the function
  
  if(year==F){
    marres=levels(x$Id_Sites) #Assign in an object with the name of the site
    
    
    pb <- txtProgressBar(min = 0, max = length(marres), style = 3) #Create the progress bar
    
    temp_count=matrix(nrow = length(marres),ncol = length(levels(x$Taxon))+1)
    colnames(temp_count)=c(as.character(as.data.frame(table(x[x[,2],7]))[,1]),"Site")
    
    descriptors=as.data.frame(matrix(nrow=length(levels(x$Id_Sites)),ncol=8)) #Create an empty data frame for the descriptors
    colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level","Mean_Omnivory_Index", "Richesse_spe") #Specified the name of the column for the descriptos object
    name_site=levels(x$Id_Sites) #Assign in an object with the name of the site
    name_matrix=1:length(levels(x$Id_Sites)) #Assign a numeric level for the name of matrix
    
    for (o in 1:length(marres)) {#Loop through each ponds
      Sys.sleep(0.1) # Define refresh rate for progress bar
      setTxtProgressBar(pb, o) #Set the current status of the progress bar
      
      temp_count[o,1:length(levels(x$Taxon))]=as.data.frame(table(x[x[,2]==marres[o],7]))[,2]
      temp_count[o,length(levels(x$Taxon))+1]=name_site[o]
      
      temp_site=x[x[,2]==marres[o],] #Subset the taxon occurrence data for a define pond
      temp_site=droplevels(temp_site) #Forget old levels
      
      
      nm=paste("Matrix_site_",marres[o],sep="") #Define the name of the matrix returned for the current pond
      
      
      
      level=as.character(c(levels(temp_site$Taxon)))
      
      temp=y[level,level]
      
      if(descriptor==T){ # If arg "descriptor" is true, calculate the descriptor
        g=graph_from_adjacency_matrix(temp, mode="undirected")
        descriptors[o,1]=name_site[o]
        name_matrix[o]=paste("Matrix_site_",name_site[o],sep="")
        descriptors[o,2:4]=enaStructure(temp)$ns[,c(1,2,3)]
        descriptors[o,5]=modularity(g,membership(cluster_walktrap(g)))
        descriptors[o,6]=mean(TrophInd(temp)[,1])
        descriptors[o,7]=mean(TrophInd(temp)[,2])
        descriptors[o,8]=length(temp[,1])
      }
      if(Matrix==T){assign(nm,temp, env = .GlobalEnv)} # If arg "Matrix" is true, put in global environnement adjancy matrix of the current pond
    }
    
    descriptors$Site=as.factor(descriptors$Site)
    
    if(descriptor==T & Count==F){assign(paste(name_tab),
                                        descriptors, env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
    col_temp=temp_count[,length(levels(x$Taxon))+1]
    temp_count[temp_count>0]=1
    temp_count[,length(levels(x$Taxon))+1]=col_temp
    if(Count==T){assign(paste(name_tab,"count",sep=""),
                        merge(merge(descriptors,z,by.x = c("Site"),by.y = c("Mares")),temp_count,by.x = c("Site"),by.y = c("Site")), env = .GlobalEnv)}
  }
  
  
    }




marais_mares=matrix(NA,ncol=2,nrow = length(data_marais$Mean_Omnivory_Index)+length(data$Mean_Omnivory_Index))
marais_mares=as.data.frame(marais_mares)
colnames(marais_mares)=c("Omn","Site")

marais_mares[1:length(data_marais$Mean_Omnivory_Index),1]=data_marais$Mean_Omnivory_Index
marais_mares[1:length(data_marais$Mean_Omnivory_Index),2]="Marais"
marais_mares[(1+length(data_marais$Mean_Omnivory_Index)):(length(data$Mean_Omnivory_Index)+length(data_marais$Mean_Omnivory_Index)),1]=data$Mean_Omnivory_Index
marais_mares[(1+length(data_marais$Mean_Omnivory_Index)):(length(data$Mean_Omnivory_Index)+length(data_marais$Mean_Omnivory_Index)),2]="Mares"
marais_mares$Site=as.factor(marais_mares$Site)

wilcox.test(marais_mares$Omn~marais_mares$Site)





BOXPLOT=ggplot(marais_mares, aes(x=Site,y=Omn))+
  geom_boxplot(aes(fill=Site), fatten = 5)+
  xlab("") + ylab("Indice d'omnivorie moyen") +
  expand_limits(y=1)+
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

df2 <- data.frame(a = c(1, 1,2, 2), b = c(1.05, 1.07, 1.07, 1.05))

par(mar = c(6,9,4,2) + 0.1)

BOXPLOT+ geom_line(data = df2, aes(x = a, y = b))+geom_text(aes(x=1.5,y=1.1), label = "*",size=8)

