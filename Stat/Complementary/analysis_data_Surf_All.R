setwd("C:/Users/pglem/Documents/Master/Stage M1/Rendu/Data/Data")
source('C:/Users/pglem/Documents/Master/Stage M1/Rendu/Data/Script_Function/utils_NAR.R') #This code is used to compute all the network metrics
#setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/Data")

require("enaR")
require("igraph")
require("NetIndices")
require("nlme")
require("mgcv")
library("ggplot2")
source('utils_NAR.R') #This code is used to compute all the network metrics

# Import trophic link database
tab_trophic=read.delim("tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
# Import taxon occurrence data
tab_releve=read.delim("tab_releve.txt")
which(tab_releve$Id_Sites=="")
tab_releve=tab_releve[-which(tab_releve$Id_Sites==""),]
tab_releve=droplevels(tab_releve)

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



tab_releve=tab_releve[tab_releve$Taxon!="Couleuvre a collier",]

  # Import environmental variables
envi_data=read.table("envi_data.txt",header=T)
# Import distance from marsh
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
# Import pond surface
Surf=read.delim2("Surface.txt")
# Import manual matrix
Matrix_Reg_Manual=read.delim2("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
#data=read.delim("data.txt")

# Import function to create matrix
#source("S:/2018/pglemasle/Data/Script/matrix_function.R")

##########Marais/Mare

Surf[,1]=gsub("M_", "", Surf[,1], fixed = TRUE)
Surf[,1]=sub("^0+", "", Surf[,1])
Surf[,1]=factor(Surf[,1])

  
dist[,1]=gsub("Mare_", "", dist[,1], fixed = TRUE)
dist[,1]=sub("^0+", "", dist[,1])
dist[,1]=as.factor(dist[,1])
dist=droplevels(dist)

## Keep only ponds
#tab_releve=tab_releve[tab_releve[,1]=="Mare",]

## Reset factor levels
tab_releve=droplevels(tab_releve)
## Reset format
tab_releve[,2]=as.character(tab_releve[,2], trim = TRUE)
tab_releve[,2]=as.factor(tab_releve[,2])

test=as.factor(tab_releve[!duplicated(x[,2]),2])	# Keep only one line per site
test=as.factor(tab_releve[,2])


#### ADD non-observed taxa at all sites
doublonPhyto=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonPhyto[,c(5,6,8,9)]=NA				# Remove useless data but keep matrix size		
doublonPhyto[,7]="Phytoplankton"
doublonZoo=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonZoo[,c(5,6,8,9)]=NA
doublonZoo[,7]="Zooplankton"
doublonDetri=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonDetri[,c(5,6,8,9)]=NA
doublonDetri[,7]="Detritus"
doublonMacrop=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonMacrop[,c(5,6,8,9)]=NA
doublonMacrop[,7]="Macrophyte"
## Add these taxa to the dataset
tab_releve=rbind(tab_releve,doublonPhyto,doublonZoo,doublonDetri,doublonMacrop)

rm(doublonZoo,doublonPhyto,doublonDetri,doublonMacrop)

tab_releve$Id_Sites=as.factor(tab_releve$Id_Sites)
length(levels(tab_releve$Id_Sites))
length(levels(droplevels(tab_releve[tab_releve$Id_Sites %in% levels(Surf$Id_Mare),2])))


x=tab_releve
y=Matrix_Reg_Manual
z=Surf
Matrix=F
descriptor=T
Count=T
year=F

create_matrix_marre=function(x,y,z,Matrix=F,descriptor=T,Count=F, year=F, name_tab = "data"){
  
  start.time = Sys.time() #Establish the time when running the function
  
  if(year==F){
    marres=levels(x$Id_Sites) #Assign in an object with the name of the site
    
    
    pb <- txtProgressBar(min = 0, max = length(marres), style = 3) #Create the progress bar
    
    temp_count=matrix(nrow = length(marres),ncol = length(levels(x$Taxon))+1)
    colnames(temp_count)=c(as.character(as.data.frame(table(x[x[,2],7]))[,1]),"Site")
    
    descriptors=as.data.frame(matrix(nrow=length(levels(x$Id_Sites)),ncol=18)) #Create an empty data frame for the descriptors
    colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level",
                            "Mean_Omnivory_Index", "Richesse_spe","Links_per_species","indegree",
                            "outdegree","Percentage_basal","Percentage_top","Percentage_intermediate"
                            ,"Omnivory","S_Top",
                            "S_intermediate","S_Basal") #Specified the name of the column for the descriptos object
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
      
      
      
      level=as.character(c(levels(temp_site$Taxon))) #Take species species column name of species present in ponds
      
      temp=y[level,level] #subset the regionnal network with only present species
      
      if(descriptor==T){ # If arg "descriptor" is true, calculate the descriptor
        g=graph_from_adjacency_matrix(temp, mode="undirected")
        descriptors[o,1]=name_site[o]
        name_matrix[o]=paste("Matrix_site_",name_site[o],sep="")
        descriptors[o,2:4]=enaStructure(temp)$ns[,c(1,2,3)]
        descriptors[o,5]=modularity(g,membership(cluster_walktrap(g)))
        descriptors[o,6]=mean(TrophInd(temp)[,1])
        descriptors[o,7]=mean(TrophInd(temp)[,2])
        descriptors[o,8]=length(temp[,1])
        descriptors[o,9]=sum(temp)/length(temp[,1]) # Number of links per species
        descriptors[o,10]= MeanGenerality(temp) #indegree
        descriptors[o,11]= MeanVulnerability(temp) #outdegree
        descriptors[o,12]= FractionOfBasal(temp) # Percentage of basal species
        descriptors[o,13]= FractionOfTop(temp) # Percentage of top species
        descriptors[o,14]= FractionOfIntermediate(temp) # Percentage of intermediate species
        descriptors[o,15]= Omnivory(temp)
        descriptors[o,16]= NumberOfTop(temp)# Number of top species
        descriptors[o,17]= NumberOfIntermediate(temp) # Number of intermediate species
        descriptors[o,18]= NumberOfBasal(temp) # Number of basal species
      }
      if(Matrix==T){assign(nm,temp, env = .GlobalEnv)} # If arg "Matrix" is true, put in global environnement adjancy matrix of the current pond
    }
    
    descriptors$Site=as.factor(descriptors$Site)
    
    if(descriptor==T & Count==F){assign(paste(name_tab),
                                        merge(descriptors,Surf,by.x = c("Site"),by.y = c("Id_Mare")), env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
    col_temp=temp_count[,length(levels(x$Taxon))+1]
    temp_count[temp_count>0]=1
    temp_count[,length(levels(x$Taxon))+1]=col_temp
    if(Count==T){assign(paste(name_tab,"count",sep=""),
                        merge(merge(descriptors,Surf,by.x = c("Site"),by.y = c("Id_Mare")),temp_count,by.x = c("Site"),by.y = c("Site")), env = .GlobalEnv)}
  }
  
  if(year==T){
    marres=levels(x$Id_Sites) #Assign in an object with the name of the site
    years=levels(x$Annee)
    
    pb <- txtProgressBar(min = 0, max = length(marres), style = 3) #Create the progress bar
    
    temp_count=matrix(nrow = length(marres),ncol = length(levels(x$Taxon))+1)
    colnames(temp_count)=c(as.character(as.data.frame(table(tab_releve[tab_releve[,2],7]))[,1]),"Site")
    
    count_year=as.matrix(table(x[,c(2,3)]))[,1:6]
    count_year[count_year>0]=1
    
    
    
    descriptors=as.data.frame(matrix(nrow=sum(count_year),ncol=19)) #Create an empty data frame for the descriptors
    colnames(descriptors)=c("Site","n","L","C","Modularity",
                            "Mean_trophic_level","Mean_Omnivory_Index", 
                            "Richesse_spe","Annee","Links_per_species","indegree",
                            "outdegree","Percentage_basal","Percentage_top","Percentage_intermediate"
                            ,"Omnivory","S_Top",
                            "S_intermediate","S_Basal") #Specified the name of the column for the descriptos object
    name_site=levels(x$Id_Sites) #Assign in an object with the name of the site
    name_matrix=1:length(levels(x$Id_Sites)) #Assign a numeric level for the name of matrix
    year_matrix=1:length(levels(x$Annee))
    
    for (o in 1:length(marres)) {#Loop through each ponds
      Sys.sleep(0.1) # Define refresh rate for progress bar
      setTxtProgressBar(pb, o) #Set the current status of the progress bar
      
      for(h in 1:length(levels(x$Annee))){
        if(years[h] %in% x[x[,2]==marres[o],3]){
          
          temp_count[o,1:length(levels(x$Taxon))]=as.data.frame(table(x[x[,2]==marres[o] & x[,3]==years[h],7]))[,2]
          temp_count[o,length(levels(x$Taxon))+1]=name_site[o]
          
          temp_site=x[x[,2]==marres[o] & x[,3]==years[h],] #Subset the taxon occurrence data for a define pond
          temp_site[,2]=droplevels(temp_site[,2]) #Forget old levels
          temp_site[,7]=droplevels(temp_site[,7]) #Forget old levels
          
          nm=paste("Matrix_site_",marres[o],"_",years[h],sep="") #Define the name of the matrix returned for the current pond
          
          
          level=as.character(c(levels(temp_site$Taxon)))
          
          
          if(length(level)==1){level=as.character(c(levels(temp_site$Taxon),"Zooplankton","Phytoplankton","Detritus","Macrophyte"))}
          
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
            descriptors[o,9]=years[h]
            descriptors[o,10]=sum(M)/length(temp[,1]) # Number of links per species
            descriptors[o,11]= MeanGenerality(temp) #indegree
            descriptors[o,12]= MeanVulnerability(temp) #outdegree
            descriptors[o,13]= FractionOfBasal(temp) # Percentage of basal species
            descriptors[o,14]= FractionOfTop(temp) # Percentage of top species
            descriptors[o,15]= FractionOfIntermediate(temp) # Percentage of intermediate species
            descriptors[o,16]= Omnivory(temp)
            descriptors[o,17]= NumberOfTop(temp)# Number of top species
            descriptors[o,18]= NumberOfIntermediate(temp) # Number of intermediate species
            descriptors[o,19]= NumberOfBasal(temp) # Number of basal species
          }
          if(Matrix==T){assign(nm,temp, env = .GlobalEnv)} # If arg "Matrix" is true, put in global environnement adjancy matrix of the current pond
        }
        
        
        if(descriptor==T & Count==F){assign(paste(name_tab),
                                            merge(descriptors,envi_data,by.x = c("Site"),by.y = c("Mares"),all.x = T), env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
        col_temp=temp_count[,length(levels(x$Taxon))+1]
        temp_count[temp_count>0]=1
        temp_count[,length(levels(x$Taxon))+1]=col_temp
        if(Count==T){assign(paste(name_tab,"count",sep = ""),
                            merge(merge(descriptors,Surf,by.x = c("Site"),by.y = c("Id_Mare")),temp_count,by.x = c("Site"),by.y = c("Site")), env = .GlobalEnv)}
        
      }
    }
  }
  
  
  
  close(pb) # Close the progress bar
  end.time <- Sys.time() #Establish the time when ending the function
  time.taken <- end.time - start.time # Calculate the time the function worked
  time.taken # Return the time taking for the function to work
  
}

##########


data=datacount


########################## RELATION AMONG DESCRIPTORS
# Environmental
pairs(data[,c(13:15,17:20,25)])
# Topological
pairs(data[,c(2:7)])

cor_envi=cor(data[,c(13:15,17:20,25)])
cor_envi>=0.4 | cor_envi<=-0.4

cor_topo=cor(data[,c(2:7)])
cor_topo>=0.4 | cor_topo<=-0.4

#########Connectance

m0 <- lm(C ~ scale(log(Surface)) 
         , data = data)

anova(m0)

summary(m0)

write.table(summary(m0)$coefficients,"clipboard",sep="\t")

#########SPecies Richness

m1=glm(n ~ scale(log(Surface))+scale(I(log(Surface)^2)) 
      , data = data, family = quasipoisson(link = "log"))

anova(m1, test="Chisq")

summary(m1)

write.table(summary(m1)$coefficients,"clipboard",sep="\t")

#########Modularity
m2=lm(Modularity ~ scale(log(Surface))+scale(I(log(Surface)^2))
      , data = data)

anova(m2, test="Chisq")

summary(m2)

write.table(summary(m2)$coefficients,"clipboard",sep="\t")

m2=gam(n ~ s(scale(log(Surface)))
       , data = data)


summary(m2)

plot(m2)

#########Mean Trophic Level

m3=lm(Mean_trophic_level ~ scale(log(Surface))+scale(I(log(Surface)^2))
      , data = data)

anova(m3, test="Chisq")

summary(m3)

write.table(summary(m3)$coefficients,"clipboard",sep="\t")

m3=gam(Mean_trophic_level ~ s(scale(log(Surface)))
       , data = data)


summary(m3)

plot(m3)



#########Mean_Omnivory_Index

m4=lm(Mean_Omnivory_Index ~ scale(log(Surface))+scale(I(log(Surface)^2))
      , data = data)

anova(m4, test="Chisq")

  summary(m4)

write.table(summary(m4)$coefficients,"clipboard",sep="\t")

m4=gam(Mean_Omnivory_Index ~ s(scale(log(Surface)))
    , data = data)

summary(m4)

plot(m4)


#####Surface

p4=ggplot(datacount, aes(x=log(Surface), y=n)) +
  geom_point(size = 3, alpha=0.8) + geom_smooth(method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
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
  labs(x = "", y="Species richness")







p6=ggplot(datacount, aes(x=log(Surface), y=Mean_trophic_level)) +
  geom_point( size = 3, alpha=0.8) + geom_smooth( method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
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
  labs(x = "", y="Mean trophic level")




p8=ggplot(datacount, aes(x=log(Surface), y=Mean_Omnivory_Index)) +
  geom_point( size = 3, alpha=0.8) + geom_smooth( method = "lm",formula = y ~ x + I(x^2), se=F,size=2)+
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
  labs(x = "", y="Mean Omnivory Index")


p2=ggplot(datacount, aes(x=log(Surface), y=Modularity)) +
  geom_point( size = 3, alpha=0.8) + 
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
  labs(x = "log(size area)",y="Modularity")






p10=ggplot(datacount, aes(x=log(Surface), y=C)) +
  geom_point( size = 3, alpha=0.8) +
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
  labs(x = "", y="Connectance")



########*


ggarrange(p4,p6,p8,p10,p2,ncol=1,nrow=5,labels=c("A"),font.label=c(size=32), vjust=1
          ,legend=c("bottom", "right"),common.legend = T)




