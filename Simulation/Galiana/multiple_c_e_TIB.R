require("enaR")
require("igraph")
require("NetIndices")

setwd('C:/Users/pglemasle/Downloads/Code_Galiana_EtAl.2018')
source('utils_NAR.R') #This code is used to compute all the network metrics

Matrix_Reg_Manual=read.delim2("S:/2018/pglemasle/Data/Script/Matrix_Manual_Chekck.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)

#Here we define the colonisation and extinction parameters. 
#The ratio between c and e (c/e) determine the area size for this model.
colonization <- seq(0.1,0.99,0.025)   #colonisation rate
extinction <- seq(0.1,0.99,0.025) #extinction rate



timesteps <- 1000 # number of time steps used in the colonisation-extinction dynamics
replicates <- 1 # Number of replicates desired for the simulations. In the paper we used 100.

output <- NULL


for (rep in 1:100) {
  pb <- txtProgressBar(min = 0, max = length(colonization), style = 3)
  repet=0
  start.time = Sys.time()
  for(c in colonization){
    Sys.sleep(0.1)
    repet=repet+1
    setTxtProgressBar(pb, repet)
    for(e in extinction){
      for(j in 1:replicates){
        food_web=Matrix_Reg_Manual # This loads the networks previously generated with the niche model that will be used as the metaweb
        # For this to work you need to save the networks generated in a folder called networks within your working directory.
        
        S <- length(Matrix_Reg_Manual[,1])
        presences <- array(0, S)
        presences_mainland <- array(1, S)
        
        
        for(i in 1:timesteps){ # This loop runs the colonisation-extinction dynamics with the trophic constraint of the TTIB
          #basal species
          basal_sps <- which(colSums(food_web) == 0)
          cur_com <- presences
          
          #this is what species are present
          present <- which(cur_com != 0)
          
          if(sum(cur_com) != 0){
            
            potential <- !(presences %in% presences_mainland)
            
          }else{
            potential <- !(presences %in% presences_mainland)
          }
          
          if(length(potential) != 0){
            colonisers <- potential[which(runif(length(potential), 0, 1) < c)]
            presences[colonisers] <- 1
          }
          
          extinctions <- FALSE
          if(length(present) != 0){
            extinct <- present[which(runif(length(present), 0, 1) < e)]
            if(length(extinct) != 0){
              presences[extinct] <- 0
              extinctions <- TRUE
            }
          }
          
          
          
        }
        
        results <- presences
        net <- food_web
        
        current_species <- which(results > 0) # Number of species
        
        if(length(current_species) < 2){
          next
        }
        
        M <- food_web[current_species, current_species] # This is the assembled community after the colonization-extinction dynamics
        

        # Now we compute all the network metrics using the code in utils_NAR.R
        
        g=graph_from_adjacency_matrix(M, mode="undirected")
        
        C <- sum(M)/((dim(M)[1])*(dim(M)[1]-1)) #Connectance if we prevent cannibalism
        mfcl <- 0
        mfcls <- mfcl # Mean food chain length
        links <- sum(M) # Number of links
        modularity_t=modularity(g,membership(cluster_walktrap(g)))
        Mean_Trophic_Level=mean(TrophInd(M)[,1])
        Mean_Omnivory=mean(TrophInd(M)[,2])
        links.species <- sum(M)/length(current_species) # Number of links per species
        indegree <- MeanGenerality(M)
        outdegree <- MeanVulnerability(M)
        basal <- FractionOfBasal(M) # Percentage of basal species
        top <- FractionOfTop(M) # Percentage of top species
        intermediate <- FractionOfIntermediate(M) # Percentage of intermediate species
        omnivory <- Omnivory(M)
        S_top <- NumberOfTop(M)# Number of top species
        S_intermediate <- NumberOfIntermediate(M) # Number of intermediate species
        S_basal <- NumberOfBasal(M) # Number of basal species
        ratio <- c/e # Area
        
        colonization_rate <- c
        extinction_rate <- e
        
        
        cur_out <- data.frame(replicate=rep, S=length(current_species), C, links,modularity_t,
                              Mean_Trophic_Level,Mean_Omnivory,links.species, 
                              indegree, outdegree, basal, top, intermediate, omnivory, S_top, 
                              S_intermediate, S_basal, ratio, colonization_rate, extinction_rate)
        
        
        if(is.null(output)){
          output <- cur_out
        }else{
          output <- rbind(output, cur_out)
        } 
        
      }
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time # Calculate the time the function worked
  print(time.taken)
}

write.csv(output, file='output_tib_test.csv') # Save the output


### Plotting Network-Area relationships; In the x-axis we plot the ratio because as said
### before in this model the ratio between c and e is used as a proxy of area. 

plot(output$extinction_rate, output$S)

plot(output$ratio, output$basal)

plot(output$ratio, output$omnivory)


