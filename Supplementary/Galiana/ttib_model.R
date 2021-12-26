

setwd('C:/Users/pglemasle/Downloads/Code_Galiana_EtAl.2018')
source('utils_NAR.R') #This code is used to compute all the network metrics


#Here we define the colonisation and extinction parameters. 
#The ratio between c and e (c/e) determine the area size for this model.
c <- 0.2   #colonisation rate
extinction <- seq(0.5,0.95,0.01) #extinction rate


timesteps <- 1000 # number of time steps used in the colonisation-extinction dynamics
replicates <- 1 # Number of replicates desired for the simulations. In the paper we used 100.

output <- NULL


for(e in extinction){
  for(j in 1:replicates){
    load(paste('./networks/network-',j,sep='')) # This loads the networks previously generated with the niche model that will be used as the metaweb
    # For this to work you need to save the networks generated in a folder called networks within your working directory.
    S <- food_web$S
    presences <- array(0, S)
    
    print(j)
    for(i in 1:timesteps){ # This loop runs the colonisation-extinction dynamics with the trophic constraint of the TTIB
      #basal species
      basal_sps <- which(colSums(food_web$M) == 0)
      cur_com <- presences
      
      #this is what species are present
      present <- which(cur_com != 0)
      
      if(sum(cur_com) != 0){
        
        #this is the set of species in the network that have prey in the current community
        with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
        #this are the potential colonisers
        potential <- append(basal_sps , with_prey)
        #to which we remove the ones that are already present
        potential <- setdiff(potential, present)
        
      }else{
        potential <- basal_sps
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
      
      if(extinctions){
        cur_com <- presences
        #this is what species are present
        present <- which(cur_com != 0)
        #this is the set of species in the network that have prey in the current community
        with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
        legal <- union(with_prey, basal_sps)
        illegal <- present[!(present %in% legal)]
        presences[illegal] <- 0  
      }
      
    }
    
    results <- presences
    net <- food_web
    
    current_species <- which(results > 0) # Number of species
    
    if(length(current_species) < 2){
      next
    }
    
    M <- net$M[current_species, current_species] # This is the assembled community after the colonization-extinction dynamics
    
    #we can check if the network is connected
    graf <- igraph::graph.adjacency(M); # Making sure the network is connected
    y <- igraph::no.clusters(graf, mode='weak');
    attempts <- 0
        if(y > 1 ){
          print('network not found')
          next
        }
    # Now we compute all the network metrics using the code in utils_NAR.R
    C <- sum(M)/((dim(M)[1])*(dim(M)[1]-1)) #Connectance if we prevent cannibalism
    
    mfcl <- tryCatch({
      MeanFoodChainLength(M)
    }, warning = function(w) {
      'NA'
    }, error = function(e) {
      'NA'
    }, finally = {

    })

   mfcls <- mfcl # Mean food chain length
   links <- sum(M) # Number of links
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
  

    cur_out <- data.frame(replicate=j, S=length(current_species), C, links, links.species, 
                          indegree, outdegree, basal, top, intermediate, omnivory, S_top, 
                          S_intermediate, S_basal, ratio, colonization_rate, extinction_rate)


    if(is.null(output)){
      output <- cur_out
    }else{
      output <- rbind(output, cur_out)
    } 
    
  }
}

write.csv(output, file='output_ttib_test.csv') # Save the output


### Plotting Network-Area relationships; In the x-axis we plot the ratio because as said
### before in this model the ratio between c and e is used as a proxy of area. 

plot(output$ratio, output$S)
