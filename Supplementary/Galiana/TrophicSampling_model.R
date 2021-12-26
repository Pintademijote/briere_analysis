setwd('./Code_Galiana_EtAl.2018')

source('utils_NAR.R')


## We define the power function describing the species area relationship
sar <- function(area){
  return(10*(area**0.27))
}

connectance <- function(S){
  return(0.8*(S**(-0.5)))
}



sps <- seq(20,200,1) # This determines the minimun and the maximum number of species we want to have in the local networks
c_s <- connectance(sps)

simulations <- 1 # Number of replicates desired for the simulations. In the paper we used 100.


output <- NULL

for(j in 1:simulations){ 
  print(j)
  load(paste('./networks/network-',j,sep='')) # This loads the networks previously generated with the niche model that will be used as the metaweb.
  # For this to work you need to save the networks generated in a folder called networks within your working directory.
  
  n_regional <- food_web #this is the regional network
  
  species <- c()
  area <- c()
  C <- c()
  links <- c()
  links.species <- c()
  mfcls <- c()
  indegree <- c()
  outdegree <- c()
  basal <- c()
  top <- c()
  intermediate <- c()
  omnivory <- c()
  S_top <- c()
  S_intermediate <- c()
  S_basal <- c()

  for(i in 1:length(sps)){ #this loop is to generate the local networks from the regional one based on SAR and C
    #this networks are the ones we want to analise. 
    print(paste('species',sps[i])) 
    
    connectance<-c_s[i] #Connectance associated to the given number of species 
    cols <- sample(ncol(n_regional$M), sps[i]) # We pick randomly the number of species we want to create the local network from the regional one
   
    n <- n_regional$M[cols,cols] #this is the local network

    graf <- igraph::graph.adjacency(n); # Making sure the network is connected
    y <- igraph::no.clusters(graf, mode='weak');
    while(y > 1){
      cols <- sample(ncol(n_regional$M), sps[i]) # We pick randomly the number of species we want to create the local network from the regional one
      n <- n_regional$M[cols,cols] #this is the local network
      graf <- igraph::graph.adjacency(n); # Making sure the network is connected
      y <- igraph::no.clusters(graf, mode='weak');
    }
    
# Here we compute all network metrics using the utils_NAR.R code
   C <- append(C,  sum(n)/(dim(n)[1])^2) # Connectance
   
   mfcl <- tryCatch({  
     MeanFoodChainLength(n)
   }, warning = function(w) {
     'NA'
   }, error = function(e) {
     'NA'
   }, finally = {
     
   })
   
   species <- append(species, dim(n)[1]) # Number of species
   area <- append(area, round(exp(log(sps[i]/10)/0.27))) 
   links <- append(links, sum(n)) # Number of links
   links.species <- append(links.species, sum(n)/dim(n)[1]) # Number of links per species
   mfcls <- append(mfcls, mfcl) # Mean food chain length
   indegree <- append(indegree, MeanGenerality(n))
   outdegree <- append(outdegree, MeanVulnerability(n))
   basal <- append(basal, FractionOfBasal(n)) # Percentage of basal species
   top <- append(top, FractionOfTop(n)) # Percentage of top species
   intermediate <- append(intermediate, FractionOfIntermediate(n)) # Percentage of intermediate species
   omnivory <- append(omnivory, Omnivory(n))
   S_basal <- append(S_basal, NumberOfBasal(n)) # Number of basal species
   S_top <- append(S_top, NumberOfTop(n)) # Number of top species
   S_intermediate <- append(S_intermediate, NumberOfIntermediate(n)) # Number of intermediate species
  }
  
  cur_out <- data.frame(rep(i,length(species)), species, area, C, links, links.species, mfcls, indegree, outdegree, basal, top, intermediate, omnivory, S_intermediate, S_top, S_basal)

  if(is.null(output)){
    output <- cur_out
  }else{
    output <- rbind(output, cur_out)
  }
   
}

write.csv(output, file='output_trophicSampling.csv')

### Plotting Network-Area relationships; 
### In this model area is determined by the number of species following the power function defined for the Species-Area Relationship.

plot(output$area, output$species)
