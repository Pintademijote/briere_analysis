

library(hash)
library(igraph)
library(cheddar)
library(sna)


Generality <- function(M){
  return(sum(colSums(M))/sum((colSums(M)!=0)));
}

Vulnerability <- function(M){
  return(sum(rowSums(M))/sum((rowSums(M)!=0)));
}


InDegree <- TrophicGenerality <- NumberOfResources <- function(M){
  return(colSums(M));
}

OutDegree <- TrophicVulnerability <- NumberOfCosumers <- function(M){
  return(rowSums(M));
}


MeanGenerality <- function(M){
  return(mean(colSums(M)));
}

MeanVulnerability <- function(M){
  return(mean(rowSums(M)));
}


FractionOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
  
  return(b_sps / dim(M)[1]);
}

NumberOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
  
  return(b_sps);
}

FractionOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
  
  return(t_sps / dim(M)[1]);
}

NumberOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
  
  return(t_sps);
}

FractionOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
  
  return(i_sps / dim(M)[1]);
}

NumberOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
  
  return(i_sps);
}

FractionOfCannibalism <- function(M){
  return(sum(diag(M)) / dim(M)[1]);
}


Omnivory <- function(M){
  #Use PredationMatrixToLinks() to create a Cheddar community from a predation
  # matrix
  M <- t(M)
  node <- 1:dim(M)[1];
  for(n in 1:length(node)){
    node[n] <- paste(node[n],'-');
  }
  
  pm <- matrix(M, ncol=dim(M)[2], dimnames=list(node, node), byrow=TRUE);
  
  community <- Community(nodes=data.frame(node=node), trophic.links=PredationMatrixToLinks(pm), properties=list(title='Community'));
  
  community <- RemoveCannibalisticLinks(community, title='community');
  
  #community is a cheddar community
  
  Fractionomnivory <- FractionOmnivorous(community)
  
  return(Fractionomnivory)
}

MeanFoodChainLength <- function(M){
  
  # Use PredationMatrixToLinks() to create a Cheddar community from a predation
  # matrix
  M <- t(M)
  node <- 1:dim(M)[1];
  for(n in 1:length(node)){
    node[n] <- paste(node[n],'-');
  }
  
  pm <- matrix(M, ncol=dim(M)[2], dimnames=list(node, node), byrow=TRUE);
  
  community <- Community(nodes=data.frame(node=node), trophic.links=PredationMatrixToLinks(pm), properties=list(title='Community'));
  
  community <- RemoveCannibalisticLinks(community, title='community');
  
  chain.stats <- TrophicChainsStats(community)
  ch_lens <- (chain.stats$chain.lengths + 1)
  
  return(sum(ch_lens)/length(ch_lens));
}






