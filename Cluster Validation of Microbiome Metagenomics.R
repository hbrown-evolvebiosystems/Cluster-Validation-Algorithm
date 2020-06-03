#load libraries
library(readr)
library(readxl)
library(tidyverse)
library(cluster)
library(clusterSim)
library(ade4)
library(matlib)
library(data.table)


#set working directory  
setwd("~/SCOOP/Enterotype Validation")

#import data

##species metagenomics for all samples (293 samples over 77 subjects)
SCOOP_species_overall <- read_csv("SCOOP_species.csv")

##metadata
antibiotic_last_7_days <- read_csv("antibiotic_last_7_days.csv")
Corrected_Evivo_Ever <- read_csv("Corrected_Evivo_Ever.csv")

##1000 rows of 1 randomly selected sample per subject (77 total subjects)
Random_Sampling <- read_csv("Random_Sampling.csv")

#Enterotyping functions

##distance matrix function
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

##partitioning around mediods clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

##noise removal
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

#################################################
#Start Simulation

#subsetting species dataset
for(sub in c(0)){
  summary_data = data.frame(matrix(nrow = 0, ncol = 10))
  colnames(summary_data) = c("X1", "X2", "X3", "X4", "X5", "X6", "length", "species", "nclusters", "sim")
  
  SCOOP_species = SCOOP_species_overall
  
  obs.silhouette = NULL
  mean.obs.silhouette = NULL
  CH_vals = NULL
  
  for (iter in 1:1000){ 
    nam <- paste("obs.bet", iter, sep = "")
    
    
    
    df = SCOOP_species[SCOOP_species$SampleID %in% Random_Sampling[iter,],]
    
    df = column_to_rownames(df, var = "SampleID") %>% t(.) %>% data.frame(.) %>% rownames_to_column(var = "#OTU ID")
    df$species = gsub(".*s__", "", df$`#OTU ID`)
    df = df[-1]
    
    row.names(df) = df$species
    
    #df = df[df$species %in% filter(species_summary_overall, `n()` >= sub)$species,]
    
    
    species.dist = dist.JSD(df[,-78])
    
    #species.cluster = pam.clustering(species.dist, k = 3)
    
    #nclusters = index.G1(t(test_species[,-78]), species.cluster, d = species.dist, centrotypes = "medoids")
    
    nclusters=NULL
    
    for (k in 2:20) { 
      if (k==1) {
        nclusters[k]=NA 
      } else {
        data.cluster_temp=pam.clustering(species.dist, k)
        nclusters[k]=index.G1(t(df[,-78]),data.cluster_temp,  d = species.dist,
                              centrotypes = "medoids")
      }
    }
    
    #plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
    
    CH_vals[[paste("iter", iter, sep = "")]] <- nclusters
    
    species.cluster=pam.clustering(species.dist, k= which.max(nclusters))
    
    mean.obs.silhouette[iter]=mean(silhouette(species.cluster, species.dist)[,3])
    
    obs.silhouette[[paste("iter",iter, sep = "")]] <- silhouette(species.cluster, species.dist)[,3]
    
    data.denoized=noise.removal(df[,-78], percent=0.01)
    
    obs.pca=dudi.pca(data.frame(t(df[,-78])), scannf=F, nf=10)
    rownames(obs.pca$co) = df$species
    obs.bet=bca(obs.pca, fac=as.factor(species.cluster), scannf=F, nf=k-1)
    rownames(obs.bet$co) = df$species
    
    assign(nam, obs.bet)
    saveRDS(object = get(x = nam),file = paste(sub,"obs.bet", iter, ".RDS", sep = ""))
    
    obs.bet$co$length <- if(ncol(obs.bet$co) == 1) abs(obs.bet$co$Comp1) else sqrt(obs.bet$co$Comp1^2 + obs.bet$co$Comp2^2)
    
    top_nclusterx2_species = obs.bet$co %>%
      rownames_to_column(var = "species") %>%
      arrange(desc(length)) %>%
      top_n(which.max(nclusters)*2) %>%
      column_to_rownames(var = "species")
    
    t = data.frame(matrix(nrow = nrow(top_nclusterx2_species), ncol = 24), row.names = row.names(top_nclusterx2_species))
    colnames(t) = c(paste("X",1:20, sep = ""), "length", "species", "nclusters", "sim")
    
    
    for (i in 1:nrow(obs.bet$li)) {
      for (j in 1:nrow(top_nclusterx2_species)) {
        t[j,i] = ifelse(angle(as.numeric(obs.bet$li[i,]), as.numeric(top_nclusterx2_species[j,-ncol(top_nclusterx2_species)])) < 22.5,
                        angle(as.numeric(obs.bet$li[i,]), as.numeric(top_nclusterx2_species[j,-ncol(top_nclusterx2_species)])), NA)    
        
      }
    }
    
    t$species = row.names(t)
    t$length = top_nclusterx2_species$length
    t$nclusters = which.max(nclusters)
    t$sim = iter
    rownames(t) = c()
    
    summary_data = rbind(summary_data, t)
    
    write.csv(summary_data, paste("summary_data", sub, ".csv", sep = ""), row.names = FALSE)
    
    
    rm(list = paste("obs.bet",iter,sep = ""))
  }
  
  obs.silhouette = data.frame(obs.silhouette)
  write.csv(obs.silhouette, paste("silhouette_data", sub, ".csv", sep = ""), row.names = FALSE)
  
  CH_vals = data.frame(CH_vals)
  write.csv(CH_vals, paste("CH_vals", sub, ".csv", sep = ""), row.names = FALSE)
  
  
}

cluster_summary = summary_data %>% distinct(sim, nclusters)
species_summary = summary_data %>% group_by(species) %>% summarize(n())
#write.csv(species_summary, "species_summary.csv")
species_cluster_summary = summary_data %>% group_by(species, nclusters) %>% summarize(n())


###################################
#silhouette obs entire dataset
#df = SCOOP_species[SCOOP_species$SampleID %in% Random_Sampling[iter,],]
df = SCOOP_species
df = column_to_rownames(df, var = "SampleID") %>% t(.) %>% data.frame(.) %>% rownames_to_column(var = "#OTU ID")
df$species = gsub(".*s__", "", df$`#OTU ID`)
df = df[-1]

df = df[df$species %in% filter(species_summary_overall, `n()` >= 800)$species,]

row.names(df) = df$species

species.dist = dist.JSD(df[,-294])

nclusters=NULL

for (k in 2:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(species.dist, k)
    nclusters[k]=mean(silhouette(data.cluster_temp, species.dist)[,3])
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="Average Silhouette Value")

#species.cluster=pam.clustering(species.dist, k= which.max(nclusters))

#mean.obs.silhouette[iter]=mean(silhouette(species.cluster, species.dist)[,3])

#obs.silhouette[[paste("iter",iter, sep = "")]] <- silhouette(species.cluster, species.dist)[,3]
#CH index
nclusters=NULL

for (k in 2:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(species.dist, k)
    nclusters[k]=index.G1(t(df[,-294]),data.cluster_temp,  d = species.dist,
                          centrotypes = "medoids")
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

#######################################################################
#silhouette obs subset
df = SCOOP_species[SCOOP_species$SampleID %in% Random_Sampling[iter,],]
df = SCOOP_species
df = column_to_rownames(df, var = "SampleID") %>% t(.) %>% data.frame(.) %>% rownames_to_column(var = "#OTU ID")
df$species = gsub(".*s__", "", df$`#OTU ID`)
df = df[-1]

row.names(df) = df$species

species.dist = dist.JSD(df[,-294])

nclusters=NULL

for (k in 2:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(species.dist, k)
    nclusters[k]=mean(silhouette(data.cluster_temp, species.dist)[,3])
  }
}

###################################################
#sampling/cross-validation
df = SCOOP_species_overall
df = column_to_rownames(df, var = "SampleID") %>% t(.) %>% data.frame(.) %>% rownames_to_column(var = "#OTU ID")
df$species = gsub(".*s__", "", df$`#OTU ID`)
df = df[-1]

df = df[df$species %in% filter(species_summary_overall, `n()` >= 800)$species,]

row.names(df) = df$species
df = df[,-ncol(df)]
pct_0 = NULL
mean.obs.silhouette = NULL
summary_data = NULL
cluster_validation = NULL

for(iter in 1:1000) {
  
  train_index = sample(1:ncol(df), 0.8*ncol(df))
  df_train = df[,train_index]
  pct_0[iter] = mean(rowSums(df_train == 0)/ncol(df_train))
  
  species.dist = dist.JSD(df_train)
  
  species.cluster=pam.clustering(species.dist, k= 3)
  
  mean.obs.silhouette[iter]=mean(silhouette(species.cluster, species.dist)[,3])
  
  obs.pca=dudi.pca(data.frame(t(df_train)), scannf=F, nf=10)
  rownames(obs.pca$co) = df$species
  obs.bet=bca(obs.pca, fac=as.factor(species.cluster), scannf=F, nf=3-1)
  rownames(obs.bet$co) = rownames(df_train)
  
  obs.bet$co$length <- sqrt(obs.bet$co$Comp1^2 + obs.bet$co$Comp2^2)
  
  top_nclusterx2_species = obs.bet$co %>%
    rownames_to_column(var = "species") %>%
    arrange(desc(length)) %>%
    #top_n(which.max(nclusters)*2) %>%
    column_to_rownames(var = "species")
  
  t = data.frame(matrix(nrow = 4, ncol = 7), row.names = row.names(top_nclusterx2_species))
  colnames(t) = c(paste("X",1:3, sep = ""), "length", "species", "nclusters", "sim")
  
  
  for (i in 1:nrow(obs.bet$li)) {
    for (j in 1:nrow(obs.bet$co)) {
      t[j,i] = ifelse(angle(as.numeric(obs.bet$li[i,]), as.numeric(top_nclusterx2_species[j,-ncol(top_nclusterx2_species)])) < 22.5,
                      angle(as.numeric(obs.bet$li[i,]), as.numeric(top_nclusterx2_species[j,-ncol(top_nclusterx2_species)])), NA)    
      
    }
  }
  
  t$species = row.names(t)
  t$length = top_nclusterx2_species$length
  t$nclusters = 3
  t$sim = iter
  rownames(t) = c()
  
  summary_data = rbind(summary_data, t)
  write.csv(summary_data, paste("summary_data_3_clusters", ".csv", sep = ""), row.names = FALSE)
  
  
  cluster_assignment = obs.bet$ls
  cluster_assignment$SampleID = rownames(cluster_assignment)
  cluster_assignment$cluster = species.cluster
  cluster_assignment$sim  = iter
  rownames(cluster_assignment) = c()
  
  cluster_validation = rbind(cluster_validation, cluster_assignment)
  write.csv(cluster_validation, paste("cluster_validation", ".csv", sep = ""), row.names = FALSE)
  
}

cluster_names = summary_data %>% gather(cluster,angle,1:3) %>% 
  mutate(cluster = as.numeric(substr(cluster,2,2))) %>% drop_na() %>%
  group_by(sim, cluster) %>% mutate(cluster_name = paste0(species, collapse = " & ")) %>%
  distinct(sim, cluster, cluster_name)

cluster_mapping = merge(cluster_validation[,-c(1:2)], cluster_names, by = c("sim", "cluster"))


cluster_agreement = cluster_mapping %>% group_by(SampleID, cluster_name) %>% 
  summarize(n = n()) %>% mutate(N = sum(n),freq = n/sum(n))

mean(cluster_agreement$N)

accuracy = round(1 - sum(filter(cluster_agreement, freq < 0.8)$n)/sum(filter(cluster_agreement, freq > 0.8)$N),4)*100
tot = sum(filter(cluster_agreement, freq > 0.8)$N)
miss = sum(filter(cluster_agreement, freq < 0.8)$n)

validated_cluster_assignments  =  filter(cluster_agreement, freq > 0.8)  %>% distinct(SampleID,cluster_name)
write.csv(validated_cluster_assignments, paste("validated_cluster_assignments", ".csv", sep = ""), row.names = FALSE)

##################################
