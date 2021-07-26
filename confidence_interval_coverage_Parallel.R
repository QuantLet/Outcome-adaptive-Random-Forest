library(ranger)
library(scales)
library(MASS)
library(ggplot2)
library(cowplot)
library(glmnet)
library(RRF)


install.packages(vec.pac)
vec.pac= c("mvtnorm","clusterGeneration",
           "ranger","scales","MASS","glmnet","RRF","lqa","foreach","doParallel")
lapply(vec.pac, require, character.only = TRUE) 






cl   <- makeCluster(5, outfile="")
registerDoParallel(cl)



# Generate new data.

# Setting works quite ok with p <- 50
# Setting works quite ok with n <- 400 

B <- 250






r <- foreach(sim = 1:B, .combine='rbind', .inorder=FALSE, .packages=vec.pac, .errorhandling ="remove" ) %dopar% { 
  
  coverage <- matrix(NA,1,5)
  colnames(coverage) <- c("RF_full", "OARF", "RRF", "OAL", "IPTW")
  
  
  DGP <- datagen(setting=1, n=1000,p=20)
  data <- DGP[[1]]
  var.list <- DGP[[2]]
  bA <- 0.5
  
  # Normalize covariates to have mean 0 and standard deviation 1
  temp.mean = colMeans(data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(data),byrow=TRUE)
  data[,var.list] = data[,var.list] - Temp.mean
  temp.sd = apply(data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(data),byrow=TRUE)
  data[var.list] = data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  ####################################################################
  
  
  
  ###############
  # Bootstrap iterations
  nboot <- 500
  
  # RF, OARF and RRF
  RF_ci <- all_RF_boot(data=data,verbose=FALSE)
  
  for(method in 1:3){
    coverage[1,method]  <- truth_in_ci(truth=bA,ci=RF_ci[[1]][,method])
  }
  # OAL 
  OAL_ci <- shortreed_boot(data=data,nboot=nboot)
  coverage[1,4]  <- truth_in_ci(truth=bA,ci=OAL_ci[[1]])
  # IPTW
  IPTW_ci <- iptw_boot_ci(data=data,nboot=nboot)
  coverage[1,5]  <- truth_in_ci(truth=bA,ci=IPTW_ci)
  
  #cat("This is iteration", sim, "out of", B, "\n")

  return(coverage)
}


colMeans(r)
