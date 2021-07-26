### This generates the Simulation Data as a dataframe
# by Daniel Jacob (daniel.jacob@hu-berlin.de) 

# Arguments to specify are: 

# N = Number of observations (real number)
# k = Number of covariates (real number)
# random_d = treatment assignment: (Either T for random assignment or F for confounding on X)
# theta = treatment effect: (Either real number for only one theta, or "binary" {0.1,0.3} or "con" for continuous values (0.1,0.3))
# var = Size of the variance (Noise-level)

#Required Packages
if(!require("clusterGeneration")) install.packages("clusterGeneration"); library("clusterGeneration")
if(!require("mvtnorm")) install.packages("mvtnorm"); library("mvtnorm")


datagen <- function(p=20, n=1000,bA=0.5,setting,rho=0) {
  

  
  mean_x = 0 
  sig_x = 1
  rho = 0
  # sample size
  pC = pI = pP = 2
  pS = p - (pC+pI+pP)
  var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
  # Set strength of relationship between covariates and outcome
  beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )
  # Set strength of relationship between covariates and treatment
  alpha_v = c( 1, 1,   0,   0, 1, 1,  rep(0,p-6))
  
  # Settings with more covariates 
  if(setting==8 |setting==9 | setting==10){ 
    pC = 6
    pI = pP = 2
    pS = p - (pC+pI+pP)
    var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
    # Set strength of relationship between covariates and outcome
    beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,0.6,0.6, rep(0,p-8) )
    # Set strength of relationship between covariates and treatment
    alpha_v = c( 1, 1,   1,   1, 1, 1,0,0,1,1,  rep(0,p-10) )}
  
  
  names(beta_v) = names(alpha_v) = var.list
  ### set true average treatment effect
  
  
  
  
  ### simulate data
  Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list)) 
  diag(Sigma_x) = sig_x^2
  Mean_x = rep(mean_x,length(var.list))
  data = as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
  names(data) = var.list
  
 

  
  ### Options for D (m_(X))

  if(setting==1) {alpha_v = c( 1, 1,   0,   0, 1, 1,  rep(0,p-6))
    gA_x = rowSums(data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))} # Setting 1
  if(setting==2) {alpha_v = c( 0.4, 0.4,   0,   0, 0.4, 0.4,  rep(0,p-6) )
    gA_x = rowSums(data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))} # Setting 2
  if(setting==3) {alpha_v = c( 1, 1,   0,   0, 1, 1,  rep(0,p-6) )
    gA_x = rowSums(data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))} # Setting 3
  if(setting==4) {gA_x = data[,1]*(1-data[,2]) + data[,5]*(1-data[,6])} # Setting 4
  if(setting==5) {gA_x = 0.8*(data[,1]*data[,2]) + 0.8*data[,5]*data[,6]} # Setting 5
  if(setting==6) {alpha_v = c( 1, 1,   0,   0, 1, 1,  rep(0,p-6))
    gA_x = 2*cos(data[,2]) + rowSums(data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list) ,byrow=TRUE))} # Setting 6
  if(setting==7) {gA_x = 2*(data[,1]>0)*(data[,2]>1) + 2*(data[,5]>0)*(data[,6]>1) + data[,1]*data[,6]} # Setting 7
  
  
  
  if(setting==8) {gA_x = 2*data[,2]*(1-data[,6])  + 2*data[,1]*(data[,9]>1) + 0.5*rowSums(data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))} # Setting 8
  if(setting==9) {gA_x = 0.5*(data[,1]^2) + 0.5*data[,2] - data[,3]*data[,4] + 0.5*data[,6] + 0.5*(data[,9]^2) + 0.5*data[,10]} # Setting 9
  if(setting==10){gA_x = -exp(data[,1]) + 0.4*data[,2] + exp(data[,3]) + 0.4*data[,4] + 0.5*data[,5]^2 + data[,6]*data[,9] + 0.4*data[,10]} # Setting 10
  
  # settings for gY_xA
  if(setting==1 | setting==2){gY_xA = rowSums(data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list) ,byrow=TRUE))}
  if(setting==3 | setting==4 | setting==5 | setting==6 | setting==7){ 
    gY_xA =  0.8*(data[,1]*data[,2]) + 0.8*data[,3]*data[,4] }
  if(setting==8 |setting==9 | setting==10){
    gY_xA =  0.8*(data[,1]*data[,2]) + 0.8*data[,3]*data[,4] + 0.8*(data[,5]*data[,6]) + 0.8*data[,7]*data[,8]}

  
  pA = pnorm(gA_x)
  data$A = as.numeric(rbinom(n,1,pA)) # simulate A 
  
  
  data$Y = gY_xA + rnorm(n=n,sd=sig_x)
  data$Y = data$Y + data$A*bA
  
  out <- list()
  out[[1]] <- data
  out[[2]] <- var.list

  return(out)
}


### Example
dataset <- datagen(setting=1, n=500)

summary(dataset)




