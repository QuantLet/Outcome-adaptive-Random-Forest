install.packages(vec.pac)
vec.pac= c("mvtnorm","clusterGeneration",
           "ranger","scales","MASS","glmnet","RRF","lqa","ggplot2","cowplot")
lapply(vec.pac, require, character.only = TRUE) 

sv_all <- list()
S <- 10

for(s in 1:S){


B<- 200
p <- 20

select_var_full <- matrix(NA,p,B)
select_var_reg <- matrix(NA,p,B)
select_var_reg_only <- matrix(NA,p,B)
select_var_OAL <- matrix(NA,p,B)

for( b in 1:B) {
  
  DGP <- datagen(setting=s, n=500,p=p)
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
  RF_out <- all_RF(data)
  
  # Covariate selection
  select_var_full[,b] <- RF_out[[2]][,1]
  #select_var_reg[,b] <- rep(1:p) %in% all_var 
  select_var_reg[,b] <- RF_out[[2]][,2]
  select_var_reg_only[,b] <-  RF_out[[2]][,3]
  #}
  
 
  
  # OAL 
  OAL_results <- OAL(data,var.list)
  select_var_OAL[,b] <- abs(OAL_results[[2]])>0.00001
  
  
  #cat("This is iteration", b, "out of", B, "\n")

}

sv_nor=apply(select_var_full,1,mean)
sv_reg=apply(select_var_reg,1,mean)
sv_reg_only <- apply(select_var_reg_only,1,mean)
sv_OAL=apply(select_var_OAL,1,mean)

sv <- as.data.frame(c(sv_nor,sv_reg,sv_reg_only,sv_OAL))
colnames(sv) <- c("value")
sv$Method <- factor(rep(c("RF full","OARF","RRF","OAL"),each=p),levels=c("OAL","RF full","OARF","RRF"))
sv$Var <- rep(c(1:p))

#cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


sv_all[[s]] <- sv_all
cat("This is setting", s, "out of", S, "\n")

}

ggplot(sv,aes(y=value,x=Var,color=Method))+
  geom_line(size=1) +
  geom_point(aes(shape=Method)) +
  theme_cowplot() +
  labs(y="Proportion of times covariates selected ",x="Covariates") +
  scale_color_manual(values=c("#000000", "#E69F00","#56B4E9","#009E73"))


