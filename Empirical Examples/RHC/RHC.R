library(ranger)
library(scales)
library(MASS)
library(ggplot2)
library(cowplot)
library(glmnet)
library(RRF)



var.list <- colnames(rhc.dt[,!(names(rhc.dt) %in% c("Y","A"))])

# Point estimates

# RF, OARF and RRF
B <- 500
RF_est_B <- matrix(NA,B,3)
for(b in 1:B){
RF_est <- all_RF(data=rhc.dt)
RF_est_B[b,] <- RF_est

}
RF_est_med <- apply(RF_est_B,2,median)

# OAL 
data_m <- rhc.dt %>% 
  mutate_if(is.logical, ~as.numeric(.))
data_m[,var.list] <-  rapply(data_m[,var.list],scale,c("numeric","integer"),how="replace")
OAL_est <- shortreed_est(data=data_m)

# IPTW
IPTW_est <- iptw(data=data_m)



###############
# Bootstrap iterations
nboot <- B

# RF, OARF and RRF
RF_ci <- all_RF_boot(data=rhc.dt)

# OAL 
data_m <- rhc.dt %>% 
  mutate_if(is.logical, ~as.numeric(.))
data_m[,var.list] <-  rapply(data_m[,var.list],scale,c("numeric","integer"),how="replace")
OAL_ci <- shortreed_boot(data=data_m,nboot=nboot)

# IPTW
IPTW_ci <- iptw_boot_ci(data=data_m,nboot=nboot)

# ATE Estimates 
c(IPTW_est,IPTW_ci)
c(OAL_est,OAL_ci$iptw)
c(RF_est_med[1],RF_ci[[1]][,1]) # RF full
c(RF_est_med[3],RF_ci[[1]][,3]) # RF RRF
c(RF_est_med[2],RF_ci[[1]][,2]) # OARF

# CI width
IPTW_ci[2] - IPTW_ci[1]
OAL_ci$iptw[2] - OAL_ci$iptw[1]
RF_ci[[1]][,1][2] - RF_ci[[1]][,1][1]
RF_ci[[1]][,3][2] - RF_ci[[1]][,3][1]
RF_ci[[1]][,2][2] - RF_ci[[1]][,2][1]


# Covariate selection
sv_OAL <- OAL_ci$var_sel
sv_RF_full <- as.data.frame(RF_ci[[2]][,1])
sv_OARF <- RF_ci[[2]][,2]
sv_RRF <- RF_ci[[2]][,3]



  


sv <- as.data.frame(c(sv_RF_full,sv_OARF,sv_RRF,sv_OAL))
colnames(sv) <- c("value")
sv$Method <- factor(rep(c("RF full","OARF","RRF","OAL"),each=length(var.list)),levels=c("OAL","RF full","OARF","RRF"))
sv$Var <- rep(c(1:length(var.list)))

#cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(sv,aes(y=value,x=Var,color=Method))+
  geom_line(size=1) +
  theme_cowplot() +
  labs(y="Proportion of times covariates selected ",x="Covariates") +
  scale_color_manual(values=c("#000000", "#E69F00","#56B4E9","#009E73"))





# Check covariates

all_var_table <- cbind(var.list,sv_RF_full,sv_RRF,sv_OARF,sv_OAL)
all_var_table[,-1] <- round(all_var_table[,-1],3)
xtable(all_var_table,digits=c(0,0,1,1,1,1))
