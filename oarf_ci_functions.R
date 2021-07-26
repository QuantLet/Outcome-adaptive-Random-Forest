# Bootstrapping
all_RF <- function(data){
  
  var.list <- colnames(data[,!(names(data) %in% c("Y","A"))])
  p <- length(var.list)
  
  ranger_full <- matrix(NA,2,1)
  ranger_reg <- matrix(NA,2,1)
  ranger_reg_only <- matrix(NA,2,1)
  

  
  train.idx <- sample(nrow(data), 0.5 * nrow(data))
  
  for(i in 1:2){ 
    
    if(i==1){
      data_train <- data[train.idx, ]
      data_test <- data[-train.idx, ]
    }
    if(i==2){
      data_train <- data[-train.idx, ]
      data_test <- data[train.idx, ]
    }
    
    # ranger full
    p_ranger_full <- ranger(y=as.factor(data_train$A),x=data_train[,var.list],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10)
    p_hat_ranger_full <- predict(p_ranger_full,data=data_test[,var.list])$predictions[,2]
    
    
    #########################
    # ranger outcome model
    y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
    y_ranger_out <- ranger(y.form,data=data,importance="impurity_corrected",num.trees = 1000, always.split.variables = "A")
    # standardize var.imp
    var.imp_0 <- ifelse(y_ranger_out$variable.importance<0,0,y_ranger_out$variable.importance)
    var.imp <- abs(round(var.imp_0/max(var.imp_0[-1]),3))
    var.init <- ifelse(y_ranger_out$variable.importance >= mean(y_ranger_out$variable.importance[-1]),1,0)
    #select_outcome[,b] <- var.init[-1]
    # ranger guided and regularized
    #p_ranger_reg <- ranger(y=as.factor(data$A),x=data[,var.list],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10,
    #                      regularization.factor = (var.imp[-1]),regularization.usedepth=FALSE, always.split.variables = c(var.list[var.init[-1]==1]))
    
    # RRF
    var.init_count <- which(var.init[-1]==1)
    p_rrf <- RRF(y=as.factor(data$A),x=data[,var.list],flagReg=1,feaIni = var.init_count,coefReg = var.imp[-1],ntree = 1000 )
    
    # final ranger model with selected features
    #fea_sel <- names(which(p_ranger_reg$variable.importance>0))
    fea_sel <- var.list[p_rrf$feaSet]
    #print(fea_sel)
    if(length(fea_sel)<2){
      sub_train <- as.data.frame(data_train[,fea_sel])
      colnames(sub_train) <- fea_sel
      sub_test <- as.data.frame(data_test[,fea_sel])
      colnames(sub_test) <- fea_sel
      p_ranger_select <- ranger(y=as.factor(data_train$A),x=sub_train,importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10)
      p_hat_ranger_reg <- predict(p_ranger_select,data=sub_test)$predictions[,2]
    }
    if(length(fea_sel)>1){
      p_ranger_select <- ranger(y=as.factor(data_train$A),x=data_train[,fea_sel],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10)
      p_hat_ranger_reg <- predict(p_ranger_select,data=data_test[,fea_sel])$predictions[,2]
    }
    ##########################
    
    ## Regularized RF without initial feature space
    p_ranger_reg_only <- ranger(y=as.factor(data_train$A),x=data_train[,var.list],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10,
                                regularization.factor = (var.imp[-1]),regularization.usedepth=FALSE)
    p_hat_ranger_reg_only <- predict(p_ranger_reg_only,data=data_test[,var.list])$predictions[,2]
    
    
    
    # Overlap bounding
    
    p_hat_ranger_full <-ifelse(p_hat_ranger_full<0.025, 0.025, ifelse(p_hat_ranger_full>.975,.975, p_hat_ranger_full)) # Overlap bounding
    p_hat_ranger_reg <-ifelse(p_hat_ranger_reg<0.025, 0.025, ifelse(p_hat_ranger_reg>.975,.975, p_hat_ranger_reg)) # Overlap bounding
    p_hat_ranger_reg_onl <-ifelse(p_hat_ranger_reg_only<0.025, 0.025, ifelse(p_hat_ranger_reg_only>.975,.975, p_hat_ranger_reg_only)) # Overlap bounding
    
    ranger_full[i,1] <- ATE_est(data_test$Y,p_hat_ranger_full,data_test$A)
    ranger_reg[i,1] <- ATE_est(data_test$Y,p_hat_ranger_reg,data_test$A)
    ranger_reg_only[i,1] <- ATE_est(data_test$Y,p_hat_ranger_reg_only,data_test$A)
    
    
    select_var_full <- p_ranger_full$variable.importance/max(p_ranger_full$variable.importance)>0.05
    #select_var_reg[,b] <- rep(1:p) %in% all_var 
    select_var_reg <- rep(1:p) %in% p_rrf$feaSet
    select_var_reg_only <-  p_ranger_reg_only$variable.importance/max(p_ranger_reg_only$variable.importance)>0.05
    #}
    
  }
  

  ite_ranger_full <- mean(ranger_full)
  ite_ranger_reg <- mean(ranger_reg)
  ite_ranger_reg_only <- mean(ranger_reg_only)
  
  rf_est <- cbind(ite_ranger_full,ite_ranger_reg,ite_ranger_reg_only)
  colnames(rf_est) <- c("est_RF_full","est_OARF","est_RRF" )
  
 
  res <- list(rf_est,cbind(select_var_full,select_var_reg, select_var_reg_only))
  
  return(res)
  
}





# Bootstrapping
all_RF_boot <- function(data,verbose){
  
  var.list <- colnames(data[,!(names(data) %in% c("Y","A"))])
  p <- length(var.list)
  
  ite_ranger_full <- matrix(NA,nboot,1)
  ite_ranger_reg <- matrix(NA,nboot,1)
  ite_ranger_reg_only <- matrix(NA,nboot,1)
  
  
  select_var_full <- matrix(NA,p,nboot)
  select_var_reg <- matrix(NA,p,nboot)
  select_var_reg_only <- matrix(NA,p,nboot)
  
  for(b in 1:nboot){
    
    set.seed(123+b)
    data_boot <- createbootstrappedData(data)
    
    
    ranger_full <- matrix(NA,2,1)
    ranger_reg <- matrix(NA,2,1)
    ranger_reg_only <- matrix(NA,2,1)
    
    
    train.idx <- sample(nrow(data_boot), 0.5 * nrow(data_boot))
    
    for(i in 1:2){ 
      
      if(i==1){
        data_train <- data_boot[train.idx, ]
        data_test <- data_boot[-train.idx, ]
      }
      if(i==2){
        data_train <- data_boot[-train.idx, ]
        data_test <- data_boot[train.idx, ]
      }
      
      # ranger full
      p_ranger_full <- ranger(y=as.factor(data_train$A),x=data_train[,var.list],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10)
      p_hat_ranger_full <- predict(p_ranger_full,data=data_test[,var.list])$predictions[,2]
      
      
      #########################
      # ranger outcome model
      y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
      y_ranger_out <- ranger(y.form,data=data,importance="impurity_corrected",num.trees = 1000, always.split.variables = "A")
      # standardize var.imp
      var.imp_0 <- ifelse(y_ranger_out$variable.importance<0,0,y_ranger_out$variable.importance)
      var.imp <- abs(round(var.imp_0/max(var.imp_0[-1]),3))
      var.init <- ifelse(y_ranger_out$variable.importance >= mean(y_ranger_out$variable.importance[-1]),1,0)
      #var.init
      #select_outcome[,b] <- var.init[-1]
      # ranger guided and regularized
      #p_ranger_reg <- ranger(y=as.factor(data$A),x=data[,var.list],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10,
      #                      regularization.factor = (var.imp[-1]),regularization.usedepth=FALSE, always.split.variables = c(var.list[var.init[-1]==1]))
      
      # RRF
      var.init_count <- which(var.init[-1]==1)
      p_rrf <- RRF(y=as.factor(data$A),x=data[,var.list],flagReg=1,feaIni = var.init_count,coefReg = var.imp[-1],ntree = 1000 )
      
      # final ranger model with selected features
      #fea_sel <- names(which(p_ranger_reg$variable.importance>0))
      fea_sel <- var.list[p_rrf$feaSet]
      #print(fea_sel)
      if(length(fea_sel)<2){
        sub_train <- as.data.frame(data_train[,fea_sel])
        colnames(sub_train) <- fea_sel
        sub_test <- as.data.frame(data_test[,fea_sel])
        colnames(sub_test) <- fea_sel
        p_ranger_select <- ranger(y=as.factor(data_train$A),x=sub_train,importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10)
        p_hat_ranger_reg <- predict(p_ranger_select,data=sub_test)$predictions[,2]
      }
      if(length(fea_sel)>1){
        p_ranger_select <- ranger(y=as.factor(data_train$A),x=data_train[,fea_sel],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10)
        p_hat_ranger_reg <- predict(p_ranger_select,data=data_test[,fea_sel])$predictions[,2]
      }
      ##########################
      
      ## Regularized RF without initial feature space
      p_ranger_reg_only <- ranger(y=as.factor(data_train$A),x=data_train[,var.list],importance="impurity",probability = TRUE,num.trees = 500,min.node.size = 10,
                                  regularization.factor = (var.imp[-1]),regularization.usedepth=FALSE)
      p_hat_ranger_reg_only <- predict(p_ranger_reg_only,data=data_test[,var.list])$predictions[,2]
      
      
      
      # Overlap bounding
      
      p_hat_ranger_full <-ifelse(p_hat_ranger_full<0.025, 0.025, ifelse(p_hat_ranger_full>.975,.975, p_hat_ranger_full)) # Overlap bounding
      p_hat_ranger_reg <-ifelse(p_hat_ranger_reg<0.025, 0.025, ifelse(p_hat_ranger_reg>.975,.975, p_hat_ranger_reg)) # Overlap bounding
      p_hat_ranger_reg_onl <-ifelse(p_hat_ranger_reg_only<0.025, 0.025, ifelse(p_hat_ranger_reg_only>.975,.975, p_hat_ranger_reg_only)) # Overlap bounding
      
      ranger_full[i,1] <- ATE_est(data_test$Y,p_hat_ranger_full,data_test$A)
      ranger_reg[i,1] <- ATE_est(data_test$Y,p_hat_ranger_reg,data_test$A)
      ranger_reg_only[i,1] <- ATE_est(data_test$Y,p_hat_ranger_reg_only,data_test$A)
      
      # Safe selected variables
      #rep(1:p) %in% p_rrf$feaSet
      
      #if(i==1){sel_var <- p_ranger_reg$variable.importance/max(p_ranger_reg$variable.importance)>0.05}
      #if(i==2){all_var <- unique(c(which(sel_var),which(p_ranger_reg$variable.importance/max(p_ranger_reg$variable.importance)>0.05)))
      
      select_var_full[,b] <- p_ranger_full$variable.importance/max(p_ranger_full$variable.importance)>0.05
      #select_var_reg[,b] <- rep(1:p) %in% all_var 
      select_var_reg[,b] <- rep(1:p) %in% p_rrf$feaSet
      select_var_reg_only[,b] <-  p_ranger_reg_only$variable.importance/max(p_ranger_reg_only$variable.importance)>0.05
      #}
      
    }
    
    ite_ranger_full[b,1] <- mean(ranger_full)
    ite_ranger_reg[b,1] <- mean(ranger_reg)
    ite_ranger_reg_only[b,1] <- mean(ranger_reg_only)
    
    
    if(verbose==T)
    {cat("This is iteration", b, "out of", nboot, "\n")}
    
  }
  rf_boot <- cbind(ite_ranger_full,ite_ranger_reg,ite_ranger_reg_only)
  #boot_sds <- apply(rf_boot, 2, function(x){ sd(x) })
  #boot_ests <- colMeans(rf_boot, na.rm = TRUE) # smoothed bootstrap
  boot_q_lower <- apply(rf_boot, 2, function(x){ quantile(x,0.05) })
  boot_q_upper <- apply(rf_boot, 2, function(x){ quantile(x,0.95) })
  ci_all <- rbind(boot_q_lower,boot_q_upper)
  
  # Covariate selection
  sv_nor=apply(select_var_full,1,mean)
  sv_reg=apply(select_var_reg,1,mean)
  sv_reg_only <- apply(select_var_reg_only,1,mean)
  
  res <- list(ci_all,cbind(sv_nor,sv_reg,sv_reg_only))
  
  return(res)
  
}


