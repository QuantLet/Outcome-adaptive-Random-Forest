install.packages("H:/Nonpara adaptive RF/lqa_1.0-3.tar", type="source", repos=NULL)
library(lqa)





ATE_est = function(fY,fp,fA){
  
  
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}


OAL <- function(data,var.list){
  
  create_weights = function(fp,fA,fw){
    fw = (fp)^(-1)
    fw[fA==0] = (1 - fp[fA==0])^(-1)
    return(fw)
  }
  
  
  wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
    trt = untrt = diff_vec = rep(NA,length(beta)) 
    names(trt) = names(untrt) = names(diff_vec) = varlist
    for(jj in 1:length(varlist)){ 
      this.var = paste("w",varlist[jj],sep="") 
      DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
      trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
      untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
      diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
    } 
    wdiff_vec = diff_vec * abs(beta) 
    wAMD = c( sum(wdiff_vec))
    ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
    return(ret) 
  }
  
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  lm.Y = lm(y.form,data=data)
  betaXY = coef(lm.Y)[var.list] 
  
  
  
  lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  
  
  
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  
  
  
  
  
  p_hat_OAL <- as.data.frame(matrix(NA,nrow(data),length(lambda_vec)))
  colnames(p_hat_OAL) <- names(lambda_vec)
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    
    
    
    ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = adaptive.lasso(lambda=nrow(data)^(il),al.weights = abs(betaXY)^(-ig) )
    ### run outcome-adaptive lasso model with appropriate penalty
    logit_oal = lqa.formula( w.full.form, data=data, penalty=oal_pen, family=binomial(logit) )
    
    data[,paste("f.pA",lil,sep="")] = predict.lqa(logit_oal)$mu.new
    
    # save propensity score coefficients
    coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
    # create inverse probability of treatment weights
    data[,paste("w",lil,sep="")] = create_weights(fp=data[,paste("f.pA",lil,sep="")],fA=data$A)
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=data,varlist=var.list,trt.var="A",
                                  wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
    # save ATE estimate for this lambda value
    ATE[lil] = ATE_est(fY=data$Y,fp=data[,paste("f.pA",lil,sep="")],fA=data$A)
  } # close loop through lambda values
  
  # print out wAMD for all the lambda values tried
  wAMD_vec
  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)
  # print out ATE corresponding to smallest wAMD value
  ATE[tt]
  # print out the coefficients for the propensity score that corresponds with smalles wAMD value 
  coeff_XA[,tt]
  
  res <- list(ATE[tt],coeff_XA[,tt])
  
  return(res)
  
}



createbootstrappedData <- function(df_boot) {
  
  smpl_0 <- sample((1:nrow(df_boot))[df_boot$A == 0],
                   replace = TRUE,
                   size = sum(1 - df_boot$A))
  smpl_1 <- sample((1:nrow(df_boot))[df_boot$A == 1],
                   replace = TRUE,
                   size = sum(df_boot$A))
  smpl <- sample(c(smpl_0, smpl_1))
  
  return(df_boot[smpl,])
}


shortreed_est <- function(data, family = binomial()){
  OAL_res <- OAL(data,var.list)
  return(OAL_res[[1]])
}


#' Function to do one bootstrap iteration of Shortreed-based estimators
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param family Family for outcome regression for glm
one_shortreed_boot <- function(data, family = binomial()){
  n <- nrow(data)
  idx <- sample(1:n, replace = TRUE)
  Yij_vec <- sapply(1:n, function(x,idx){
    sum(idx == x)
  }, idx = idx)
  data_boot <- data[idx,]
  OAL_res <- OAL(data_boot,var.list)
  t_star <- OAL_res[[1]]
  coef_var <- abs(OAL_res[[2]])>0.00001
  return(list(Yij = Yij_vec, t_star = t_star,var_sel = coef_var))
}


shortreed_boot <- function(data, nboot = 5e2, family = binomial()){
  rslt <- replicate(nboot, one_shortreed_boot(data, family = family))
    all_Yijs <- Reduce(rbind, rslt[1,])
  all_tstars <- apply(Reduce(rbind, rslt[2,]), 2, unlist, use.names= FALSE)
  all_covs <- apply(all_tstars, 2, function(tstar){
    apply(all_Yijs, 2, function(x){
      cov(x, tstar, use = "complete.obs")
    })
  })
  boot_sds <- apply(all_covs, 2, function(x){ sqrt(sum(x^2)) })
  boot_ests <- colMeans(all_tstars, na.rm = TRUE)
  cis <- rbind(boot_ests - 1.96*boot_sds, boot_ests + 1.96*boot_sds)
  # Covariate selection
  var <- Reduce(rbind, rslt[3,])
  var_mean <- apply(var,2,mean,na.rm=T)
  return(list(iptw = cis[,1],var_sel = var_mean ))
}



iptw <- function(data, formula = "A ~ .", family = binomial()){
  ps_fit <- glm(as.formula(formula), data = data[,c("A",var.list)], family = family)
  g1hat <- predict(ps_fit, type = "response")
  est <- sum(data$A * data$Y /(g1hat))/sum(data$A/g1hat) - sum((1-data$A)*data$Y/(1- g1hat))/sum((1-data$A)/(1-g1hat))
  return(est)
}


iptw_one_boot <- function(data, formula = "A ~ .", family = binomial()){
  resamp_idx <- sample(1:nrow(data), replace = TRUE)
  return(iptw(data=data[resamp_idx,],
              formula = formula, family = family))
}


iptw_boot_ci <- function(data, formula = "A ~ .", nboot = 5e2, family = binomial()){
  ate_star <- replicate(nboot, iptw_one_boot(data, formula, family = family))
  return(as.numeric(quantile(ate_star, p = c(0.025, 0.975))))
}

#' Help function to check whether true value in a CI
#' @param truth The true value
#' @param ci A two-length vector CI
truth_in_ci <- function(truth, ci){
  truth > min(ci) & truth < max(ci)
}

