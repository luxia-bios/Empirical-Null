### This R file contains functions that will be called for the empirical null project ###

## split data frame 
df_split <- function(df, number){
  sizedf <- length(df[,1])
  bound <- floor(sizedf/number)
  fold_idx <- cut(seq(1,nrow(df)), breaks=number, labels=F)
  list <- list() 
  for (i in 1:number){
    list[[i]] <- df[fold_idx==i,]
  }
  return(list)
}

## use R function, optim, for direct optimization using Nelder-Mead ##
modifEN.optim1 <- function(z, xlim, p.grid) {
  N = length(z)
  if (missing(xlim)) {
    if (N>500000) b = 1
    else b = 4.3 * exp(-0.26*log(N,10)) 
    xlim = c(median(z),b*diff(quantile(z)[c(2,4)])/(2*qnorm(.75)))
  }
  aorig = xlim[1]-xlim[2]
  borig = xlim[1]+xlim[2]
  z0 = z[which(z>=aorig & z<=borig)]
  N0 = length(z0)
  N1 = N - N0
  niter <- length(p.grid)
  eval.p.grid <- array(NA, niter)
  for(i in 1:niter){
    p0 <- p.grid[i]
    negloglik <- function(arg) {
      mu0 <- arg[1]
      sig0 <- arg[2]
      Q <- pnorm((borig - mu0)/sig0) - pnorm((aorig - mu0)/sig0)
      loglik <- N0*log(p0) - N0*log(sig0) - sum((z0-mu0)**2)/(2*sig0**2) + N1*log(1 - p0*Q)
      return((-1)*loglik)
    }
    res.optim <- optim(par=c(mean(z0), sd(z0)), fn=negloglik, hessian = FALSE, method="BFGS") #"Nelder-Mead"
    eval.p.grid[i] <- res.optim$value
  } # end for niter
  p0 <- p.grid[which.min(eval.p.grid)]
  negloglik <- function(arg) {
    mu0 <- arg[1]
    sig0 <- arg[2]
    Q <- pnorm((borig - mu0)/sig0) - pnorm((aorig - mu0)/sig0)
    loglik <- N0*log(p0) - N0*log(sig0) - sum((z0-mu0)**2)/(2*sig0**2) + N1*log(1 - p0*Q)
    return((-1)*loglik)
  }
  res.optim <- optim(par=c(mean(z0), sd(z0)), fn=negloglik, hessian = FALSE, method="BFGS") #"Nelder-Mead"
  return(list(est=c(res.optim$par, p0), conv=res.optim$convergence, message=res.optim$message))
}

# after data frame splitting, compute group-wise mean and variance
estVar_split <- function(df, number, method, scale.est="MAD", psi=psi.huber, p.grid=NULL) ## methodd=c("lm", "m", "mm")
{
  list_split <- df_split(df, number)
  mean_est <- rep(NA, length(list_split))
  var_est <- rep(NA, length(list_split))
  median_size <- rep(NA, length(list_split))
  mean_size <- rep(NA, length(list_split))
  group_size <- rep(NA, length(list_split))
  for(i in 1:length(list_split))
  {
    size_tmp <- list_split[[i]]$size
    zscore_tmp <- list_split[[i]]$zscore
    median_size[i] <- median(list_split[[i]]$size)
    mean_size[i] <- mean(list_split[[i]]$size)
    group_size[i] <- length(size_tmp)
    if(method=="lm") # linear regression
    {
      model_tmp <- lm(zscore ~ 1, data=list_split[[i]]) 
      mean_est[i] <- as.numeric(coef(model_tmp))
      var_est[i] <- as.numeric(summary(model_tmp)$sigma)**2
    } else if(method=="m") # robust M-estimation
    {
      model_tmp <- rlm(zscore ~ 1, data=list_split[[i]], method="M", maxit=200, scale.est=scale.est, psi=psi)
      mean_est[i] <- as.numeric(coef(model_tmp))
      var_est[i] <- as.numeric(model_tmp$s)**2
    } else if(method=="mm") # robust MM-estimation
    {
      model_tmp <- rlm(zscore ~ 1, data=list_split[[i]], method="MM", maxit=200, scale.est=scale.est, psi=psi)
      mean_est[i] <- as.numeric(coef(model_tmp))
      var_est[i] <- as.numeric(model_tmp$s)**2
    } else if(method=="locfdr.optim") { # local MLE fitting, modified algorithm
      if(is.null(p.grid)) {
        stop("p.grid needs to be specified!")
      }
      y <- list_split[[i]]$zscore
      model_tmp <- rlm(zscore ~ 1, data=list_split[[i]], method="MM", maxit=200, scale.est=scale.est, psi=psi)
      aaa <- 1.65
      res_modifEN_optim <- modifEN.optim1(z=y, p.grid=p.grid)
      var_est[i] <- (res_modifEN_optim$est[2])**2
      mean_est[i] <- res_modifEN_optim$est[1]
    } else
    {
      stop("Invalid input for method!")
    }
  } 
  return(data.frame(mean_est=mean_est, var_est=var_est, median_size=median_size,
                    mean_size=mean_size, group_size=group_size))
}


# iteratively reweighted estimation of variance as a linear function
irw.group.var <- function(estVar_obj, niter, tol=1.0e-5, fixedInt, xlabel) {
  est_coef <- NULL
  iter <- 1
  diff <- Inf
  conv <- FALSE
  yvar <- estVar_obj$var_est
  if(fixedInt==TRUE) yvar <- (estVar_obj$var_est - 1)
  xvar <- NULL
  if(xlabel=="mean") {
    xvar <- estVar_obj$mean_size
  } else if(xlabel=="median") {
    xvar <- estVar_obj$median_size
  } else {
    stop("Invalid input for x-axis label.")
  }
  lm_old <- NULL
  tmp_data <- data.frame(yvar=yvar, xvar=xvar)
  if(fixedInt==TRUE) {
    lm_old <- lm(yvar ~ xvar - 1, data=tmp_data)
  } else {
    lm_old <- lm(yvar ~ xvar, data=tmp_data)
  }
  repeat { # repeat
    pre_weights <- predict(lm_old, data.frame(xvar=xvar))
    pre_weights <- estVar_obj$group_size/(pre_weights**2)
    if(fixedInt==TRUE) {
      lm_new <- lm(yvar ~ xvar - 1, data=tmp_data, weights=pre_weights)
    } else {
      lm_new <- lm(yvar ~ xvar, data=tmp_data, weights=pre_weights)
    }
    iter <- iter + 1
    diff <- sqrt(sum((coef(lm_old)-coef(lm_new))**2))
    if(iter > niter) {
      break
    } else if(diff<tol) {
      conv <- TRUE
      break
    } else {
      lm_old <- lm_new
    }
  } # end repeat
  
  return(list(coef=coef(lm_new), converge=conv))
}


