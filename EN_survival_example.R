###################################
## simulation for empirical null ##
## survival outcome, accounting for total variation ##

rm(list=ls())

library(MASS)
library(fmsb)
library(data.table)
library(survival)
library(mvtnorm)

source("EN_lib.R")

set.seed(101010)


### parameter setup ###

nloop <- 500 # number of simulation replications
num_f <- 2000 # number of facilities
size_f <- as.integer(sample(10:200, size=num_f, replace=T)) # number of patients in each provider
size_f <- size_f[order(size_f)] # important to order the providers by their size
N <- sum(size_f)
p <- 2
g_seq <- c(5,20,40,60) # try different numbers of groups G
niter <- 10
tol <- 1.0e-5
sig_level <- 0.05
v <- qnorm(sig_level, lower.tail = F)
p.grid <- seq(0.8,1,0.001)

tert <- quantile(size_f, probs=c(1/3, 2/3))
size_idx <- 1*as.numeric(size_f <= tert[1]) + 2*as.numeric((size_f > tert[1])&(size_f <= tert[2])) + 3*as.numeric(size_f > tert[2])
sum(size_idx==1); sum(size_idx==2); sum(size_idx==3)

## results ##

for(g in g_seq) 
  for(j in c("int", "noint"))
    for(k in c("mean", "median"))
    {
      eval(parse(text=paste("coef_G", g, "_", j, "_", k, " <- NULL", sep="")))
      eval(parse(text=paste("flag_G", g, "_", j, "_", k, " <- NULL", sep="")))
    }
est_slope <- NULL # ratio of sigmas
flag_fe <- NULL


## computation ##

loop <- 1
repeat{# repeat1
  gamma_subject=rnorm(num_f, mean=0, sd=0.2)
  gamma_subject=rep(gamma_subject, size_f)
  F_pre=1:num_f
  facility=rep(F_pre, size_f)
  Sigma_z1=diag(p)
  z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  U=runif(N, 0,1)
  pre_time=-log(U)/(0.1*exp(gamma_subject+z[,1]-z[,2]))
  
  pre_censoring=runif(N,10,30)
  tcens=(pre_censoring<pre_time) # censoring indicator
  delta=1-tcens; mean(tcens)
  time=pre_time*(delta==1)+pre_censoring*(delta==0)
  
  ### order data
  delta = delta[order(time)]
  facility=facility[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  
  cox_stratify=coxph(Surv(time,delta)~z+strata(facility))
  
  beta=cox_stratify$coefficients
  
  S0=rev(cumsum(rev(exp(z%*%beta))))
  Lambda=cumsum(delta/S0)
  
  # martigale calculation
  cox_unstratify=coxph(Surv(time,delta)~offset(z%*%beta))
  residuals=residuals(cox_unstratify, type="martingale")
  E=exp(z%*%beta)*cumsum(delta/S0)
  
  facility_sum=function(x){
    as.vector(sapply(split(x,factor(facility)), sum))
  }
  
  O_F=as.vector(sapply(split(delta,factor(facility)), sum))
  E_F=as.vector(sapply(split(E,factor(facility)), sum))
  SMR=O_F/E_F
  
  p_oneside <- (ppois(O_F-1, lambda=E_F, lower.tail=F) +
                  ppois(O_F, lambda=E_F, lower.tail=F))/2
  p_twoside <- 2*ifelse(p_oneside>0.5, 1-p_oneside, p_oneside)
  
  flag_fe <- rbind(flag_fe, as.numeric(p_oneside < sig_level))
  
  # Z-score
  
  zscore <- ifelse(SMR > 1, 1, -1)*ifelse(p_twoside != 0, 
                                          qnorm(p_twoside/2, lower.tail=F), 4.2649)  
  
  
  my_data <- data.frame(zscore=zscore, size=size_f)
  for(g in g_seq) {
    eval(parse(text=paste("estVar_G", g, " <- estVar_split(df=my_data, number=", g, ", method='locfdr.optim', p.grid=p.grid)", sep="")))
  }
  
  for(g in g_seq) 
    for(j in c("int", "noint"))
      for(k in c("mean", "median"))
      {
        fixed <- "FALSE"
        if(j=="noint") fixed <- "TRUE"
        eval(parse(text=paste("res_", j, " <- irw.group.var(estVar_G", g, ", niter, tol, ",
                              fixed, ", '", k, "')", sep="")))
        eval(parse(text=paste("coef_G", g, "_", j, "_", k, " <- rbind(coef_G", g, "_", j, "_",
                              k, ", res_", j, "$coef)",  sep="")))
      }
  
  for(g in g_seq)
  {
    k <- "median"
    j <- "int"
    fixed <- "FALSE"
    eval(parse(text=paste("res_", j, " <- irw.group.var(estVar_G", g, ", niter, tol, ",
                          fixed, ", '", k, "')", sep="")))
    eval(parse(text=paste("weights1 <- cbind(rep(1,", g, "), estVar_G", g, "$", k, "_size)%*%res_int$coef", sep="")))
    weights1 <- as.vector(1/weights1)
    eval(parse(text=paste("max_median_size <- max(estVar_G", g, "$median_size)", sep="")))
    eval(parse(text=paste("min_median_size <- min(estVar_G", g, "$median_size)", sep="")))
    size_range2 <- size_f[(size_f >= min_median_size)&(size_f <= max_median_size)]
    size_range1 <- size_f[(size_f < min_median_size)]
    size_range3 <- size_f[(size_f > max_median_size)]
    eval(parse(text=paste("mspline <- smooth.spline(x=estVar_G", g, "$", k, "_size, y=estVar_G", g, "$mean_est, w=weights1)", sep="")))
    eval(parse(text=paste("pred_spline <- predict(mspline, data.frame(median_size=size_range2))", sep=""))) 
    pred_mean <- as.vector(pred_spline$y[,1])
    pred_mean <- c(rep(pred_mean[1], length(size_range1)), pred_mean, rep(pred_mean[length(pred_mean)], length(size_range3)))
    eval(parse(text=paste("pred_var <- as.vector(cbind(rep(1,", length(size_f), "), size_f)%*%res_int$coef)", sep="")))
    eval(parse(text=paste("flag_G", g, "_", j, "_", k, " <- rbind(flag_G", g, "_", j, "_", k, 
                          ", zscore > (pred_mean + v*sqrt(pred_var)))", sep="")))
  }
  
  loop <- loop + 1
  if(loop > nloop) break
  
}# end repeat1


### plotting results ###
# take G=20 for example

# plot mean function of Z-scores #
g <- 20; k <- "median"
plot(x=range(size_f), y=c(-4,4), col="white", xlab=paste("Median", " size", sep=""),
     ylab="Est. Mean")
eval(parse(text=paste("points(x=estVar_G", g, "$", k, "_size, y=estVar_G", g, "$mean_est, pch=20)", sep="")))
eval(parse(text=paste("weights1 <- cbind(rep(1,", g, "), estVar_G", g, "$", k, "_size)%*%res_int$coef", sep="")))
weights1 <- as.vector(1/weights1)
eval(parse(text=paste("mspline <- smooth.spline(x=estVar_G", g, "$", k, "_size, y=estVar_G", g, "$mean_est, w=weights1)", sep="")))
eval(parse(text=paste("pred_spline <- predict(mspline, data.frame(median_size=estVar_G", g, "$median_size))", sep="")))
eval(parse(text=paste("lines(x=estVar_G", g, "$median_size, y=as.vector(pred_spline$y[,1]), col='orange', lwd=1.5)", sep="")))
eval(parse(text=paste("lines(x=c(size_f[size_f<=estVar_G", g,
                      "$median_size[1]], estVar_G", g, "$median_size[1]), y=rep(as.numeric(head(pred_spline$y, n=1)), sum(size_f<=estVar_G", g,
                      "$median_size[1])+1), col='orange', lwd=1.5)", sep="")))
eval(parse(text=paste("lines(x=c(estVar_G", g, "$median_size[", g, "], size_f[size_f>=estVar_G", g,
                      "$median_size[", g, "]]), y=rep(as.numeric(tail(pred_spline$y, n=1)), sum(size_f>=estVar_G", g,
                      "$median_size[", g, "])+1), col='orange', lwd=1.5)", sep="")))


# plot variance function of Z-scores #
plot(x=range(size_f), y=c(0,10), col="white", xlab=paste("Median", " size", sep=""),
     ylab="Est. Variance")
eval(parse(text=paste("points(x=estVar_G", g, "$", k, "_size, y=estVar_G", g, "$var_est, pch=20)", sep="")))
for(j in c("int")) {
  if(j=="int") {
    eval(parse(text=paste("abline(a=", "coef_G", g, "_", j, "_", k,
                          "[", 1, ",1], b=", "coef_G", g, "_", j, "_", k,
                          "[", 1, ",2], col='steelblue3', xlim=range(size_f), lwd=1.5)", sep="")))
  } else {
    eval(parse(text=paste("abline(a=1", ", b=", "coef_G", g, "_", j, "_", k,
                          "[", 1, ",1], col='orange', xlim=range(size_f), lwd=1.5)", sep="")))
  }
}


# plot flagging rates (figures will look different due to data generating mechanism) #
lambda_seq <- c(0.5, 0.75, 1)
mean_by_idx <- function(x) {
  sapply(split(x, factor(size_idx)), mean)
}
strata_flag_fe <- apply(flag_fe, 2, mean)
g <- 20
for(i in 1:length(lambda_seq)) {
  lam <- lambda_seq[i]
  eval(parse(text=paste("strata_flag_G", g, "_lam", 100*lam, " <- apply(flag_G", g, "_int_median_lam", 100*lam, ", 2, mean)", sep="")))
}
stacked_data <- data.frame(mean=c(strata_flag_fe, strata_flag_G20_lam50, strata_flag_G20_lam75,
                                  strata_flag_G20_lam100))
stacked_data$method <- c(rep("FE", num_f), rep("EN_lam3", num_f), rep("EN_lam2", num_f),
                         rep("EN_lam1", num_f))  # 3: lam50; 2: lam75; 1: lam100
stacked_data$size <- rep(c(rep("1Small", sum(size_idx==1)), rep("2Medium", sum(size_idx==2)),
                           rep("3Large", sum(size_idx==3))), 4)
stacked_data <- melt(stacked_data, id=c("method", "size"))
stacked_data <- stacked_data[,-3]

boxplot(value ~ method + size, data=stacked_data, at = c(1:4, 6:9, 11:14),
        col=c("orange", "blue", "white", "gray"), ylim=c(0,0.35), lwd=1,
        ylab="Probability of Signal", xaxt="n")
axis(side=1, at=c(2.5, 7.5, 13.5), labels=c("Small", "Medium", "Large"), line=0.5, lwd=0)
legend("topleft", fill=c("orange", "blue", "white", "gray"),
       pch=15, legend=c(expression(paste("EN (", lambda, "=1)", sep="")),
                        expression(paste("EN (", lambda, "=0.75)", sep="")),
                        expression(paste("EN (", lambda, "=0.5)", sep="")),
                        expression(paste("FE / EN (", lambda, "=0)", sep=""))), cex=0.9)