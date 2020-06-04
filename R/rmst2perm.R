#' @name rmst2perm
#' @aliases rmst2perm
#' @title Permutation Test for Comparing Restricted Mean Survival Time
#' @description Performs the permutation test using difference in the restricted mean
#' survival time (RMST) between groups as a summary measure of the survival time distribution.
#' When the sample size is less than 50 per group, it has been shown that there is non-negligible
#' inflation of the type I error rate in the commonly used asymptotic test for the RMST comparison.
#' Generally, permutation tests can be useful in such a situation.
#' However, when we apply the permutation test for the RMST comparison, particularly in small sample
#' situations, there are some cases where the survival function in either group cannot be defined
#' due to censoring in the permutation process. Horiguchi and Uno (2020) <doi:10.1002/sim.8565>
#' have examined six workable solutions to handle this numerical issue.
#' It performs permutation tests with implementation of the six methods outlined in the paper
#' when the numerical issue arises during the permutation process.
#' The result of the asymptotic test is also provided for a reference.
#' @usage rmst2perm(time, status, arm, tau=NULL, mperm=c(1:6), nperm=100000,
#' seed=NULL, asy="greenwood", test="2_side")
#' @param time The follow-up time for right censored data.
#' @param status The event indicator, 1=event, and 0=censor.
#' @param arm The group indicator for comparison with a value of either 0 or 1. Normally, 0=control group, 1=active treatment group. The three vectors (\code{time}, \code{status}, \code{arm}) need to have the same length.
#' @param tau A scaler value to specify the truncation time point for the RMST calculation.
#' It needs to be smaller than the minimum value of the largest observed time in each of the two groups.
#' @param mperm A vector with the numbers from 1 to 6 to specify the method for conducting the permutation test when the last observation time from either group does not reach the specified \code{tau}.
#' It supports: 1=ignoring the inestimable cases (Method 1),
#' 2=extending the survival curve to tau (Method 2),
#' 3=switching the last censored observation to the event observation (Method 3),
#' 4=averaging RMSTs derived from Methods 2 and 3 (Method 4),
#' 5=fitting a Weibull distribution to each inestimable case (Method 5), and
#' 6=utilizing pseudo-observations (Method 6). Please see Horiguchi and Uno (2020) <doi:10.1002/sim.8565> for details.
#' @param nperm The number of iterations for the resampling. It is recommended to specify at least 100,000 (default) or larger.
#' @param seed An integer value, used for random number generation in the resampling procedure. Default is \code{NULL}.
#' @param asy Specify the asymptotic variance estimator for the difference in RMST.
#' \code{asy} supports \code{"greenwood"} for Greenwood plug-in estimator (default) and \code{"aj"} for Aalen-Johansen plug-in estimators. Please see Horiguchi and Uno (2020) <doi:10.1002/sim.8565> for details.
#' @param test Specify \code{"1_side"} for the one-sided test where the alternative hypothesis is that the treatment effect in arm=1 is superior to that in arm=0 with respect to survival.
#' Specify \code{"2_side"} for the two-sided test where the alternative hypothesis is that the treatment effect in arm=1 is not equal to that in arm=0 with respect to survival.
#' Default is \code{"2_side"}.
#' @return An object of class rmst2perm.
#' @return \item{point_estimate}{Estimated RMST values for arm=1, arm=0, and their difference}
#' @return \item{asymptotic_test_pval}{P-value of the asymptotic test for the difference in RMST}
#' @return \item{permutation_test_methodX_pval}{P-value of the permutation test for the difference in RMST using Method X (X is the number specified in \code{mperm})}
#' @return \item{methodX_number_applied}{The number of times Method X was applied during the permutation process. (X is the number specified in \code{mperm} except 6.) For X=1 (Method 1), this returns how many additional permutations were performed in order to obtain \code{nperm} of realizations.}
#' @return \item{method5_number_exponential_used}{The number of times the exponential distribution was used for the parametric fit during the permutations. Normally, the Weibull distribution is used for Method 5.
#' However, when the maximum likelihood estimate (MLE) for the Weibull distribution cannot be derived or the hessian of the covariance matrix for the MLE is singular, the exponential distribution will be used.}
#' @return \item{tau}{The truncation time used in the analysis}
#' @return \item{mperm}{The method used to conduct the permutation test}
#' @return \item{nperm}{The number of iterations for the resampling}
#' @return \item{asy}{The type of the asymptotic variance estimator used for the asymptotic test}
#' @return \item{test}{The type of test used in the analysis}
#' @references
#' Horiguchi M, Uno H. On permutation tests for comparing restricted mean survival time
#' with small sample from randomized trials. Statistics in Medicine 2020. doi:10.1002/sim.8565.
#' @author Miki Horiguchi, Hajime Uno
#' @examples
#' #--- sample data ---#
#' D      = rmst2perm.sample.data()
#' time   = D$time
#' status = D$status
#' arm    = D$arm
#' tau    = 34
#' mperm  = c(1:6)
#' nperm  = 100 #--This number is only for the example.
#'              #--It is recommended to specify at least 100K (default) or larger.
#' seed   = 123
#'
#' a = rmst2perm(time=time, status=status, arm=arm,
#'               tau=tau, mperm=mperm, nperm=nperm, seed=seed)
#' print(a)


#'@export
#########################################
# rmst2perm (main function)
#########################################
rmst2perm <- function(time, status, arm, tau=NULL, mperm=c(1:6), nperm=100000, seed=NULL, asy="greenwood", test="2_side"){

#====================
#--- 1: Settings ---
#====================
if(is.null(tau)){
  NOTE=(paste("Please specify a truncation time, tau."))
}

#===== tau =====
idx=arm==0; tt=time[idx]; tt0max=max(tt); ss=status[idx]; ss0max=min(ss[tt==tt0max]);
idx=arm==1; tt=time[idx]; tt1max=max(tt); ss=status[idx]; ss1max=min(ss[tt==tt1max]);

#--check--
chk0 = as.numeric(tt0max>=tau | ss0max==1)
chk1 = as.numeric(tt1max>=tau | ss1max==1)

if(chk0==1 & chk1==1 & max(tt0max, tt1max)>=tau){
  #--OK--
  NULL
}else{
  stop("The difference between RMSTs was inestimable with the observed data. A smaller tau should be specified to avoid this error.")
}


#===== Method =====
mperm = sort(unique(mperm))

#--
if(!max(mperm) %in% c(1:6)){
  stop(paste("Inappropriate numbers were specified in 'mperm'. Please specify a vector with the numbers from 1 to 6."))
}

#--method labels used in the analysis--
mlab = mperm
if(4 %in% mlab){mlab = c(mlab, 2, 3)}  #--need Methods 2 and 3 to calculate Method 4--
mlab = sort(unique(mlab))

#--Need to classify methods because Method 5 needs more strict conditions than the others--
class1  = mlab[mlab!=6]
nclass1 = length(class1)

#===== Output =====
out = list()

#===========================
#--- 2: Asymptotic test  ---
#===========================
wk_asy = rmst2_update(time=time, status=status, arm=arm, tau=tau, var.method=asy, test=test)

out$point_estimate = wk_asy$unadjusted.result[1,][c("rmst1", "rmst0", "Est.")]
names(out$point_estimate) = c("rmst1", "rmst0", "difference")

out$asymptotic_test_pval = as.numeric(wk_asy$unadjusted.result[1,]["p"])

#============================
#--- 3: Permutation tests ---
#============================
obs_rmstdiff = wk_asy$unadjusted.result[1,]["Est."]
indata  = data.frame(time=time, status=status, arm=arm)
tmppval = c()

#===== Methods 1 to 5 =====
if(sum(c(1:5) %in% class1)>0){

  #--resampling--
  out_perm=c()
  k=0

  if(!is.null(seed)){set.seed(seed)}

  while(
    if(1 %in% class1){
      sum(out_perm[,1])!=nperm
    }else{
      length(out_perm[,1])!=nperm
    }
  ){

    k=k+1
    perm_dat = shuffle.flat(indata, "arm")

    tmpout = matrix(0, 1, 3+nclass1)
    colnames(tmpout) = c("chk", "chk_converge", "chk_hessian", paste0("method", class1))

    perm_d0 = perm_dat[perm_dat$arm==0,]
    perm_d1 = perm_dat[perm_dat$arm==1,]

    #--check at risk at tau--
    perm_d0_rs = sum(perm_d0$time>=tau)
    perm_d1_rs = sum(perm_d1$time>=tau)

    #--check if the method is needed--
    idx0 = as.numeric(perm_d0_rs<1 & min(perm_d0$status[perm_d0$time==max(perm_d0$time)])==0)
    idx1 = as.numeric(perm_d1_rs<1 & min(perm_d1$status[perm_d1$time==max(perm_d1$time)])==0)

    idx = sum(idx0,idx1)

    if(idx==0){
      tmpout[,"chk"] = 1
    }

    #--needed for Methods 1 to 5--
    b = rmst2_update(time=perm_dat$time, status=perm_dat$status, arm=perm_dat$arm, tau=tau, test=test)
    tmpout[,paste0("method", class1)] = b$unadjusted.result[1]

    #--needed for methods 3 to 5--
    if(tmpout[,"chk"]==0 & sum(c(3,4,5) %in% class1)>0){
      perm_tmp0 = perm_d0
      perm_tmp1 = perm_d1

      #--arm0--
      if(idx0==1){
        if(3 %in% class1){
          perm_tmp0[perm_tmp0$time==max(perm_tmp0$time),]$status = 1
        }

        if(5 %in% class1){
          if(sum(perm_d0$status)!=0){
            #--mle (exponential) to obtain the start value--
            tmpmle = expmle1(perm_d0$time, perm_d0$status)

            #--mle (weibull)--
            mle = weimle1(perm_d0$time, perm_d0$status, start_b=-log(tmpmle$lambda))
            tmpfunc <- function(x){1-pweibull(x, shape=mle$shape, scale=mle$scale)}
            tmprmst0 = integrate(tmpfunc, lower=0, upper=tau)$value
          }else{
            #--if there is no event in arm0--
            tmprmst0 = tau
          }
          tmprmst1 = rmst1_update(perm_d1$time, perm_d1$status, tau)$rmst[1]
        }
      }

      #--arm1--
      if(idx1==1){
        if(3 %in% class1){
          perm_tmp1[perm_tmp1$time==max(perm_tmp1$time),]$status = 1
        }

        if(5 %in% class1){
          if(sum(perm_d1$status)!=0){
            #--mle (exponential) to obtain the start value--
            tmpmle = expmle1(perm_d1$time, perm_d1$status)

            #--mle (weibull)--
            mle = weimle1(perm_d1$time, perm_d1$status, start_b=-log(tmpmle$lambda))
            tmpfunc <- function(x){1-pweibull(x, shape=mle$shape, scale=mle$scale)}
            tmprmst1 = integrate(tmpfunc, lower=0, upper=tau)$value
          }else{
            #--if there is no event in arm1--
            tmprmst1 = tau
          }
          tmprmst0 = rmst1_update(perm_d0$time, perm_d0$status, tau)$rmst[1]
        }
      }


      #--
      perm_tmp = rbind(perm_tmp0, perm_tmp1)
      d = rmst2_update(time=perm_tmp$time, status=perm_tmp$status, arm=perm_tmp$arm, tau=tau, test=test)

      if(3 %in% class1){
        tmpout[,"method3"] = d$unadjusted.result[1]
      }

      if(4 %in% class1){
        tmpout[,"method4"] = mean(c(b$unadjusted.result[1], d$unadjusted.result[1]))
      }

      if(5 %in% class1){
        if(sum(perm_d0$status)!=0 & sum(perm_d1$status)!=0){
          tmpout[,"chk_converge"] = mle$convergence
          if(mle$hessian_det==0){tmpout[,"chk_hessian"] = 1}
          if(mle$hessian_det!=0){tmpout[,"chk_hessian"] = 0}
        }else{
          tmpout[,"chk_converge"] = 0
          tmpout[,"chk_hessian"]  = 0
        }

        if(tmpout[,"chk_converge"]==0 & tmpout[,"chk_hessian"]==0){
          #-- use weibull --
          tmpout[,"method5"] = tmprmst1 - tmprmst0
        }

        if(tmpout[,"chk_converge"]!=0 | tmpout[,"chk_hessian"]==1){
          #-- use exponential --
          if(idx0==1){
            mle=expmle1(perm_d0$time, perm_d0$status)
            tmpfunc <- function(x){1-pexp(x, rate=mle$lambda)}
            tmprmst0 = integrate(tmpfunc, lower=0, upper=tau)$value
            tmprmst1 = rmst1_update(perm_d1$time, perm_d1$status, tau)$rmst[1]
          }
          if(idx1==1){
            mle=expmle1(perm_d1$time, perm_d1$status)
            tmpfunc <- function(x){1-pexp(x, rate=mle$lambda)}
            tmprmst1 = integrate(tmpfunc, lower=0, upper=tau)$value
            tmprmst0 = rmst1_update(perm_d0$time, perm_d0$status, tau)$rmst[1]
          }
          tmpout[,"method5"] = tmprmst1 - tmprmst0
        }
      }


    }
    out_perm=rbind(out_perm,tmpout)
  }


  #--p-values--
  rmstdiff = matrix(0, nperm, nclass1)
  colnames(rmstdiff) = paste0("method", class1)

    #--Method 1--
    if(1 %in% class1){
      rmstdiff[,"method1"] = out_perm[,"method1"][out_perm[,"chk"]==1]
      out["method1_number_applied"]   = sum(out_perm[,"chk"]==0)
    }

    #--Methods 2 to 5--
    idxq=c(2,3,4,5) %in% class1
    idxr=c(2,3,4,5) %in% mperm
    if(sum(idxq)>0){
      rmstdiff[,paste0("method", c(2,3,4,5)[idxq])] = out_perm[1:nperm,paste0("method", c(2,3,4,5)[idxq])]
      out[paste0("method", c(2,3,4,5)[idxr], "_number_applied")] = sum(out_perm[1:nperm,"chk"]==0)
    }

    #--Method 5--
    if(5 %in% class1){
      out[paste0("method5_number_exponential_used")]  = sum(out_perm[1:nperm,"chk_converge"]!=0 | out_perm[1:nperm,"chk_hessian"]==1)
    }

    #--
    tmppval = apply(rmstdiff, 2, func_pval, obs=obs_rmstdiff, nperm=nperm, test=test)
  }


#===== Methods 6 =====
if(6 %in% mlab){

  #--create pseudo values--
  indata2 = indata
  indata2$pseudo = rep(0, nrow(indata2))

  indata2[indata2$arm==1,]$pseudo = pseudo_rmst1(indata2[indata2$arm==1,]$time, indata2[indata2$arm==1,]$status, tau)
  indata2[indata2$arm==0,]$pseudo = pseudo_rmst1(indata2[indata2$arm==0,]$time, indata2[indata2$arm==0,]$status, tau)

  #--resampling--
  out_perm2=c()
  set.seed(seed)

  for(k in 1:nperm){
    perm_dat = shuffle.flat(indata2, "arm")

    tmpout2 = matrix(0, 1, 1)
    colnames(tmpout2) = c("method6")
    perm_d0 = perm_dat[perm_dat$arm==0,]
    perm_d1 = perm_dat[perm_dat$arm==1,]

    #--
    tmpout2[1,1] = mean(perm_d1[,"pseudo"]) -  mean(perm_d0[,"pseudo"])

    #--
    out_perm2=rbind(out_perm2,tmpout2)
  }

  colnames(out_perm2) = c("method6")
  tmppval = c(tmppval, apply(out_perm2, 2, func_pval, obs=obs_rmstdiff, nperm=nperm, test=test))
}

#===== Output =====
out[paste0("permutation_test_method", mperm, "_pval")] = tmppval[paste0("method", mperm)]

out$tau     = tau
out$mperm   = mperm
out$nperm   = nperm
out$asy     = asy
out$test    = test
class(out)  = "rmst2perm"

out
}
NULL


#' @name print.rmst2perm
#' @aliases print.rmst2perm
#' @title print.rmst2perm
#' @description S3 method for class 'rmst2perm'
#' @param x Object to be printed
#' @param digits Integer indicating the number of decimal places
#' @param ... Further arguments ignored in this function
#' @return returns summary output for class 'rmst2perm'

#' @export
######################################
# print.rmst2perm (hidden)
######################################
print.rmst2perm <- function(x, digits=3, ...){

  cat("\n")

  taus = x$tau
  cols = as.numeric(x$mperm)
  est  = x$point_estimate

  #---summary---
  tmp8 = matrix(0, (length(cols)+1), 1)
  rownames(tmp8) = c("asymptotic", paste0("method ", cols))
  colnames(tmp8) = "p-value"
  tmp8["asymptotic",]  = as.numeric(x["asymptotic_test_pval"])
  tmp8[c(paste0("method ", cols)),] = as.numeric(x[paste0("permutation_test_method", cols, "_pval")])
  result = round(tmp8, digits=digits)

  #---print out---
  cat ("<RMST estimation>", "\n")
  cat("Tau:", x$tau)
  cat("\n")
  cat("\n")
  print(est)

  cat("\n")
  cat("\n")

  cat ("<Test result> \n")
  cat("\n")
  print(result)

  invisible(x)
}
NULL
