#=======================================
#--- rmst1_update (one-arm) --- hidden
#========================================
rmst1_update=function(time, status, tau, alpha=0.05, var.method="greenwood"){

  ft  = survfit(Surv(time, status)~1)
  idx = ft$time<=tau

  wk.time = sort(c(ft$time[idx],tau))
  wk.surv = ft$surv[idx]
  wk.n.risk  = ft$n.risk[idx]
  wk.n.event = ft$n.event[idx]

  time.diff = diff(c(0, wk.time))
  areas = time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst

  #--- asymptotic variance ---
  if(var.method=="greenwood"){
    #--Greenwood plug-in estimator
    wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0, wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  }
  if(var.method=="aj"){
    #--Aalen-Johansen plug-in estimators
    wk.var <- ifelse( wk.n.risk==0, 0, wk.n.event /(wk.n.risk *wk.n.risk))
  }

  wk.var = c(wk.var,0)
  rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  = sqrt(rmst.var)

  #--- output ---
  out=matrix(0,2,4)
  out[1,]=c(rmst, rmst.se, rmst-qnorm(1-alpha/2)*rmst.se, rmst+qnorm(1-alpha/2)*rmst.se)
  out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out)=c("RMST","RMTL")
  colnames(out)=c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))

  Z=list()
  Z$result=out
  Z$rmst = out[1,]
  Z$rmtl = out[2,]
  Z$tau=tau
  Z$rmst.var = rmst.var
  Z$fit=ft
  class(Z)="rmst1_update"

  return(Z)
}
NULL

#=======================================
#--- rmst2_update (two-arm) --- hidden
#=======================================
rmst2_update=function(time, status, arm, tau=NULL, alpha=0.05, var.method="greenwood", test){

  Z=list()

  #--
  wk1=rmst1_update(time[arm==1], status[arm==1], tau, alpha, var.method)
  wk0=rmst1_update(time[arm==0], status[arm==0], tau, alpha, var.method)

  #--- contrast (RMST difference) ---
  rmst.diff.10     = wk1$rmst[1]-wk0$rmst[1]
  rmst.diff.10.se  = sqrt(wk1$rmst.var + wk0$rmst.var)
  rmst.diff.z      = rmst.diff.10/rmst.diff.10.se

  if(test=="1_side"){
    rmst.diff.z.1side      = rmst.diff.z
    rmst.diff.pval.1side   = 1-pnorm(rmst.diff.z.1side) # one-sided test (upper)
    rmst.diff.result.1side = cbind(rmst.diff.10, rmst.diff.10.se, 1, rmst.diff.z.1side, rmst.diff.pval.1side, wk1$rmst[1], wk0$rmst[1])
    #--
    out = rmst.diff.result.1side
  }else{
    #test=="2_side"
    rmst.diff.z.2side      = abs(rmst.diff.z)
    rmst.diff.pval.2side   = pnorm(-rmst.diff.z.2side)*2 # two-sided test
    rmst.diff.result.2side = cbind(rmst.diff.10, rmst.diff.10.se, 2, rmst.diff.z.2side, rmst.diff.pval.2side, wk1$rmst[1], wk0$rmst[1])
    #--
    out = rmst.diff.result.2side
  }

  #--- results ---
  rownames(out)=c("RMST (arm=1)-(arm=0)")
  colnames(out) = c("Est.", "S.E.", "test-side", "z", "p", "rmst1", "rmst0")

  #--- output ---
  Z$unadjusted.result = out
  class(Z)="rmst2_update"

  Z
}
NULL

#======================================
#--- pseudo_rmst1 (one-arm) --- hidden
#======================================
pseudo_rmst1 <- function(time, status, tau){
	pseudo = c()
	all    = rmst1_update(time, status, tau)$rmst[1]
	nn     = length(time)

	for(j in 1:nn){
	  if(max(time[-j])<tau & min(status[time==max(time[-j])])==0){
	    stop("To implement Method 6, the observed data from each group need to satisfy one of the three conditions described in the paper. See Horiguchi and Uno (2020) for details. Two options to avoid this error: 1) delete 6 from mperm, 2) specify a smaller tau.")
	  	pseudo[j] = NA
	  }else{
	  	pseudo[j] = nn*all - (nn-1)*rmst1_update(time[-j], status[-j], tau)$rmst[1]
	  }
	}
	return(pseudo)
}
NULL


#============================
#--- shuffle.flat --- hidden
#============================
shuffle.flat <- function(data, key.var, seed=NULL){

  if(!is.null(seed)) set.seed(seed)

  #--- figure out the patterns ---
  tmp = data[,-1]
  n   = nrow(tmp)
  k   = ncol(tmp)
  pattern = rep(0, n)
  for (i in 1:k){
    pattern = pattern + as.numeric(!is.na(tmp[,i]))*(10^(i-1))
  }

  unique_pattern = sort(unique(pattern))
  npatterns      = length(unique_pattern)

  tmp2 = cbind(data, pattern)

  #--- permuate entire (ver.1)---
  D = data
  n = nrow(D)
  key.org = D[,key.var]
  key.new = sample(key.org, size=n)
  D[,key.var] = key.new

  D
}
NULL

#===============================
#--- obtain p-value --- hidden
#===============================
func_pval <- function(x, obs, nperm, test){
  if(test=="1_side"){
    tmp = sum(x>obs)/nperm
  }else{
    tmp = sum(abs(x)>abs(obs))/nperm
  }
  tmp
}


#============================================
#--- MLE: Weibull and Exponential --- hidden
#============================================
mle_edit<-function(minuslogl, start = formals(minuslogl), method = "BFGS",
                   fixed = list(), nobs, ...){
  call <- match.call()
  n <- names(fixed)
  fullcoef <- formals(minuslogl)
  if (any(!n %in% names(fullcoef)))
    stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
  fullcoef[n] <- fixed
  if (!missing(start) && (!is.list(start) || is.null(names(start))))
    stop("'start' must be a named list")
  start[n] <- NULL
  start <- sapply(start, eval.parent)
  nm <- names(start)
  oo <- match(nm, names(fullcoef))
  if (anyNA(oo))
    stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
  start <- start[order(oo)]
  nm <- names(start)
  f<-function(p){
    l <- as.list(p)
    names(l) <- nm
    l[n] <- fixed
    do.call("minuslogl", l)
  }
  oout <- if(length(start))
    optim(start, f, method = method, hessian = TRUE, ...)
  else list(par = numeric(), value = f(start))
  coef <- oout$par

  new("mle", call = call, coef = coef, details = oout)
}
NULL

#===weibull mle for 1 sample====
weimle1<-function(time, status, maxit=400, start_b){
  #--- exp(a): shape
  #--- exp(b): scale

  #--minuslog likelihood function
  ll<-function(a, b){
    -sum(status*( log(exp(a)) + (exp(a)-1)*log(time) - exp(a)*log(exp(b)) ) -(time/exp(b))^exp(a) )
  }

  #--
  est = mle_edit(minuslogl=ll, start=list(a=0, b=start_b), control=c(maxit=maxit))
  shape = as.numeric(exp(attributes(est)$coef[1]))
  scale = as.numeric(exp(attributes(est)$coef[2]))
  convergence = attributes(est)$details$convergence
  hessian_det = det(attributes(est)$details$hessian)
  return(list(out=est, shape=shape, scale=scale, convergence=convergence, hessian_det=hessian_det))
}
NULL


#===exponential mle for 1 sample====
expmle1<-function(time, status){
  lambda=sum(status)/sum(time)
  scale=1/lambda
  return(list(lambda = lambda, scale = scale))
}
NULL
