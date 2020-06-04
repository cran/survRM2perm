#' @name survRM2perm-package
#' @aliases  survRM2perm-package
#' @docType  package
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
#'
#' @author Miki Horiguchi, Hajime Uno
#'
#' Maintainer: Miki Horiguchi <horiguchimiki@gmail.com>
#' @references
#' Horiguchi M, Uno H. On permutation tests for comparing restricted mean survival time
#' with small sample from randomized trials. Statistics in Medicine 2020.doi:10.1002/sim.8565.
#' @keywords
#' survival
#' @seealso
#' survRM2 survival
#' @import survival
#' @import methods
#' @import stats4
#' @importFrom stats integrate optim pexp pnorm pweibull qnorm
#' @importFrom utils data
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
NULL
