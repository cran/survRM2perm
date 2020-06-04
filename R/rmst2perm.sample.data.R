#' @name rmst2perm.sample.data
#' @aliases  rmst2perm.sample.data
#' @title Sample Dataset from Ovarian Data
#' @description Generates a sample dataset of 26 randomized patients from the ovarian data.
#' @usage rmst2perm.sample.data(t.unit="month")
#' @param t.unit Specify the time unit. It supports \code{"month"} (default) and \code{"day"}.
#' @details The function creates a sample dataset to illustrate the usage of the function \code{rmst2perm()} in this package.
#' This function loads the ovarian data from the survival package, deriving three variables.
#' The variables in the dataset are as follows: \code{time}, survival time in months;
#' \code{status}, event indicator (0=censor, 1=event); \code{arm}, treatment arm (0=cyclophosphamide, 1=cyclophosphamide+adriamycin).
#' @return returns a data frame
#' @references Collett, D. Modelling Survival Data in Medical Research. Chapman and Hall/CRC. 2015; page 213.
#'
#' Edmonson JH, Fleming TR, Decker DG, et al. Different chemotherapeutic sensitivities and host factors affecting prognosis in
#' advanced ovarian carcinoma versus minimal residual disease. Cancer Treat Rep. 1979;63(2):241-247.
#' @seealso \code{ovarian} in survival package
#' @examples
#' D = rmst2perm.sample.data()
#' head(D)


#' @export
#######################################
# rmst2perm sample data
######################################
rmst2perm.sample.data=function(t.unit="month"){
  tmp <- survival::ovarian

  D=tmp[,c("futime", "fustat", "rx")]

  if(t.unit=="month"){
    D$time=D$futime/365.25*12
  }else{
    D$time=D$futime
  }

  D$status=D$fustat
  D$arm=as.numeric(D$rx==2)
  DA=D[,c("time", "status", "arm")]
  DA

}
NULL
