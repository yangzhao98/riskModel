# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' @title Calculate the cumulative incidence function in the presence of competing risks
#' @description This function is used to calculate the cumulative incidence function based on the time-varying weights that corrects for right censored and left truncated data via a counting process format data. With the counting process format data, sub-distribution proportional hazards model can be performed. More details for the time-varying weights can be found at \code{\link[mstate]{crprep}}.
#'
#' @param dat Data set used in the analysis
#' @param time Time to event
#' @param status Status at the end of study
#' @param Xs Predictors used in the study
#' @param Yeoi Value of \code{status} indicate the event of interest
#' @param Ycrs Value of \code{status} indicate the competing risks
#' @param id ID for each individual
#' @param timePoints Time points of interest for calculating cumulative incidence
#'
#' @return A data frame with the estimated cumulative incidence rate at time points of interest
#' @seealso \code{\link[mstate]{crprep}}, \code{\link[survival]{survfit}}
#'
#' @export
#'
#' @examples
#' library(riskModel)
#' head(rdata)
#' getCIF(dat=rdata,time="time",
#'        status="status_num",Yeoi=1,Ycrs=2,id="id",
#'        Xs=c("age","size"),timePoints=1:5)
#'
getCIF <- function(dat,time,status,Xs,Yeoi,Ycrs,id,timePoints=c(1:5)) {
  weight.cens <- NULL
  dat <- as.data.frame(dat)
  trans <- c(Yeoi,Ycrs)
  datMS <- mstate::crprep(Tstop = time, status = status, trans= trans,
                          id= id, keep = Xs, data=dat)
  datMS.w <- datMS[datMS$failcode==Yeoi,]
  mfit <- survival::survfit(
    survival::Surv(Tstart,Tstop,status==Yeoi) ~ 1,
    data = datMS.w,
    weights = weight.cens
  )
  mfit_summary <- summary(mfit,times=timePoints,extend=TRUE)
  res <- list()
  res$mfit <- mfit
  res$CIF <- with(mfit_summary,
                  data.frame(
                    time=mfit_summary$time,
                    n.risk=round(n.risk,0),
                    n.event=round(n.event,0),
                    n.censor=round(n.censor,0),
                    survival=1-round(surv,3),
                    std.err=round(std.err,3),
                    Lower.95=round(lower,3),
                    Upper.95=round(upper,3)
                  ))
  class(res) <- "getCIF"
  return(res)
}

