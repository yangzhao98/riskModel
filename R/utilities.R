#' @title Calculate R2
#' @param lp Linear predictors
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param ties Allow for tied time. By default, \code{ties=TRUE}
#' @keywords internal
#' @export
Rsq <- function(lp, time, status, ties = TRUE) {
  # https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Functions/Rsquared.R
  # Function to calculate R2 and also Royston's D
  # @author Terry Therneau, Daniele Giardiello
  pi <- 3.141592653589793 # a user might have overwritten the constant
  phat <- lp
  y2 <- survival::Surv(time, status)
  n <- length(phat)

  if (ties && any(duplicated(phat))) {
    z <- stats::qnorm((1:n - 3 / 8) / (n + .25))
    # per the paper, take the average z over any ties
    index1 <- match(phat, sort(unique(phat)))
    index2 <- rank(phat, ties.method = "first")
    z2 <- tapply(z[index2], index1, mean)
    qhat <- z2[index1]
  }
  else {
    qhat <- qnorm((rank(phat) - 3 / 8) / (n + .25))
  } # simple case of no ties

  rfit <- coxph(y2 ~ qhat)
  beta <- unname(coef(rfit))
  D <- beta * sqrt(8 / pi)
  Dse <- sqrt(rfit$var[1,1]) * sqrt(8/pi)
  R2 <- beta^2 / ((pi^2 / 6) + beta^2)
  c(D = D, Dse=Dse, R2 = R2) # return vector
}

#' @title Calculate the brier score
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param predictYear Predicted t year of disease's risk
#' @param survival Predicted survival at time t
#' @keywords internal
#' @export
brier_score <- function(time, status, predictYear, survival) {
  # https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Functions/brier_score.R
  # Calculate Brier score for survival data without development data and model
  # @author Daniele Giardiello
  # library(dplyr)
  # library(riskRegression)
  # library(pec)

  db <- data.frame(time = time,status = status)

  db$predsurv <- survival
  db$predmort <- 1 - db$predsurv

  db <- db %>% dplyr::mutate(cat = case_when(
    time <= predictYear & status == 1 ~ 1,
    time > predictYear ~ 2,
    time <= predictYear & status == 0 ~ 3
  ))


  sfit <- db %>% survival::survfit(survival::Surv(time, status == 0) ~ 1, data = .)
  sfit_df <- data.frame(
    time = c(0, sfit$time),
    surv = c(1, sfit$surv),
    weight = c(1, 1 / sfit$surv)
  )

  db2 <- db %>%
    dplyr::left_join(sfit_df, by = "time") %>%
    dplyr::select(time, status, predsurv, predmort, cat, weight)

  db2$weight[db2$time > predictYear] <- max(db2$weight[db2$weight != max(db2$weight)])
  # summary(data2$weight)
  # tail(data2)

  db2$weight[db2$cat == 3] <- 0
  db2$contrib[db2$cat == 1] <- (-db$predsurv[db2$cat == 1])**2
  db2$contrib[db2$cat == 2] <- (1 - db$predsurv[db2$cat == 2])**2
  db2$contrib[db2$cat == 3] <- 0
  db2$bs <- db2$contrib * db2$weight

  brier <- (1 / sum(db2$weight)) * sum(db2$bs)

  # IPA
  # Null model
  cox_null <- db %>% survival::coxph(
    survival::Surv(time, status) ~ 1, data = ., x = T, y = T)

  db_null <- db2
  db_null$predsurv_null <- db_null %>%
    pec::predictSurvProb(cox_null, times = predictYear, newdata = .)

  sfit_null <- db %>%
    survival::survfit(survival::Surv(time, status) ~ 1, data = .)
  sfit_null_df <- data.frame(
    time = c(0, sfit_null$time),
    surv = c(1, sfit_null$surv)
  )

  db_null2 <- db_null %>%
    dplyr::left_join(sfit_null_df, by = "time") %>%
    dplyr::select(time, status, predsurv_null, predmort, weight, cat, surv)

  db_null2$weight[db_null2$cat == 3] <- 0
  db_null2$contrib[db_null2$cat == 1] <-
    (-db_null2$predsurv_null[db_null2$cat == 1])**2
  db_null2$contrib[db_null2$cat == 2] <-
    (1 - db_null2$predsurv_null[db_null2$cat == 2])**2
  db_null2$contrib[db_null2$cat == 3] <- 0
  db_null2$bs <- db_null2$contrib * db_null2$weight

  brier_null <- (1 / sum(db_null2$weight)) * sum(db_null2$bs)
  IPA <- 1 - (brier / brier_null)

  res <- c(brier, brier_null, IPA)
  names(res) <- c("Brier", "Null Brier", "IPA")
  return(res)
}

#' @title Conduct the GND Chi-squre test
#' @param pred Predicted risk
#' @param tvar Survival time
#' @param out Value for recording the event
#' @param cens.t Value for recording the censoring
#' @param groups Groups assigned to each observation
#' @param adm.cens End of follow-up time
#' @keywords internal
#' @export
GNDCalibrate = function(pred, tvar, out, cens.t, groups, adm.cens){
  # Source: https://ncook.bwh.harvard.edu/assets/GND_w_practical_example.R
  kmdec=function(dec.num,dec.name, datain, adm.cens){
    stopped=0
    data.sub=datain[datain[,dec.name]==dec.num,]
    if (sum(data.sub$out)>1){
      avsurv=survfit(Surv(tvar,out) ~ 1, data=datain[datain[,dec.name]==dec.num,], error="g")
      avsurv.est=ifelse(min(avsurv$time)<=adm.cens,
                        avsurv$surv[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],
                        1)

      avsurv.stderr=ifelse(min(avsurv$time)<=adm.cens,
                           avsurv$std.err[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],
                           0)
      avsurv.stderr=avsurv.stderr*avsurv.est

      avsurv.num=ifelse(min(avsurv$time)<=adm.cens,
                        avsurv$n.risk[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],
                        0)

    } else {
      return(c(0,0,0,0,stopped=-1))
    }

    if (sum(data.sub$out)<5) stopped=1
    c(avsurv.est, avsurv.stderr, avsurv.num, dec.num, stopped)
  }#kmdec

  tvar.t=ifelse(tvar>adm.cens, adm.cens, tvar); tvar.t
  out.t=ifelse(tvar>adm.cens, 0, out); out.t

  datause=data.frame(pred=pred, tvar=tvar.t,
                     out=out.t, count=1,
                     cens.t=cens.t, dec=groups)
  numcat=length(unique(datause$dec));numcat
  groups=sort(unique(datause$dec));groups

  kmtab=matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),
               ncol=5, byrow=TRUE)

  # if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  # else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\
  #                                       (see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
  #                                       Consider collapsing some groups to avoid this problem.")

  hltab=data.frame(group=kmtab[,4],
                   totaln=tapply(datause$count,datause$dec,sum),
                   censn=tapply(datause$cens.t,datause$dec,sum),
                   numevents=tapply(datause$out,datause$dec,sum),
                   expected=tapply(datause$pred,datause$dec,sum),
                   kmperc=1-kmtab[,1],
                   kmvar=kmtab[,2]^2,
                   kmnrisk=kmtab[,3],
                   expectedperc=tapply(datause$pred,datause$dec,mean))

  hltab$kmnum=hltab$kmperc*hltab$totaln
  hltab$GND_component=ifelse(hltab$kmvar==0,
                             0,
                             (hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))

  list(datTable=hltab,
       df=numcat-1,
       chi2gw=sum(hltab$GND_component,na.rm=TRUE),
       pvalgw=1-pchisq(sum(hltab$GND_component,na.rm=TRUE),numcat-1))

}

#' @title Calculate the calibration measures
#' @param dat Data set
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param riskPred Predicted risk at t-year follow-up
#' @param predictYear Predicted year of developing disease's risk
#' @keywords internal
#' @export
getCal <- function(dat,time,status,riskPred,predictYear) {
  dat$cll_pred <- log(-log(1-riskPred))
  modelFit <- as.formula(paste("Surv(",time,",",status,")~rcs(cll_pred,3)",sep=""))
  calibrate.cox <- coxph(modelFit,x=TRUE,data=dat)
  predict.grid.cox <- seq(quantile(riskPred,probs=0.01),
                          quantile(riskPred,probs=0.99),length=100)
  predict.grid.cox.cll <- data.frame(cll_pred=log(-log(1-predict.grid.cox)))
  predict.calibrate.cox <- 1 - pec::predictSurvProb(
    calibrate.cox,newdata=predict.grid.cox.cll,times=predictYear)
  dat.cal <- data.frame(predicted=predict.grid.cox,
                        observed=predict.calibrate.cox,
                        ID=seq(1:100))
  ## Naive method
  tfit <- lm(observed~predicted,data=dat.cal)
  fit_int <- coef(tfit)[1]
  fit_int.se <- sqrt(diag(vcov(tfit)))[1]
  fit_slope <- coef(tfit)[2]
  fit_slope.se <- sqrt(diag(vcov(tfit)))[2]
  res <- data.frame(
    Measure=c("Calibration-in-the-large","Calibration slope"),
    Estimate=c(round(fit_int,3),
               round(fit_slope,3)),
    Lower.95=c(round(fit_int-qnorm(1-0.05/2)*fit_int.se,3),
               round(fit_slope-qnorm(1-0.05/2)*fit_slope.se,3)),
    Upper.95=c(round(fit_int+qnorm(1-0.05/2)*fit_int.se,3),
               round(fit_slope+qnorm(1-0.05/2)*fit_slope.se,3)))
  # ## Recommended method
  # fit_cal_int <- geepack::geese(observed ~ offset(predicted),
  #                               data = dat.cal, id = ID,
  #                               scale.fix = TRUE, family = gaussian,
  #                               mean.link = "cloglog",
  #                               corstr = "independence", jack = TRUE)
  # fit_cal_slope <- geepack::geese(observed ~ offset(predicted)+predicted,
  #                                 data = dat.cal, id = ID,
  #                                 scale.fix = TRUE, family = gaussian,
  #                                 mean.link = "cloglog",
  #                                 corstr = "independence", jack = TRUE)
  # fit_int <- summary(fit_cal_int)$mean; fit_int
  # fit_slope <- summary(fit_cal_slope)$mean["predicted",]; fit_slope
  # res <- data.frame(
  #     Measure=c("Calibration-in-the-large","Calibration slope"),
  #     Estimate=c(round(fit_int$estimate,3),
  #                round(1+fit_slope$estimate,3)),
  #     Lower.95=c(round(fit_int$estimate-qnorm(1-0.05/2)*fit_int$san.se,3),
  #                round(1+(fit_slope$estimate-qnorm(1-0.05/2)*fit_slope$san.se),3)),
  #     Upper.95=c(round(fit_int$estimate+qnorm(1-0.05/2)*fit_int$san.se,3),
  #                round(1+(fit_slope$estimate+qnorm(1-0.05/2)*fit_slope$san.se),3)))
  res
}


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


#' @title data scaling
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Xs Candidate predictors
#' @param newdat New data set or the external validation cohort
#' @export
datZScale <- function(dat,time,status,Xs,newdat) {
  ## predictors considered
  nameList <- names(dat); nameList
  idxXs <- nameList %in% Xs
  idxSurv <- nameList %in% c(time,status)
  meanXs <- apply(datTrain[,idxXs],2,mean); meanXs
  sdXs <- apply(datTrain[,idxXs],2,sd); sdXs
  ZScale <- function(x) (x-mean(x))/sd(x)
  datTrain <- cbind(dat[,idxSurv],apply(dat[,idxXs],2,ZScale))
  datTmp <- newdat[,idxXs]
  xTestZ <- as.data.frame(
    do.call("cbind",
            lapply(1:length(Xs),
                   FUN=function(i) (datTmp[,i]-meanXs[i])/sdXs[i] )))
  names(xTestZ) <- Xs
  datTest <- cbind(newdat[,idxSurv],xTestZ)
  res <- list()
  res$meanXs <- meanXs
  res$sdXs <- sdXs
  res$datTrain <- datTrain
  res$datValid <- datTest
  res
}


#' @title Calculate reclassification measures of NRI and IDI
#' @param dat Data set used
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param predrisk1 Predicted risk from Model 1
#' @param predrisk2 Predicted risk from Model 2
#' @param cutoff Cutoffs used for caculating NRI and IDI
#' @export
getNRIIDI <- function(dat, time, status, predrisk1, predrisk2, cutoff) {
  nameList <- names(dat)
  cOutcome <- nameList %in% status
  c1 <- cut(predrisk1, breaks = cutoff, include.lowest = TRUE, right = FALSE)
  c2 <- cut(predrisk2, breaks = cutoff, include.lowest = TRUE, right = FALSE)
  tabReclas <- table(`Initial Model` = c1, `Updated Model` = c2)
  # cat(" _________________________________________\n")
  # cat(" \n     Reclassification table    \n")
  # cat(" _________________________________________\n")
  ta <- table(c1, c2, dat[, cOutcome])
  # cat("\n Outcome: absent \n  \n")
  TabAbs <- ta[, , 1]
  tab1.Outcome_absent <- cbind(
    TabAbs,
    ` % reclassified` =  round((rowSums(TabAbs) - diag(TabAbs))/rowSums(TabAbs), 2) * 100)
  names(dimnames(tab1.Outcome_absent)) <- c("Initial Model", "Updated Model")
  # print(tab1)
  # cat("\n \n Outcome: present \n  \n")
  TabPre <- ta[, , 2]
  tab2.Outcome_present <- cbind(
    TabPre,
    ` % reclassified` = round((rowSums(TabPre) - diag(TabPre))/rowSums(TabPre), 2) * 100)
  names(dimnames(tab2.Outcome_present)) <- c("Initial Model", "Updated Model")
  # print(tab2)
  # cat("\n \n Combined Data \n  \n")
  Tab <- tabReclas
  tab.Outcome_combined <- cbind(
    Tab,
    ` % reclassified` = round((rowSums(Tab) - diag(Tab))/rowSums(Tab), 2) * 100)
  names(dimnames(tab.Outcome_combined)) <- c("Initial Model", "Updated Model")
  # print(tab)
  # cat(" _________________________________________\n")
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  x <- improveProb(x1 = as.numeric(c11) * (1/(length(levels(c11)))),
                   x2 = as.numeric(c22) * (1/(length(levels(c22)))),
                   y = dat[, cOutcome])
  y <- improveProb(x1 = predrisk1, x2 = predrisk2, y = dat[, cOutcome])
  datResult <- data.frame(
    Measure=c("NRI (Categorical)", "NRI (Continuous)", "IDI"),
    Estimate=c(round(x$nri, 4), round(y$nri, 4), round(y$idi, 4)),
    Lower.95=c(round(x$nri - 1.96 * x$se.nri, 4),
               round(y$nri - 1.96 * y$se.nri, 4),
               round(y$idi - 1.96 * y$se.idi, 4)),
    Upper.95=c(round(x$nri + 1.96 * x$se.nri, 4),
               round(y$nri + 1.96 * y$se.nri, 4),
               round(y$idi + 1.96 * y$se.idi, 4)),
    pvalue=c(round(2 * pnorm(-abs(x$z.nri)), 5),
             round(2 * pnorm(-abs(y$z.nri)), 5),
             round(2 * pnorm(-abs(y$z.idi)), 5)))
  res <- list()
  res$dat.ReclassTable_Outcome_absent <- tab1.Outcome_absent
  res$dat.ReclassTable_Outcome_present <- tab2.Outcome_present
  res$dat.ReclassTable_Outcome_combined <- tab.Outcome_combined
  res$dat.Result <- datResult
  res
}


#' @title Check the non-linear association of predictor on the outcome
#' @param dat Data set used or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param x The continuous predictor of interest
#' @param Yeoi Event value of the outcome of interest
#' @export
checkNonlinearX <- function(dat,time,status,x,Yeoi) {
  par <- axis <- as.formula <- NULL
  ## original form
  fitCox <- as.formula(paste("Surv(",time,",",status,"==",Yeoi,")~",x,sep=""))
  fitCoxNL <- as.formula(paste("Surv(",time,",",status,"==",Yeoi,")~",
                               paste("rms::rcs(",x,",3)",collapse="+",sep=""),
                               sep=""))
  CoxM1 <- rms::cph(fitCox,x=TRUE,y=TRUE,surv=TRUE,data=dat)
  CoxM2 <- rms::cph(fitCoxNL,x=TRUE,y=TRUE,surv=TRUE,data=dat)
  ## Calculate the AIC and perform log-likelihood test
  pval <- rms::lrtest(CoxM1,CoxM2)
  aicRes <- c(stats::AIC(CoxM1),stats::AIC(CoxM2),pval$stats[1],pval$stats[3])
  names(aicRes) <- c("LinearModel","NonLinearModel","likhoodTest","pval")
  ## Predict
  dd <- rms::datadist(dat)
  options(datadist="dd")
  rmsPredCoxM1 <- rms::Predict(CoxM1)
  rmsPredCoxM2 <- rms::Predict(CoxM2)
  ## Plot the non-linear effects
  par(xaxs="i",yaxs="i",las=1)
  plot(rmsPredCoxM2[[x]],rmsPredict$yhat,
       type="l",lwd=2,col="blue",bty="n",
       xlab=stringr::str_to_title(x), ylab="log Relative Hazard")
  polygon(c(rmsPredict[[x]],rev(rmsPredict[[x]])),
          c(rmsPredict$lower,rev(rmsPredict$upper)),
          col="grey70",border=FALSE)
  par(new=TRUE)
  plot(rmsPredCoxM2[[x]],rmsPredict$yhat,
       type="l",lwd=2,col="blue",bty="n",
       xlab=stringr::str_to_title(x), ylab="log Relative Hazard")

}


#' @title Plot the predicted vs observed t-year risk
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set used
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param x The candidate predictor of interest
#' @param xlim The x-axis
#' @param factorX Indicator of categorical predictor. By default, \code{factorX=FALSE}.
#' @export
plotXsRisk <- function(cox.fit,dat,Yeoi,predictYear,x,xlim=NULL,factorX=FALSE) {
  dat$pred <- riskRegression::predictRisk(cox.fit,cause=Yeoi,
                                          newdata=dat,
                                          times=predictYear)
  par(xaxs="i",yaxs="i",las=1)
  oldpar <- par(mar=c(5.1,5.1,2.8,1.2))
  if(factorX) {
    plot(dat[[x]],dat$pred,
         bty="n", xlab=stringr::str_to_title(x),
         ylim=c(0,1), ylab="Estimated risk")
  } else {
    plot(dat[[x]],dat$pred,
         bty="n", xlim=xlim, xlab=stringr::str_to_title(x),
         ylim=c(0,1), ylab="Estimated risk")
    lines(
      lowess(dat[[x]], dat$pred),
      col = "red",
      lwd = 2
    )
  }
  par(oldpar)
}


#' @title Plot the smoothed calibration plot
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param maxRisk The maximum value of the predicted risk
#' @export
plotLOESSCalibration <- function(cox.fit,
                                 dat,time,status,
                                 Yeoi,predictYear,maxRisk=NULL) {
  baseModel <- as.formula(paste("Hist(",time,",",status,")~1",sep=""))
  capture.output(
    dat.pseudo <- data.frame(
      Score(list("Model" = cox.fit),
            formula = baseModel,
            cens.method = "km", data = dat,
            cpnf.int = TRUE,
            times= predictYear, summary = c("ipa"),
            cause = Yeoi,
            plots = "calibration")$Calibration$plotframe
    )
  )

  capture.output(
    pred_dat <- riskRegression::predictRisk(cox.fit,
                                            cause=Yeoi,
                                            times=predictYear,
                                            newdata=dat)
  )

  if(is.null(maxRisk)) {
    maxRisk <- pmin(ceiling(round(max(dat.pseudo$risk),1)*10+1)/10,1)
  } else {
    maxRisk <- maxRisk
  }
  dat.pseudo <- dat.pseudo[order(dat.pseudo$risk),]
  smooth_pseudos <- predict(stats::loess(pseudovalue~risk,
                                         data=dat.pseudo,
                                         degree=1,span=0.33),
                            se=TRUE)
  spike_bounds <- c(-0.025, 0)
  bin_breaks <- seq(0, maxRisk, length.out=100+1)
  freqs <- table(cut(pred_dat,breaks=bin_breaks)); freqs
  bins <- bin_breaks[-1]; bins
  freqs_valid <- freqs[freqs>0]; freqs_valid
  freqs_rescaled <- spike_bounds[1] + (spike_bounds[2]-spike_bounds[1]) *
    (freqs_valid - min(freqs_valid))/(max(freqs_valid)-min(freqs_valid));freqs_rescaled

  par(xaxs="i",yaxs="i",las=1)
  oldpar <- par(mar = c(5.1, 5.1, 3, 1))
  plot(dat.pseudo$risk, dat.pseudo$pseudovalue,
       xlim=c(0,maxRisk),ylim=c(spike_bounds[1],maxRisk),
       yaxt="n", frame.plot=FALSE, xaxt="n",
       xlab="Estimated risks",
       ylab="Observed outcome proportions",
       type="n")
  axis(2,seq(0,maxRisk,by=0.1),labels=seq(0,maxRisk,by=0.1))
  axis(1,seq(0,maxRisk,by=0.1),labels=seq(0,maxRisk,by=0.1))
  polygon(x=c(dat.pseudo$risk, rev(dat.pseudo$risk)),
          y=c(
            pmax(smooth_pseudos$fit-qt(0.975,smooth_pseudos$df)*smooth_pseudos$se,0),
            rev(smooth_pseudos$fit+qt(0.975,smooth_pseudos$df)*smooth_pseudos$se)),
          border=FALSE,
          col="lightgray")
  abline(a=0, b=1, col="gray", lty=2)
  lines(x=dat.pseudo$risk, y=smooth_pseudos$fit, lwd=2)
  segments(
    x0 = bins[freqs > 0], y0 = spike_bounds[1],
    x1 = bins[freqs > 0], y1 = freqs_rescaled)
  par(oldpar)
}


#' @title Plot the grouped calibration plot
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param riskThreshold User-defined thresholds used for calculating calibration measures
#' @param maxRisk The maximum value of the predicted risk
#' @export
plotGroupCalibration <- function(cox.fit,dat,
                                 time,status,Yeoi,predictYear,
                                 riskThreshold,maxRisk=0.3) {
  baseModel <- as.formula(paste("Hist(",time,",",status,")~1",sep=""))
  capture.output(
    dat.pseudo <- data.frame(
      Score(list("Model" = cox.fit),
            formula = baseModel,
            cens.method = "km",
            data = dat,
            cpnf.int = TRUE,
            times= predictYear,
            summary = c("ipa"),
            cause = Yeoi,
            plots = "calibration")$Calibration$plotframe
    )
  )

  ylim <- c(0,maxRisk)
  dat.pseudo <- dat.pseudo[order(dat.pseudo$risk),]
  dat.pseudo$smooth_pseudo <- predict(stats::loess(pseudovalue~risk,
                                                   data=dat.pseudo, degree=1,
                                                   span=0.33),se=TRUE)$fit
  dat.pseudo$g <- with(dat.pseudo,cut(risk,breaks=riskThreshold,right=FALSE))
  glabel <- levels(dat.pseudo$g)
  glabel[length(glabel)] <- gsub(")","]",glabel[length(glabel)])
  levels(dat.pseudo$g) <- glabel;glabel
  datCal <- dat.pseudo %>%
    group_by(g) %>%
    summarise(obs=mean(smooth_pseudo),
              pred=mean(risk))
  names(datCal) <- c("Risk_groups","Observed","Predicted")
  data.table::setDT(datCal)
  datCalNew <- data.table::melt(datCal,
                                id.vars=c("Risk_groups"),
                                measure.vars=c("Observed","Predicted"),
                                variable.name=c("group"),
                                value.name="risk")
  data.table::setDF(datCalNew)
  gNames <- unique(datCalNew$Risk_groups)
  p <- ggplot(datCalNew,
              aes(x=Risk_groups,y=risk,group=group,fill=group)) +
    geom_bar(stat="identity",position="dodge") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = c(0.2,0.8)) +
    scale_y_continuous(name="Probability of Event",expand=c(0,0),
                       limits=ylim, breaks=seq(ylim[1],ylim[2],0.1)) +
    scale_fill_manual("fill",values=c("Observed"="grey","Predicted"="black")) +
    labs(x="Risk groups")
  res <- list()
  res$CalibrationPlot <- p
  res$datCalibration <- datCalNew
  res
}


#' @title Calcuate the ratio of observed and predicted risk
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @export
getObsExpRatio <- function(cox.fit,dat,time,status,Yeoi,predictYear) {
  dat$statusNew <- "2"
  dat$statusNew[dat[[status]]==0] <- "0"
  dat$statusNew[dat[[status]]==Yeoi] <- "1"
  dat$statusNew <- factor(dat$statusNew)
  fitMS <- as.formula(paste("Surv(",time,",statusNew)~1",sep=""))
  aj_dat <- summary(survfit(fitMS,data=dat),times=predictYear)
  idx <- which(aj_dat$states == "1")
  obsTRisk <- data.frame(obs=aj_dat$pstate[,idx],se=aj_dat$std.err[,idx],
                         lower=aj_dat$lower[,idx],upper=aj_dat$upper[,idx])
  ## Reference: Debray et al. (2017) doi:10.1136/bmj.i6460
  predRisk <- riskRegression::predictRisk(cox.fit,cause=Yeoi,
                                          times=predictYear,
                                          newdata=dat)
  OE <- obsTRisk$obs/mean(predRisk)
  res <- list()
  res$dat.ModelPrediction <- data.frame(time=dat[[time]],
                                        status=dat[[status]],
                                        Yeoi=(dat[[status]]==Yeoi)+0,
                                        predRisk=predRisk)
  res$dat.ObsExpRatio <-  data.frame(
    Measure=c("Avg observed risk","Avg predicted risk","Observed/Expected ratio"),
    Estimate=c(round(obsTRisk$obs,3),round(mean(predRisk),3),round(OE,3)),
    Lower.95=c(round(obsTRisk$lower,3),"NA",
               round(exp(log(OE)-qnorm(1-0.05/2) *obsTRisk$se/obsTRisk$obs ),3)),
    Upper.95=c(round(obsTRisk$upper,3),"NA",
               round(exp(log(OE)+qnorm(1-0.05/2) *obsTRisk$se/obsTRisk$obs ),3)))
  res
}


#' @title Calcuate the calibration slope
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @export
getCalSlopeInt <- function(cox.fit,dat,time,status,Yeoi,predictYear) {
  baseModel <- as.formula(paste("Hist(",time,",",status,")~1",sep=""))
  dat.pseudo <- data.frame(Score(list("csh_development" = cox.fit),
                                 formula = baseModel,
                                 cens.method = "km", data = dat,
                                 cpnf.int = TRUE,
                                 times= predictYear,
                                 summary = c("ipa"),
                                 cause = Yeoi,
                                 plots = "calibration")$Calibration$plotframe)
  dat.pseudo$cll_pred <- log(-log(1-dat.pseudo$risk))
  fit_cal_int <- geepack::geese(pseudovalue ~ offset(cll_pred),
                                data = dat.pseudo, id = ID,
                                scale.fix = TRUE, family = gaussian,
                                mean.link = "cloglog",
                                corstr = "independence", jack = TRUE)
  fit_cal_slope <- geepack::geese(pseudovalue ~ offset(cll_pred)+cll_pred,
                                  data = dat.pseudo, id = ID,
                                  scale.fix = TRUE, family = gaussian,
                                  mean.link = "cloglog",
                                  corstr = "independence", jack = TRUE)
  fit_int <- summary(fit_cal_int)$mean; fit_int
  fit_slope <- summary(fit_cal_slope)$mean["cll_pred",]; fit_slope

  # # Value, confidence interval and test for calibration slope
  # summary(fit_cal_slope)
  # with(
  #   summary(fit_cal_slope)$mean["cll_pred", ],
  #   c(
  #     "slope" = 1 + estimate,
  #     `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
  #     `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
  #   )
  # )
  #
  # # Value, confidence interval and test for calibration intercept
  # summary(fit_cal_int)
  # with(
  #   summary(fit_cal_int)$mean,
  #   c(
  #     "intercept" = estimate,
  #     `2.5 %` = estimate - qnorm(0.975) * san.se,
  #     `97.5 %` = estimate + qnorm(0.975) * san.se
  #   )
  # )

  # betas <- fit_cal_slope$beta; betas
  # vcov_mat <- fit_cal_slope$vbeta; vcov_mat
  # wald <- drop(betas %*% solve(vcov_mat) %*% betas); wald
  # p <- pchisq(wald, df = 2, lower.tail = FALSE); p
  res <- data.frame(
    Measure=c("Calibration-in-the-large","Calibration slope"),
    Estimate=c(round(fit_int$estimate,3),
               round(1+fit_slope$estimate,3)),
    Lower.95=c(round(fit_int$estimate-qnorm(1-0.05/2)*fit_int$san.se,3),
               round(1+(fit_slope$estimate-qnorm(1-0.05/2)*fit_slope$san.se),3)),
    Upper.95=c(round(fit_int$estimate+qnorm(1-0.05/2)*fit_int$san.se,3),
               round(1+(fit_slope$estimate+qnorm(1-0.05/2)*fit_slope$san.se),3)))
  res
}


#' @title Calcuate the predicted risk
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @export
getPredRisk <- function(cox.fit,dat,time,status,Yeoi,predictYear){
  dat <- as.data.frame(dat)
  baseModel <- as.formula(paste("Hist(",time,",",status,")~1",sep=""))
  ## point estimate
  capture.output(
    pred_dat <- riskRegression::predictRisk(cox.fit,cause=Yeoi,
                                            times=predictYear,
                                            newdata=dat)
  )
  return(pred_dat)
}


#' @title Calcuate R2 and C-index
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param nboot The number of bootstrapped samples
#' @export
getR2DCindexICI <- function(cox.fit,dat,time,status,Yeoi,predictYear,nboot=100) {
  baseModel <- as.formula(paste("Hist(",time,",",status,")~1",sep=""))
  ## point estimate
  capture.output(
    pred_dat <- riskRegression::predictRisk(cox.fit,cause=Yeoi,
                                            times=predictYear,
                                            newdata=dat)
  )

  score_dat <- riskRegression::Score(
    list("csh_development" = cox.fit),
    formula = baseModel,
    cens.method = "km",
    data = dat,
    cpnf.int = TRUE,
    times= predictYear,
    metrics = c("auc","brier"),
    summary = c("ipa"),
    cause = Yeoi,
    plots = "calibration"
  )
  calplot_pseudo <- plotCalibration(
    x=score_dat,
    brier.in.legend = FALSE,
    auc.in.legend = FALSE,
    cens.method = "pseudo",
    bandwidth = 0.05,
    cex=1,
    plot=FALSE,
    round=FALSE)
  dat_pseudo <- calplot_pseudo$plotFrames$csh_development
  diff_dat_pseudo <- pred_dat - dat_pseudo$Obs[match(pred_dat,dat_pseudo$Pred)]
  c <- unlist(pec::cindex(object=cox.fit,
                          formula=baseModel,
                          cause=Yeoi,
                          eval.times=predictYear,
                          data=dat)$AppCindex)
  OS <- dat[[time]]
  delta <- dat[[status]]
  temp <- dat
  temp$rank <- rank(pred_dat,
                    na.last = TRUE,
                    ties.method = "first")
  Kappa <- (8 / pi) ^0.5   #kappa scaling coefficient
  rankit <- (temp$rank - 0.375) / (nrow(temp) + 0.25)
  normal_rankit <- qnorm(rankit)
  scaled_rankit <- normal_rankit / Kappa
  fg <- cmprsk::crr(ftime = OS, fstatus = delta, cov1 = scaled_rankit,
                    failcode = Yeoi, cencode = 0)
  R2 <- round(unname((fg$coef^2 / Kappa^2)/((pi^2 / 6)+(fg$coef^2 / Kappa^2))),3)
  RoystonD <- round(unname(fg$coef),3)
  pe <- c(
    R2,
    round(score_dat$Brier$score$Brier[2],3),
    round(score_dat$Brier$score$IPA[2],3),
    round(c,3),
    round(score_dat$AUC$score$AUC,3),
    round(RoystonD,3),avgPredRisk=round(mean(pred_dat),3),
    round(mean(abs(diff_dat_pseudo),na.rm=TRUE),3),
    round(quantile(abs(diff_dat_pseudo), c(0.1),na.rm=TRUE),3),
    round(quantile(abs(diff_dat_pseudo), c(0.5),na.rm=TRUE),3),
    round(quantile(abs(diff_dat_pseudo), c(0.9),na.rm=TRUE),3),
    round(max(abs(diff_dat_pseudo),na.rm=TRUE),3),
    round(sqrt(mean(diff_dat_pseudo^2,na.rm=TRUE)),3))
  ## bootstrapping for confidence interval
  set.seed(123)
  library(parallel)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {
    library(pec)
    library(riskRegression)
    library(survival)
    library(cmprsk)
  })
  clusterExport(cl,list("dat","cox.fit","time","status",
                        "baseModel","predictYear","Yeoi","nboot"),
                envir = environment())
  Boots <- setNames(parLapply(
    cl,1:nboot,
    fun=function(i) {
      tryCatch({
        set.seed(i)
        n <- nrow(dat)
        idx <- sample(1:n,size=n,replace=TRUE)
        datBoot <- dat[idx,]
        datBoot$ID <- 1:nrow(datBoot)
        pred_dat <- riskRegression::predictRisk(cox.fit,cause=Yeoi,
                                                times=predictYear,
                                                newdata=datBoot)
        score_dat <- riskRegression::Score(
          list("csh_development" = cox.fit),
          formula = baseModel,
          cens.method = "km",
          data = datBoot,
          cpnf.int = TRUE,
          times= predictYear,
          metrics = c("auc","brier"),
          summary = c("ipa"),
          cause = Yeoi,
          plots = "calibration"
        )
        calplot_pseudo <- plotCalibration(
          x=score_dat,
          brier.in.legend = FALSE,
          auc.in.legend = FALSE,
          cens.method = "pseudo",
          bandwidth = 0.05,
          cex=1,
          plot=FALSE,
          round=FALSE)
        dat_pseudo <- calplot_pseudo$plotFrames$csh_development
        diff_dat_pseudo <- pred_dat - dat_pseudo$Obs[match(pred_dat,dat_pseudo$Pred)]
        c <- unlist(pec::cindex(object=cox.fit,
                                formula=baseModel,
                                cause=Yeoi,
                                eval.times=predictYear,
                                data=datBoot)$AppCindex)
        OS <- datBoot[[time]]
        delta <- datBoot[[status]]
        temp <- datBoot
        temp$rank <- rank(pred_dat,
                          na.last = TRUE,
                          ties.method = "first")
        Kappa <- (8 / pi) ^0.5   #kappa scaling coefficient
        rankit <- (temp$rank - 0.375) / (nrow(temp) + 0.25)
        normal_rankit <- qnorm(rankit)
        scaled_rankit <- normal_rankit / Kappa
        fg <- cmprsk::crr(ftime = OS, fstatus = delta, cov1 = scaled_rankit,
                          failcode = Yeoi, cencode = 0)
        R2 <- round(unname((fg$coef^2 / Kappa^2)/((pi^2 / 6)+(fg$coef^2 / Kappa^2))),3)
        RoystonD <- round(unname(fg$coef),3)
        data.frame(
          R2=R2,
          BS=score_dat$Brier$score$Brier[2],
          ScaledBS=score_dat$Brier$score$IPA[2],
          c=c,
          timeAUC=score_dat$AUC$score$AUC,
          RoystonD=RoystonD,avgPredRisk=mean(pred_dat),
          ICI=round(mean(abs(diff_dat_pseudo),na.rm=TRUE),3),
          E10=round(quantile(abs(diff_dat_pseudo), c(0.1),na.rm=TRUE),3),
          E50=round(quantile(abs(diff_dat_pseudo), c(0.5),na.rm=TRUE),3),
          E90=round(quantile(abs(diff_dat_pseudo), c(0.9),na.rm=TRUE),3),
          Emax=round(max(abs(diff_dat_pseudo),na.rm=TRUE),3),
          RootSquaredBias=round(sqrt(mean(diff_dat_pseudo^2,na.rm=TRUE)),3))
      },warning = function(w){
        NULL
        print("warning")
      },error = function(e){
        NULL
        print("error")
      })
    }),1:nboot)
  stopCluster(cl)
  boot <- do.call("rbind",Boots)
  boot <- boot[!boot$c %in% c("NaN","-Inf","warning","error"),]
  boot[] <- lapply(boot, as.numeric)

  res <- data.frame(
    Measure=c("R2",
              "Brier score","Scaled Brier score","C-index","Time-dependent AUC",
              "Royston D","Avg predicted risk","Integrated calibration index",
              "E10","E50","E90","Emax","Root squared bias"),
    Estimate=pe,
    Lower.95=c(apply(boot,2,FUN = function(i) round(quantile(i,probs=0.025),3))),
    Upper.95=c(apply(boot,2,FUN = function(i) round(quantile(i,probs=0.975),3))))
  rownames(res) <- NULL
  res
}


#' @title Conduct decision curve analysis
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param xstart The start value of X-axis
#' @param xstop The end value of X-axis
#' @param xby The step by the X-axis
#' @param ymin The smallest value of Y-axis
#' @param probability The Probability
#' @param harm The cost of false positive results
#' @param interventionper By default, \code{interventionper=100}
#' @param smooth Indicator of whether smoothed the net benefit curve. By default, \code{smooth=FALSE}
#' @param loess.span The bin of loess()
#' @param cmprsk Indicator of whether considering competing risk. By default, \code{cmprsk=FALSE}
#' @export
getSTDCA <- function(cox.fit, dat, time, status, Yeoi, predictYear,
                     xstart = 0.01, xstop = 0.99, xby = 0.01,
                     ymin = -0.05, probability = NULL,
                     harm = NULL,
                     interventionper = 100,
                     smooth = FALSE,
                     loess.span = 0.10,
                     cmprsk = TRUE) {
  dat <- as.data.frame(dat)
  dat$pred <- riskRegression::predictRisk(cox.fit,cause=Yeoi,newdata=dat,times=predictYear)
  predRisk <- "pred"
  dat <- as.data.frame(dat)
  # ONLY KEEPING COMPLETE CASES
  dat <- dat[complete.cases(dat[c(status,time,predRisk)]),c(status,time,predRisk)]

  # outcome MUST BE CODED AS 0 AND 1
  if ((length(dat[!(dat[status] == 0 | dat[status] == 1), status]) > 0) & cmprsk == FALSE) {
    stop("outcome must be coded as 0 and 1")
  }

  # data MUST BE A DATA FRAME
  # if (class(dat) != "data.frame") { stop("Input data must be class data.frame") }
  # xstart IS BETWEEN 0 AND 1
  if (xstart < 0 | xstart > 1) { stop("xstart must lie between 0 and 1") }
  # xstop IS BETWEEN 0 AND 1
  if (xstop < 0 | xstop > 1) { stop("xstop must lie between 0 and 1") }
  # xby IS BETWEEN 0 AND 1
  if (xby <= 0 | xby >= 1) { stop("xby must lie between 0 and 1") }
  # xstart IS BEFORE xstop
  if (xstart >= xstop) { stop("xstop must be larger than xstart") }
  # STORING THE NUMBER OF PREDICTORS SPECIFIED
  pred.n <- length(predRisk)
  # IF probability SPECIFIED ENSURING THAT EACH PREDICTOR IS INDICATED AS A T OR F
  if (length(probability) > 0 & pred.n != length(probability)) {
    stop("Number of probabilities specified must be the same as the number of predictors being checked.")
  }
  # IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
  if (length(harm) > 0 & pred.n != length(harm)) {
    stop("Number of harms specified must be the same as the number of predictors being checked.")
  }
  # INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
  if (length(harm) == 0) { harm <- rep(0, pred.n) }
  if (length(probability) == 0) { probability <- rep(TRUE, pred.n) }
  # THE PREDICTOR NAMES CANNOT BE EQUAL TO all OR none.
  if (length(predRisk[predRisk == "all" | predRisk == "none"])) {
    stop("Prediction names cannot be equal to all or none.")
  }
  # CHECKING THAT EACH probability ELEMENT IS EQUAL TO T OR F,
  # AND CHECKING THAT PROBABILITIES ARE BETWEEN 0 and 1
  # IF NOT A PROB THEN CONVERTING WITH A COX REGRESSION
  for (m in 1:pred.n) {
    if (probability[m] != TRUE & probability[m] != FALSE) {
      stop("Each element of probability vector must be TRUE or FALSE")
    }
    if (probability[m] == TRUE & (max(dat[predRisk[m]]) > 1 | min(dat[predRisk[m]]) < 0)) {
      stop(paste(predRisk[m], "must be between 0 and 1 OR sepcified as a non-probability in the probability option", sep = " "))
    }
    if (probability[m] == FALSE) {
      model <- NULL
      pred <- NULL
      model <- survival::coxph(Surv(data.matrix(dat[time]), data.matrix(dat[status])) ~ data.matrix(dat[predRisk[m]]))
      surv.data <- data.frame(0)
      pred <- data.frame(1 - c(summary(survival::survfit(model, newdata = surv.data), time = predictYear, extend = TRUE)$surv))
      names(pred) <- predRisk[m]
      dat <- cbind(dat[names(dat) != predRisk[m]], pred)
      print(paste(predRisk[m], "converted to a probability with Cox regression. Due to linearity and proportional hazards assumption, miscalibration may occur.", sep = " "))
    }
  }

  #########  CALCULATING NET BENEFIT   #########
  N <- dim(dat)[1]

  # getting the probability of the event for all subjects
  # this is used for the net benefit associated with treating all patients
  if (cmprsk == FALSE) {
    km.cuminc <- survival::survfit(Surv(data.matrix(dat[time]), data.matrix(dat[status])) ~ 1)
    pd <- 1 - summary(km.cuminc, times = predictYear, extend = TRUE)$surv
  } else {
    cr.cuminc <- cmprsk::cuminc(dat[[time]], dat[[status]])
    pd <- cmprsk::timepoints(cr.cuminc, times = predictYear)$est[1]
  }

  # creating dataset that is one line per threshold for the treat all and treat none strategies;
  # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
  nb <- data.frame(seq(from = xstart, to = xstop, by = xby))
  names(nb) <- "threshold"
  interv <- nb
  error <- NULL

  nb["all"] <- pd - (1 - pd) * nb$threshold / (1 - nb$threshold)
  nb["none"] <- 0

  # CYCLING THROUGH EACH PREDICTOR AND CALCULATING NET BENEFIT
  for (m in 1:pred.n) {
    nb[predRisk[m]] <- NA

    for (t in 1:length(nb$threshold)) {
      # calculating number of true and false postives;
      px <- sum(dat[predRisk[m]] > nb$threshold[t]) / N; px

      if (px == 0) {
        error <- rbind(error, paste(predRisk[m], ": No observations with risk greater than ", nb$threshold[t] * 100, "%", sep = ""))
        break
      } else {
        # calculate risk using Kaplan Meier
        if (cmprsk == FALSE) {
          km.cuminc <- survival::survfit(Surv(data.matrix(dat[dat[predRisk[m]] > nb$threshold[t], time]),
                                              data.matrix(dat[dat[predRisk[m]] > nb$threshold[t], status])) ~ 1)
          pdgivenx <- (1 - summary(km.cuminc, times = predictYear, extend = TRUE)$surv)
          if (length(pdgivenx) == 0) {
            error <- rbind(error, paste(predRisk[m], ": No observations with risk greater than ", nb$threshold[t] * 100, "% that have followup through the timepoint selected", sep = ""))
            break
          }
          # calculate risk using competing risk
        } else {
          cr.cuminc <- cmprsk::cuminc(dat[[time]][dat[[predRisk[m]]] > nb$threshold[t]],
                                      dat[[status]][dat[[predRisk[m]]] > nb$threshold[t]])
          pdgivenx <- cmprsk::timepoints(cr.cuminc, times = predictYear)$est[1]
          if (is.na(pdgivenx)) {
            error <- rbind(error,
                           paste(predRisk[m], ": No observations with risk greater than ",
                                 nb$threshold[t] * 100, "% that have followup through the timepoint selected", sep = ""))
            break
          }
        }
        # calculating NB based on calculated risk
        nb[t, predRisk[m]] <- pdgivenx * px - (1 - pdgivenx) * px * nb$threshold[t] / (1 - nb$threshold[t]) - harm[m]
      }
    }
    interv[predRisk[m]] <- (nb[predRisk[m]] - nb["all"]) * interventionper / (interv$threshold / (1 - interv$threshold))
  }
  if (length(error) > 0) {
    print(paste(error, ", and therefore net benefit not calculable in this range.", sep = ""))
  }

  # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED
  for (m in 1:pred.n) {
    if (smooth == TRUE) {
      lws <- stats::loess(data.matrix(nb[!is.na(nb[[predRisk[m]]]), predRisk[m]]) ~ data.matrix(nb[!is.na(nb[[predRisk[m]]]), "threshold"]), span = loess.span)
      nb[!is.na(nb[[predRisk[m]]]), paste(predRisk[m], "_sm", sep = "")] <- lws$fitted

      lws <- stats::loess(data.matrix(interv[!is.na(nb[[predRisk[m]]]), predRisk[m]]) ~ data.matrix(interv[!is.na(nb[[predRisk[m]]]), "threshold"]), span = loess.span)
      interv[!is.na(nb[[predRisk[m]]]), paste(predRisk[m], "_sm", sep = "")] <- lws$fitted
    }
  }

  # RETURNING RESULTS
  results <- list()
  results$N <- N
  results$predictors <- data.frame(cbind(predRisk, harm, probability))
  names(results$predictors) <- c("predictor", "harm.applied", "probability")
  results$interventions.avoided.per <- interventionper
  results$net.benefit <- nb
  results$interventions.avoided <- interv
  return(results)
}


#' @title Plot the decision curve
#' @param getSTDCA.fit The fitted results from \code{getSTDCA()}
#' @param ylim The min and max value of Y-axis
#' @export
plotSTDCA <- function(getSTDCA.fit,ylim=c(-0.1,0.2)) {
  datNB <- getSTDCA.fit$net.benefit
  datAvoid <- getSTDCA.fit$interventions.avoided
  namePred <- names(datNB)
  namePred <- namePred[! namePred %in% c("threshold","all","none")]
  pred.n <- length(namePred)
  xmax <- max(datNB$threshold)
  xseq <- seq(0.1,xmax,0.1); xseq
  xratio <- c("1:9","1:4","3:7","2:3","1:1","3:2","7:3","4:1","9:1")
  xseqlab <- xratio[xseq*10]
  xlabNew <- paste(xseq, " (",xseqlab,")",sep="")
  xlabNew <- c("0",xlabNew)
  oldpar <- par(xaxs = "i", yaxs = "i", las = 1,
                mar = c(5.1, 5.1, 2.8, 2.1))
  plot(datNB$threshold,datNB$all,
       type="l",lwd=2,lty=1,col="grey",
       xlab="Threshold probability (Harm to benefit ratio)",
       ylab="Net benefit",
       xaxt="n",ylim=ylim,bty="n", xlim=c(0,xmax))
  axis(1,at=c(0,xseq),labels=xlabNew)
  lines(datNB$threshold,datNB$none,type="l",lwd=2,lty=2)
  for (i in 1:pred.n) {
    lines(datNB$threshold,datNB[[namePred[i]]],type="l",lwd=2,col=i)
  }
  legend("topright",c("Treat all","Treat none", paste("Risk model ",1:pred.n,sep="")),
         lwd=c(2,2,rep(2,pred.n)),
         lty=c(1,2,rep(1,pred.n)),
         col=c("grey","black",1:pred.n),
         bty="n")
  # mtext("Threshold probability \ n (Harm to benefit ratio)",1,line=2)
  # mtext("Harm to benefit ratio",1,line=6)
  par(oldpar)
}
