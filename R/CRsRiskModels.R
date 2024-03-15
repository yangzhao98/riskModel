#' @title Competing risk prediction models' evaluation
#' @param cox.fit The fitted model using either \code{rms::cph()} or \code{riskRegression::CSC()}
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param nboot The number of bootstrapped samples
#' @export
CRsModelEvaluation <- function(cox.fit,dat,time,status,Yeoi,predictYear,nboot=100) {
  # dat <- as.data.frame(dat)
  ## get C-index and ICI
  cindex <- getR2DCindexICI(cox.fit=cox.fit,
                            dat=dat,time=time,status=status,Yeoi=Yeoi,
                            predictYear=predictYear,nboot=nboot)
  ## get calibration slope and calibration-in-the-large
  cal <- getCalSlopeInt(cox.fit=cox.fit,
                        dat=dat,time=time,status=status,Yeoi=Yeoi,
                        predictYear=predictYear)
  ## get observed and expected ratio
  oer <- getObsExpRatio(cox.fit=cox.fit,
                        dat=dat,time=time,status=status,Yeoi=Yeoi,
                        predictYear=predictYear)
  oer_OERatio <- oer$dat.ObsExpRatio
  datRes <- rbind(cindex,cal,oer_OERatio[oer_OERatio$Measure !="Avg predicted risk",])
  datRes$idx <- 2
  datRes$idx[grepl("R2",datRes$Measure)] <- 1
  datRes$idx[grepl("Scaled Brier score",datRes$Measure)] <- 3
  datRes$idx[grepl("C-index",datRes$Measure)] <- 4
  datRes$idx[grepl("Time-dependent AUC",datRes$Measure)] <- 5
  datRes$idx[grepl("Royston D",datRes$Measure)] <- 6
  datRes$idx[grepl("Avg observed risk",datRes$Measure)] <- 7
  datRes$idx[grepl("Avg predicted risk",datRes$Measure)] <- 8
  datRes$idx[grepl("Observed/Expected ratio",datRes$Measure)] <- 9
  datRes$idx[grepl("Calibration-in-the-large",datRes$Measure)] <- 10
  datRes$idx[grepl("Calibration slope",datRes$Measure)] <- 11
  datRes$idx[grepl("Integrated calibration index",datRes$Measure)] <- 12
  datRes$idx[grepl("E10",datRes$Measure)] <- 13
  datRes$idx[grepl("E50",datRes$Measure)] <- 14
  datRes$idx[grepl("E90",datRes$Measure)] <- 15
  datRes$idx[grepl("Emax",datRes$Measure)] <- 16
  datRes$idx[grepl("Root squared bias",datRes$Measure)] <- 17
  datRes <- datRes[order(datRes$idx),]
  row.names(datRes) <- NULL
  datRes <- datRes[,c("Measure","Estimate","Lower.95","Upper.95")]
  res <- list()
  res$dat.ModelPerformance <- datRes
  res$dat.ModelPrediction <- oer$dat.ModelPrediction
  res
}

#' @title Competing risk prediction models with predictors selections
#' @param Xs.list Candidate predictors list considered in developing risk models
#' @param datTrain Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Yeoi Event value of interest
#' @param predictYear Predicted year of developing disease's risk
#' @param nboot The number of bootstrapped samples
#' @param datValid Data set or the validation cohort
#' @export
CRsModelSelection <- function(Xs.list,
                              datTrain,time,status,Yeoi,predictYear,nboot=100,datValid) {
  model.X <- names(Xs.list)
  datResult <- do.call(
    "rbind",
    lapply(1:length(Xs.list), FUN=function(i) {
      ## Construct the fitted model
      fitModel <- as.formula(
        paste("Hist(",time,",",status,") ~",
              paste(Xs.list[[i]],collapse="+"), sep=""))
      ## Fit the model
      riskModel <- riskRegression::FGR(fitModel,data=datTrain,cause=Yeoi)
      ## Model evaluation
      capture.output(
        datModelValid <- CRsModelEvaluation(
          cox.fit=riskModel,dat=datValid,time=time,status=status,Yeoi=Yeoi,
          predictYear=predictYear, nboot=nboot)$dat.ModelPerformance
      )
      datModelValid$Model <- model.X[i]
      datModelValid
    }))
  datResult
}
