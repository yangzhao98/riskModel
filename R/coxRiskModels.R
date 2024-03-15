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



#' @title Develop risk prediction models using K-fold cross-validation method
#' @param kfold The k-fold used
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Xs Candidate continuous predictors used for developing risk prediction models.
#' @param newdat New data set or the external validation cohort
#' @param predictYear The horizon year or t-year predicted risk
#' @param KMS0 Whether calculate the baseline survival rate using the Kaplan-Meier estimator
#' @param riskThreshold User-defined thresholds used for calculating calibration measures
#' @param calGroup User-defined number of groups used for calculating calibration measures
#' @param nboot The number of bootstrapped samples used for calculating the predictive performance
#' @export
kFoldCVCoxModel <- function(kfold,dat,time,status,Xs,
                            newdat,predictYear,KMS0=FALSE,
                            riskThreshold,calGroup,nboot=100) {

  dat <- data.table::setDT(dat)
  newdat <- data.table::setDT(newdat)

  ## Citation: Demler, Paynter, Cook. Tests of Calibration and Goodness of Fit
  ##   in the Survival Setting' DOI: 10.1002/sim.6428.

  if(sum(!names(dat) %in% c(time,status,Xs)) > 0 ) stop("Please check the dataset!")

  curtail <- function(x) pmin(0.999999, pmax(x, 0.000001))

  ## Conduct the K-Fold cross-validation algorithm
  if (kfold == 1) {
    modelFits <-  lapply(1:kfold, FUN=function(i) {

      train <- as.data.frame(dat)

      meanXs <- unlist(lapply(1:length(Xs), FUN=function(j) mean(train[[Xs[j]]]))); meanXs

      coxModel <- as.formula(paste("survival::Surv(",time,",",status,")~",
                                   paste(Xs,collapse="+"),sep=""))
      coxFit <- rms::cph(coxModel,data=train,x=TRUE,y=TRUE,surv=TRUE)
      S0 <- rms::survest(coxFit,times=predictYear)$surv;S0
      if(KMS0) {
        S0 <- summary(survival::survfit(
          as.formula(paste("survival::Surv(",time,",",status,")~1",sep="")),data=train),
          times=c(predictYear))$surv;
        S0
      }
      dat.S0 <- data.frame(S0=S0,fold=i)

      ## Obtain coefficient for each risk factor
      betaXs <- matrix(coef(coxFit),ncol=1)
      dat.betaXs <- data.frame(t(betaXs))
      names(dat.betaXs) <- Xs
      dat.betaXs$fold <- i

      ## Obtain meanXB
      meanXB <- as.numeric(t(betaXs) %*% matrix(meanXs,ncol=1))
      dat.meanXB <- data.frame(meanXB=meanXB,fold=i)

      ## Obtain hazard ratio for each risk factor
      dat.hrs <- data.frame(exp(cbind(coef(coxFit),confint(coxFit))))
      names(dat.hrs) <- c("HR","LCI","UCI")
      dat.hrs$Xs <- Xs; rownames(dat.hrs) <- NULL
      dat.hrs$fold <- i

      ## Obtain the C-index
      trainXs <- train[,Xs]
      predTrainY <- curtail(S0^exp( t(t(trainXs)) %*% betaXs - meanXB )[,1])
      auc.train <- Hmisc::rcorr.cens(predTrainY, Surv(train[[time]],train[[status]]))
      dat.CIndex <- data.frame(
        CIndex.train=auc.train['C Index'],CIndexSE.train=auc.train['S.D.']/2)
      dat.CIndex$fold <- i

      ## Output results
      res <- NULL
      res$fitModel <- coxFit
      res$S0 <- dat.S0
      res$beta <- dat.betaXs
      res$meanXB <- dat.meanXB
      res$HR <- dat.hrs
      res$CIndex <- dat.CIndex

      res

    })
  } else {
    nObs <- nrow(dat)
    idx <- ceiling(1:nObs/(nObs/kfold))
    set.seed(666)
    # dat$fold <- sample(rep(1:kfold, each=nrow(dat)/kfold))
    dat$fold <- sample(idx, nrow(dat))
    modelFits <-  lapply(1:kfold, FUN=function(i) {

      train <- as.data.frame(dat[dat$fold!=i,])
      test <- as.data.frame(dat[dat$fold==i,])
      # test <- dat

      meanXs <- unlist(lapply(1:length(Xs), FUN=function(j) mean(train[[Xs[j]]]))); meanXs

      coxModel <- as.formula(paste("Surv(",time,",",status,")~",paste(Xs,collapse="+"),sep=""))
      coxFit <- cph(coxModel,data=train,x=TRUE,y=TRUE,surv=TRUE)
      S0 <- survest(coxFit,times=predictYear)$surv;S0
      if(KMS0) {
        S0 <- summary(survfit(as.formula(paste("Surv(",time,",",status,")~1",sep="")),
                              data=train),times=c(predictYear))$surv; S0
      }
      dat.S0 <- data.frame(S0=S0,fold=i)

      ## Obtain coefficient for each risk factor
      betaXs <- matrix(coef(coxFit),ncol=1)
      dat.betaXs <- data.frame(t(betaXs))
      names(dat.betaXs) <- Xs
      dat.betaXs$fold <- i

      ## Obtain MeanXB
      meanXB <- as.numeric(t(betaXs) %*% matrix(meanXs,ncol=1))
      dat.meanXB <- data.frame(meanXB=meanXB,fold=i)

      ## Obtain hazard ratio for each risk factor
      dat.hrs <- data.frame(exp(cbind(coef(coxFit),confint(coxFit))))
      names(dat.hrs) <- c("HR","LCI","UCI")
      dat.hrs$Xs <- Xs; rownames(dat.hrs) <- NULL
      dat.hrs$fold <- i

      ## Obtain the C-index
      trainXs <- train[,Xs]
      predTrainY <- curtail(S0^exp( t(t(trainXs)) %*% betaXs - meanXB )[,1])
      auc.train <- rcorr.cens(predTrainY, Surv(train[[time]],train[[status]]))
      testXs <- test[,Xs]
      predTestY <- curtail(S0^exp( t(t(testXs)) %*% betaXs - meanXB )[,1])
      auc.test <- rcorr.cens(predTestY, Surv(test[[time]],test[[status]]))
      dat.CIndex <- data.frame(
        CIndex.train=auc.train['C Index'],CIndexSE.train=auc.train['S.D.']/2,
        CIndex.test=auc.test['C Index'],CIndexSE.test=auc.train['S.D.']/2)
      dat.CIndex$fold <- i

      ## Output results
      res <- NULL
      res$fitModel <- coxFit
      res$S0 <- dat.S0
      res$beta <- dat.betaXs
      res$meanXB <- dat.meanXB
      res$HR <- dat.hrs
      res$CIndex <- dat.CIndex

      res

    })

  }

  dat.CIndex <- do.call("rbind",lapply(1:kfold, FUN=function(i) modelFits[[i]]$CIndex))
  dat.beta <- do.call("rbind",lapply(1:kfold, FUN=function(i) modelFits[[i]]$beta))
  dat.HR <- do.call("rbind",lapply(1:kfold, FUN=function(i) modelFits[[i]]$HR))
  dat.S0 <- do.call("rbind",lapply(1:kfold, FUN=function(i) modelFits[[i]]$S0))
  dat.meanXB <- do.call("rbind",lapply(1:kfold, FUN=function(i) modelFits[[i]]$meanXB))
  row.names(dat.CIndex) <- row.names(dat.beta) <- NULL
  row.names(dat.HR) <- row.names(dat.S0) <- row.names(dat.meanXB) <- NULL

  ## Select the optimal model using k-fold cross-validation
  if (kfold==1) {
    opt <- with(dat.CIndex, which.max(CIndex.train)); opt
  } else {
    opt <- with(dat.CIndex, which.max(CIndex.train + CIndex.test)); opt
    # opt <- with(dat.CIndex, which.max(CIndex.test)); opt
  }

  opt.beta <- as.matrix(dat.beta[opt,1:length(Xs)],ncol=1); opt.beta
  opt.HR <- dat.HR[dat.HR$fold==opt,]; opt.HR
  opt.S0 <- dat.S0$S0[dat.S0$fold==opt]; opt.S0
  opt.meanXB <- dat.meanXB$meanXB[dat.meanXB$fold==opt]; opt.meanXB
  opt.Model <- modelFits[[opt]]$fitModel
  base.Model <- as.formula(paste("Surv(",time,",",status,")~1",sep=""))

  ## Model evaluation
  ## linear predictors and predicted survival probability
  predTrainYLP <- (t(t(dat[,Xs,with=FALSE])) %*% t(opt.beta) - opt.meanXB)[,1]
  predTrainY <- curtail(opt.S0^exp(predTrainYLP))

  ## Parallel computing
  set.seed(123)
  requireNamespace(parallel)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {
    requireNamespace(survival)
    requireNamespace(dplyr)
    requireNamespace(riskRegression)
    requireNamespace(pec)
  })
  clusterExport(cl,list("dat","Rsq","brier_score","curtail",
                        "time","status","predictYear","Xs",
                        "opt.S0","opt.beta","opt.meanXB"),
                envir = environment())
  R2DBoots <- setNames(parLapply(
    cl,1:nboot,
    fun=function(i) {
      set.seed(i)
      n <- nrow(dat)
      idx <- sample(1:n,size=n,replace=TRUE)
      datBoot <- dat[idx,]
      datBoot$ID <- 1:nrow(datBoot)
      predTrainYLPboot <- (t(t(datBoot[,Xs,with=FALSE])) %*% t(opt.beta) - opt.meanXB)[,1]
      predTrainYboot <- curtail(opt.S0^exp(predTrainYLPboot))
      TrainR2D <- Rsq(lp=predTrainYLPboot, time=datBoot[[time]], status=datBoot[[status]], ties = TRUE)
      TrainBS <- brier_score(time=datBoot[[time]], status=datBoot[[status]],
                             predictYear=predictYear, survival=predTrainYboot)
      data.frame(R2=TrainR2D[3],D2=TrainR2D[1],BS=TrainBS[1])
    }),1:nboot)
  stopCluster(cl)
  R2Dboot <- do.call("rbind",R2DBoots)

  # C-index
  aucTrainY <- rcorr.cens(predTrainY, Surv(dat[[time]],dat[[status]]))
  TrainC <- aucTrainY['C Index']; TrainC
  TrainCse <- aucTrainY['S.D.']/2; TrainCse
  # D-statistics & R2 & Brier score
  TrainR2D <- Rsq(lp=predTrainYLP,time=dat[[time]],status=dat[[status]],ties=TRUE)
  TrainBS <- brier_score(time=dat[[time]],status=dat[[status]],
                         predictYear=predictYear,survival=predTrainY)

  # Calibration
  # Calibration plot
  if (!is.null(riskThreshold)) {
    TrainYDecileGrp <- cut(1-predTrainY,breaks=riskThreshold,right=FALSE)
  } else { TrainYDecileGrp <- Hmisc::cut2(1-predTrainY,g=calGroup) }
  GNDres <- GNDCalibrate(pred=1-predTrainY,
                         tvar=dat[[time]],out=dat[[status]],
                         cens.t=predictYear, groups=TrainYDecileGrp,
                         adm.cens=predictYear)
  gNames <- levels(TrainYDecileGrp);gNames
  gNames[length(gNames)] <- gsub(")","]",gNames[length(gNames)])
  predTrainYDecile <- GNDres$datTable$expectedperc; predTrainYDecile
  obsTrainYDecile <- GNDres$datTable$kmperc; obsTrainYDecile
  # GNDTable <- GNDres$datTable
  # HLTest <- GNDCalibrate(pred=1-predTrainY,
  #                        tvar=dat[[time]],out=dat[[status]],
  #                        cens.t=predictYear, groups=Hmisc::cut2(1-predTrainY,g=10),
  #                        adm.cens=predictYear)
  # HLChi <- function(GNDTable=HLTest$datTable) {
  #    Obs.notA <- with(GNDTable, totaln-numevents)
  #    Obs.A <- with(GNDTable, numevents)
  #    Exp.notA <- with(GNDTable, totaln-expected)
  #    Exp.A <- with(GNDTable, expected)
  #    Chi <- sum( (Obs.A-Exp.A)^2/Exp.A + (Obs.notA-Exp.notA)^2/Exp.notA )
  #    pval <- 1-pchisq(Chi,length(Obs)-2)
  # }

  ## Calibration-in-the-large and calibration-slope
  trainCal <- getCal(dat=dat,time=time,status=status,
                     riskPred=1-predTrainY,predictYear=predictYear)
  datTrainModel <- data.frame(
    Measure = c("R2",
                "Brier Score",
                "C-statistic",
                "D-statistic",
                "Calibration Slope",
                "Calibration-in-the-large",
                "Greenwood-Nam-D'agostino Chisq (Lower.95=p-value)"),
    Estimate = c(round(TrainR2D[3],3),
                 round(TrainBS[1],3),
                 round(TrainC,3),
                 round(TrainR2D[1],3),
                 round(trainCal$Estimate[trainCal$Measure=="Calibration slope"],3),
                 round(trainCal$Estimate[trainCal$Measure=="Calibration-in-the-large"],3),
                 round(GNDres$chi2gw,3)),
    Lower.95 = c(round(quantile(R2Dboot$R2,probs=0.025),3),
                 round(quantile(R2Dboot$BS,probs=0.025),3),
                 round(TrainC-qnorm(1-0.05/2)*TrainCse,3),
                 round(quantile(R2Dboot$D2,probs=0.025),3),
                 round(trainCal$Lower.95[trainCal$Measure=="Calibration slope"],3),
                 round(trainCal$Lower.95[trainCal$Measure=="Calibration-in-the-large"],3),
                 round(GNDres$pvalgw,3)),
    Upper.95 = c(round(quantile(R2Dboot$R2,probs=0.975),3),
                 round(quantile(R2Dboot$BS,probs=0.975),3),
                 round(TrainC+qnorm(1-0.05/2)*TrainCse,3),
                 round(quantile(R2Dboot$D2,probs=0.975),3),
                 round(trainCal$Upper.95[trainCal$Measure=="Calibration slope"],3),
                 round(trainCal$Upper.95[trainCal$Measure=="Calibration-in-the-large"],3),
                 "NA")
  )
  datTrainModel$Group <- "Train"
  if (!is.null(riskThreshold)) {
    datTrainCal <- data.frame(
      Decile = c(rep(1:(length(riskThreshold)-1),2)),
      Risk = c(predTrainYDecile,obsTrainYDecile),
      Method=c(rep(c("Predicted","Observed"),each=(length(riskThreshold)-1))))
  } else {
    datTrainCal <- data.frame(
      Decile = c(rep(1:length(gNames),2)),
      Risk = c(predTrainYDecile,obsTrainYDecile),
      Method=c(rep(c("Predicted","Observed"),each=length(gNames))))
  }
  datTrainCal$Group <- "Train"

  # (2) For validation data set
  predTestYLP <- (t(t(newdat[,Xs,with=FALSE])) %*% t(opt.beta) - opt.meanXB)[,1]
  predTestY <- curtail(opt.S0^exp(predTestYLP))

  set.seed(123)
  requireNamespace(parallel)
  cl <- parallel::makeCluster(parallel::detectCores())
  parallel::clusterEvalQ(cl, {
    requireNamespace(survival)
    requireNamespace(dplyr)
    requireNamespace(riskRegression)
    requireNamespace(pec)
  })
  parallel::clusterExport(cl,list("newdat","Rsq","brier_score","curtail",
                        "time","status","predictYear","Xs",
                        "opt.S0","opt.beta","opt.meanXB"),
                envir = environment())
  R2DBootsTest <- setNames(parLapply(
    cl,1:nboot,
    fun=function(i) {
      n <- nrow(newdat)
      idx <- sample(1:n,size=n,replace=TRUE)
      datBoot <- newdat[idx,]
      datBoot$ID <- 1:nrow(datBoot)
      predTrainYLPboot <- (t(t(datBoot[,Xs,with=FALSE])) %*% t(opt.beta) - opt.meanXB)[,1]
      predTrainYboot <- curtail(opt.S0^exp(predTrainYLPboot))
      TrainR2D <- Rsq(lp=predTrainYLPboot, time=datBoot[[time]], status=datBoot[[status]], ties = TRUE)
      TrainBS <- brier_score(time=datBoot[[time]], status=datBoot[[status]],
                             predictYear=predictYear, survival=predTrainYboot)
      data.frame(R2=TrainR2D[3],
                 D2=TrainR2D[1],
                 BS=TrainBS[1])
    }),1:nboot)
  stopCluster(cl)
  R2DbootTest <- do.call("rbind",R2DBootsTest)

  # C-index
  aucTestY <- rcorr.cens(predTestY, Surv(newdat[[time]],newdat[[status]]))
  TestC <- aucTestY['C Index']
  TestCse <- aucTestY['S.D.']/2
  # D-statistics & R2 & Brier score
  TestR2D <- Rsq(lp=predTestYLP, time=newdat[[time]], status=newdat[[status]], ties=TRUE)
  TestBS <- brier_score(time=newdat[[time]], status=newdat[[status]],
                        predictYear=predictYear, survival=predTestY)
  # Calibration
  if (!is.null(riskThreshold)) {
    TestYDecileGrp <- cut(1-predTestY,breaks=riskThreshold,right=FALSE)
  } else { TestYDecileGrp <- Hmisc::cut2(1-predTestY,g=calGroup) }
  gNamesTest <- levels(TestYDecileGrp);gNamesTest
  gNamesTest[length(gNamesTest)] <- gsub(")","]",gNamesTest[length(gNamesTest)])
  GNDresTest <- GNDCalibrate(pred=1-predTestY,
                             tvar=newdat[[time]],out=newdat[[status]],
                             cens.t=predictYear, groups=TestYDecileGrp,
                             adm.cens=predictYear)
  predTestYDecile <- GNDresTest$datTable$expectedperc; predTestYDecile
  obsTestYDecile <- GNDresTest$datTable$kmperc; obsTestYDecile

  ## Calibration statistics
  testCal <- getCal(dat=newdat,time=time,status=status,
                    riskPred=1-predTestY,predictYear=predictYear)
  datTestModel <- data.frame(
    Measure = c("R2",
                "Brier Score",
                "C-statistic",
                "D-statistic",
                "Calibration Slope",
                "Calibration-in-the-large",
                "Greenwood-Nam-D'agostino Chisq (Lower.95=p-value)"),
    Estimate = c(round(TestR2D[3],3),
                 round(TestBS[1],3),
                 round(TestC,3),
                 round(TestR2D[1],3),
                 round(testCal$Estimate[testCal$Measure=="Calibration slope"],3),
                 round(testCal$Estimate[testCal$Measure=="Calibration-in-the-large"],3),
                 round(GNDresTest$chi2gw,3)),
    Lower.95 = c(round(quantile(R2DbootTest$R2,probs=0.025),3),
                 round(quantile(R2DbootTest$BS,probs=0.025),3),
                 round(TestC-qnorm(1-0.05/2)*TestCse,3),
                 round(quantile(R2DbootTest$D2,probs=0.025),3),
                 round(testCal$Lower.95[testCal$Measure=="Calibration slope"],3),
                 round(testCal$Lower.95[testCal$Measure=="Calibration-in-the-large"],3),
                 round(GNDresTest$pvalgw,3)),
    Upper.95 = c(round(quantile(R2DbootTest$R2,probs=0.975),3),
                 round(quantile(R2DbootTest$BS,probs=0.975),3),
                 round(TestC+qnorm(1-0.05/2)*TestCse,3),
                 round(quantile(R2DbootTest$D2,probs=0.975),3),
                 round(testCal$Upper.95[testCal$Measure=="Calibration slope"],3),
                 round(testCal$Upper.95[testCal$Measure=="Calibration-in-the-large"],3),
                 "NA")
  )
  datTestModel$Group <- "Test"
  if (!is.null(riskThreshold)) {
    datTestCal <- data.frame(
      Decile = c(rep(1:(length(riskThreshold)-1),2)),
      Risk = c(predTestYDecile,obsTestYDecile),
      Method=c(rep(c("Predicted","Observed"),each=length(riskThreshold)-1)))
  } else {
    datTestCal <- data.frame(
      Decile = c(rep(1:length(gNamesTest),2)),
      Risk = c(predTestYDecile,obsTestYDecile),
      Method=c(rep(c("Predicted","Observed"),each=length(gNamesTest))))
  }
  datTestCal$Group <- "Test"

  ## Results
  # Model performance
  datModelEval <- rbind(datTrainModel,datTestModel);row.names(datModelEval) <- NULL
  datModelCal <- rbind(datTrainCal,datTestCal); row.names(datModelCal) <- NULL
  datModelPred <- rbind(
    data.frame(time=dat[[time]],status=dat[[status]],predRisk=1-predTrainY,group="Training"),
    data.frame(time=newdat[[time]],status=newdat[[status]],predRisk=1-predTestY,group="Validation")
  )

  # Calibration plot
  p.Train <- ggplot(datModelCal[datModelCal$Group=="Train",],
                    aes(x=Decile,y=Risk,group=Method,fill=Method)) +
    geom_bar(stat="identity",position="dodge") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    scale_y_continuous(name="Probability of Event",expand=c(0,0), limits=c(0,1)) +
    scale_fill_manual("fill",values=c("Observed"="grey","Predicted"="black"))

  if(!is.null(riskThreshold)) {
    p.Train <- p.Train +
      scale_x_continuous(name="Risk group",
                         breaks=1:(length(riskThreshold)-1),
                         labels=gNames, expand=c(0,0),
                         limits=c(0.5,length(riskThreshold)-0.5))
  } else {
    p.Train <- p.Train +
      scale_x_continuous(name="Decile of Risk",breaks=seq(1,length(gNames),1),
                         labels=c(1:length(gNames)), expand=c(0,0),
                         limits=c(0.5,length(gNames)+0.5))
  }
  p.Train

  p.Test <- ggplot(datModelCal[datModelCal$Group=="Test",],
                   aes(x=Decile,y=Risk,group=Method,fill=Method)) +
    geom_bar(stat="identity",position="dodge") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    scale_y_continuous(name="Probability of Event",expand=c(0,0), limits=c(0,1)) +
    scale_fill_manual("fill",values=c("Observed"="grey","Predicted"="black"))
  if(!is.null(riskThreshold)) {
    p.Test <- p.Test +
      scale_x_continuous(name="Risk groups",
                         breaks=1:(length(riskThreshold)-1),
                         labels=gNames, expand=c(0,0),
                         limits=c(0.5,length(riskThreshold)-0.5))
  } else {
    p.Test <- p.Test +
      scale_x_continuous(name="Decile of Risk",breaks=seq(1,length(gNamesTest),1),
                         labels=c(1:length(gNamesTest)), expand=c(0,0),
                         limits=c(0.5,length(gNamesTest)+0.5))
  }
  p.Test

  return(list(dat.CIndex=dat.CIndex,
              dat.beta=dat.beta, dat.HR=dat.HR,
              opt.beta=opt.beta, opt.HR=opt.HR,
              opt.S0=opt.S0, opt.meanXB=opt.meanXB,
              dat.ModelPerformance=datModelEval,
              dat.ModelCalibration=datModelCal,
              dat.ModelPrediction=datModelPred,
              pCalibration.Train=p.Train,
              pCalibration.Test=p.Test))
}

#' @title Model evaluation using established risk prediction models
#' @param dat Data set used or the validation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Xs Candidate predictors
#' @param S0 The baseline survival probability
#' @param betaX The coefficients used in the risk prediction model
#' @param meanXB The means of candidate predictors
#' @param predictYear The t-year predicted risk
#' @param riskThreshold User-defined thresholds used for calculating calibration measures
#' @param calGroup User-defined number of groups used for calculating calibration measures
#' @param nboot The number of bootstrapped samples used for calculating the predictive performance
#' @export
coxModelEvaluation <- function(dat=datTrainZ,
                               time,status,Xs,
                               S0,betaX,meanXB,
                               predictYear,
                               riskThreshold=NULL,
                               calGroup,nboot=100) {
  dat <- data.table::setDT(dat)

  if( sum(!names(dat) %in% c(time,status,Xs))>0 ) stop("Please check the dataset!")

  curtail <- function(x) pmin(0.999999, pmax(x, 0.000001))

  opt.beta <- as.matrix(betaX,ncol=1); opt.beta
  # opt.HR <- dat.HR[dat.HR$fold==opt,]; opt.HR
  opt.S0 <- S0; opt.S0
  opt.meanXB <- meanXB; opt.meanXB
  # opt.Model <- modelFits[[opt]]$fitModel
  base.Model <- as.formula(paste("Surv(",time,",",status,")~1",sep=""))

  ## Model evaluation
  # (1) For training data set
  # Overall
  # BS <- Score(list("CoxModel"=opt.Model),
  #             data=dat,
  #             formula = base.Model,
  #             summary = "ipa", se.fit = TRUE,
  #             metrics = "brier", contrasts = FALSE,
  #             times = predictYear)
  # TrainBS <- BS$Brier$score[BS$Brier$score$model=="CoxModel","Brier"]; TrainBS
  # TrainBSse <- (TrainBS-BS$Brier$score[BS$Brier$score$model=="CoxModel","lower"])/1.96;TrainBSse
  # TrainR2 <- BS$Brier$score[BS$Brier$score$model=="CoxModel","IPA"]; TrainR2
  #
  # TrainDstat <- royston(opt.Model,newdata=dat)

  # predTrainY <- curtail(opt.S0^exp( t(t(dat[,Xs])) %*% t(opt.beta) - opt.meanXB )[,1])
  predTrainYLP <- (t(t(dat[,Xs,with=FALSE])) %*% t(opt.beta) - opt.meanXB)[,1]
  predTrainY <- curtail(opt.S0^exp(predTrainYLP))

  # C-index
  aucTrainY <- rcorr.cens(predTrainY, Surv(dat[[time]],dat[[status]]))
  TrainC <- aucTrainY['C Index']
  TrainCse <- aucTrainY['S.D.']/2
  # D-statistics & R2 & Brier score
  TrainR2D <- Rsq(lp=predTrainYLP,time=dat[[time]],status=dat[[status]],ties=TRUE)
  TrainBS <- brier_score(time=dat[[time]],
                         status=dat[[status]],
                         predictYear=predictYear,
                         survival=predTrainY)

  set.seed(123)
  requireNamespace(parallel)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {
    requireNamespace(survival)
    requireNamespace(dplyr)
    requireNamespace(riskRegression)
    requireNamespace(pec)
  })
  clusterExport(cl,list("dat","Rsq","brier_score","curtail",
                        "time","status","predictYear","Xs",
                        "opt.S0","opt.beta","opt.meanXB"),
                envir = environment())
  R2DBoots <- setNames(parLapply(
    cl,1:nboot,
    fun=function(i) {
      n <- nrow(dat)
      idx <- sample(1:n,size=n,replace=TRUE)
      datBoot <- dat[idx,]
      datBoot$ID <- 1:nrow(datBoot)
      predTrainYLPboot <- (t(t(datBoot[,Xs,with=FALSE])) %*% t(opt.beta) - opt.meanXB)[,1]
      predTrainYboot <- curtail(opt.S0^exp(predTrainYLPboot))
      TrainR2D <- Rsq(lp=predTrainYLPboot, time=datBoot[[time]], status=datBoot[[status]], ties = TRUE)
      TrainBS <- brier_score(time=datBoot[[time]], status=datBoot[[status]],
                             predictYear=predictYear, survival=predTrainYboot)
      data.frame(R2=TrainR2D[3],
                 D2=TrainR2D[1],
                 BS=TrainBS[1])
    }),1:nboot)
  stopCluster(cl)
  R2Dboot <- do.call("rbind",R2DBoots)

  # Calibration
  # Calibration plot
  if (!is.null(riskThreshold)) {
    TrainYDecileGrp <- cut(1-predTrainY,breaks=riskThreshold,right=FALSE)
  } else { TrainYDecileGrp <- Hmisc::cut2(1-predTrainY,g=calGroup) }
  GNDres <- GNDCalibrate(pred=1-predTrainY,
                         tvar=dat[[time]],out=dat[[status]],
                         cens.t=predictYear, groups=TrainYDecileGrp,
                         adm.cens=predictYear)
  gNames <- levels(TrainYDecileGrp);gNames
  gNames[length(gNames)] <- gsub(")","]",gNames[length(gNames)])
  predTrainYDecile <- GNDres$datTable$expectedperc; predTrainYDecile
  obsTrainYDecile <- GNDres$datTable$kmperc; obsTrainYDecile

  ## Calibration-in-the-large and calibration-slope
  trainCal <- getCal(dat=dat,time=time,status=status,
                     riskPred=1-predTrainY,predictYear=predictYear)
  datTrainModel <- data.frame(
    Measure = c("R2",
                "Brier Score",
                "C-statistic",
                "D-statistic",
                "Calibration Slope",
                "Calibration-in-the-large",
                "Greenwood-Nam-D'agostino Chisq (Lower.95=p-value)"),
    # data.frame(R2=TrainR2D[3],
    #            D2=TrainR2D[1],
    #            BS=TrainBS[1])
    Estimate = c(round(TrainR2D[3],3),
                 round(TrainBS[1],3),
                 round(TrainC,3),
                 round(TrainR2D[1],3),
                 round(trainCal$Estimate[trainCal$Measure=="Calibration slope"],3),
                 round(trainCal$Estimate[trainCal$Measure=="Calibration-in-the-large"],3),
                 round(GNDres$chi2gw,3)),
    Lower.95 = c(round(quantile(R2Dboot$R2,probs=0.025),3),
                 round(quantile(R2Dboot$BS,probs=0.025),3),
                 round(TrainC-qnorm(1-0.05/2)*TrainCse,3),
                 round(quantile(R2Dboot$D2,probs=0.025),3),
                 round(trainCal$Lower.95[trainCal$Measure=="Calibration slope"],3),
                 round(trainCal$Lower.95[trainCal$Measure=="Calibration-in-the-large"],3),
                 round(GNDres$pvalgw,3)),
    Upper.95 = c(round(quantile(R2Dboot$R2,probs=0.975),3),
                 round(quantile(R2Dboot$BS,probs=0.975),3),
                 round(TrainC+qnorm(1-0.05/2)*TrainCse,3),
                 round(quantile(R2Dboot$D2,probs=0.975),3),
                 round(trainCal$Upper.95[trainCal$Measure=="Calibration slope"],3),
                 round(trainCal$Upper.95[trainCal$Measure=="Calibration-in-the-large"],3),
                 "NA")
  )
  datTrainModel$Group <- "New"
  if (!is.null(riskThreshold)) {
    datTrainCal <- data.frame(
      Decile = c(rep(1:(length(riskThreshold)-1),2)),
      Risk = c(predTrainYDecile,obsTrainYDecile),
      Method=c(rep(c("Predicted","Observed"),each=(length(riskThreshold)-1))))
  } else {
    datTrainCal <- data.frame(
      Decile = c(rep(1:length(gNames),2)),
      Risk = c(predTrainYDecile,obsTrainYDecile),
      Method=c(rep(c("Predicted","Observed"),each=length(gNames))))
  }
  datTrainCal$Group <- "New"

  datModelEval <- rbind(datTrainModel);row.names(datModelEval) <- NULL
  datModelCal <- rbind(datTrainCal); row.names(datModelCal) <- NULL
  datModelPred <- rbind(
    data.frame(time=dat[[time]],status=dat[[status]],predRisk=1-predTrainY,group="New")
  )

  # Calibration plot
  p.Train <- ggplot(datTrainCal,
                    aes(x=Decile,y=Risk,group=Method,fill=Method)) +
    geom_bar(stat="identity",position="dodge") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    scale_y_continuous(name="Probability of Event",expand=c(0,0), limits=c(0,1)) +
    scale_fill_manual("fill",values=c("Observed"="grey","Predicted"="black"))

  if(!is.null(riskThreshold)) {
    p.Train <- p.Train +
      scale_x_continuous(name="Risk group",
                         breaks=1:(length(riskThreshold)-1),
                         labels=gNames, expand=c(0,0),
                         limits=c(0.5,length(riskThreshold)-0.5))
  } else {
    p.Train <- p.Train +
      scale_x_continuous(name="Decile of Risk",breaks=seq(1,length(gNames),1),
                         labels=c(1:length(gNames)), expand=c(0,0),
                         limits=c(0.5,length(gNames)+0.5))
  }
  p.Train

  return(list(#dat.CIndex=dat.CIndex,
    #dat.beta=dat.beta, dat.HR=dat.HR,
    opt.beta=opt.beta, #opt.HR=opt.HR,
    opt.S0=opt.S0, opt.meanXB=opt.meanXB,
    dat.ModelPerformance=datModelEval,
    dat.ModelCalibration=datModelCal,
    dat.ModelPrediction=datModelPred,
    # pCalibration.Train=p.Train,
    pCalibration=p.Train))

}


#' @title Develop risk prediction models with internal validation using bootstrap method
#' @param dat Data set or the derivation cohort
#' @param time Survival time
#' @param status Event status at the end of follow-up
#' @param Xs Candidate continuous predictors used for developing risk prediction models.
#' @param predictYear The horizon year or t-year predicted risk
#' @param KMS0 Whether calculate the baseline survival rate using the Kaplan-Meier estimator
#' @param riskThreshold User-defined thresholds used for calculating calibration measures
#' @param calGroup User-defined number of groups used for calculating calibration measures
#' @param nboot The number of bootstrapped samples used for calculating the predictive performance
#' @export
CoxModelBootstrapping <- function(dat,time,status,Xs,predictYear,KMS0=FALSE,
                                  riskThreshold,calGroup,nboot=100) {

  curtail <- function(x) pmin(0.999999, pmax(x, 0.000001))

  singleModel <- function(train,test,time,status,Xs,predictYear,KMS0,riskThreshold,calGroup) {
    # cal the mean of each predictor
    meanXs <- unlist(lapply(1:length(Xs), FUN=function(j) mean(train[[Xs[j]]]))); meanXs
    # construct the risk prediction model
    coxModel <- as.formula(paste("Surv(",time,",",status,")~",paste(Xs,collapse="+"),sep=""))
    coxFit <- cph(coxModel,data=train,x=TRUE,y=TRUE,surv=TRUE)
    # cal the baseline survival probability
    S0 <- survest(coxFit,times=predictYear)$surv;S0
    if(KMS0) {
      S0 <- summary(survfit(as.formula(paste("Surv(",time,",",status,")~1",sep="")),
                            data=train),times=c(predictYear))$surv; S0
    }
    dat.S0 <- data.frame(S0=S0)
    # cal the coefs
    betaXs <- matrix(coef(coxFit),ncol=1)
    dat.betaXs <- data.frame(t(betaXs))
    names(dat.betaXs) <- Xs
    # cal the MeanXB
    meanXB <- as.numeric(t(betaXs) %*% matrix(meanXs,ncol=1))
    dat.meanXB <- data.frame(meanXB=meanXB)
    datModel <- cbind(dat.S0,dat.meanXB,dat.betaXs)
    # cal the hazard for each risk factor
    dat.hrs <- data.frame(coef(coxFit),exp(cbind(coef(coxFit),confint(coxFit))))
    names(dat.hrs) <- c("betaXs","HR","LCI","UCI")
    dat.hrs$Xs <- Xs; rownames(dat.hrs) <- NULL
    datHR <- dat.hrs

    ## Part 1. For the training data set
    # cal the C-index
    trainXs <- train[,Xs]
    predtrainYLP <- (t(t(trainXs)) %*% betaXs - meanXB)
    predtrainY <- curtail(S0^exp(predtrainYLP)[,1])
    auc.train <- rcorr.cens(predtrainY, Surv(train[[time]],train[[status]]))
    Cindex.train <- auc.train['C Index']
    # cal the R2, BS, and D-index
    trainR2D <- Rsq(lp=predtrainYLP, time=train[[time]], status=train[[status]], ties = TRUE)
    trainBS <- brier_score(time=train[[time]], status=train[[status]],
                           predictYear=predictYear, survival=predtrainY)
    R2.train <- trainR2D[3]
    Dindex.train <- trainR2D[1]
    BS.train <- trainBS[1]
    # cal the calibration-related measures
    if (!is.null(riskThreshold)) {
      trainYDecileGrp <- cut(1-predtrainY,breaks=riskThreshold,right=FALSE)
    } else { trainYDecileGrp <- Hmisc::cut2(1-predtrainY,g=calGroup) }
    gNamestrain <- levels(trainYDecileGrp)
    gNamestrain[length(gNamestrain)] <- gsub(")","]",gNamestrain[length(gNamestrain)])
    GNDrestrain <- GNDCalibrate(pred=1-predtrainY,
                                tvar=train[[time]],out=train[[status]],
                                cens.t=predictYear, groups=trainYDecileGrp,
                                adm.cens=predictYear)
    predtrainYDecile <- GNDrestrain$datTable$expectedperc; predtrainYDecile
    obstrainYDecile <- GNDrestrain$datTable$kmperc; obstrainYDecile
    GNDChisq.train <- GNDrestrain$chi2gw
    GNDPval.train <- GNDrestrain$pvalgw
    ## Calibration statistics
    trainCal <- getCal(dat=train,time=time,status=status,
                       riskPred=1-predtrainY,predictYear=predictYear)
    calLarge.train <- trainCal$Estimate[1]
    calSlope.train <- trainCal$Estimate[2]

    ## Part.2 For the validation data set
    # cal the C-index
    testXs <- test[,Xs]
    predtestYLP <- (t(t(testXs)) %*% betaXs - meanXB)
    predtestY <- curtail(S0^exp(predtestYLP)[,1])
    auc.test <- rcorr.cens(predtestY, Surv(test[[time]],test[[status]]))
    Cindex.test <- auc.test['C Index']
    # cal the R2, BS, and D-index
    testR2D <- Rsq(lp=predtestYLP, time=test[[time]], status=test[[status]], ties = TRUE)
    testBS <- brier_score(time=test[[time]], status=test[[status]],
                          predictYear=predictYear, survival=predtestY)
    R2.test <- testR2D[3]
    Dindex.test <- testR2D[1]
    BS.test <- testBS[1]
    # cal the calibration-related measures
    if (!is.null(riskThreshold)) {
      testYDecileGrp <- cut(1-predtestY,breaks=riskThreshold,right=FALSE)
    } else { testYDecileGrp <- Hmisc::cut2(1-predtestY,g=calGroup) }
    gNamestest <- levels(testYDecileGrp)
    gNamestest[length(gNamestest)] <- gsub(")","]",gNamestest[length(gNamestest)])
    GNDrestest <- GNDCalibrate(pred=1-predtestY,
                               tvar=test[[time]],out=test[[status]],
                               cens.t=predictYear, groups=testYDecileGrp,
                               adm.cens=predictYear)
    predtestYDecile <- GNDrestest$datTable$expectedperc; predtestYDecile
    obstestYDecile <- GNDrestest$datTable$kmperc; obstestYDecile
    GNDChisq.test <- GNDrestest$chi2gw
    GNDPval.test <- GNDrestest$pvalgw
    ## Calibration statistics
    testCal <- getCal(dat=test,time=time,status=status,
                      riskPred=1-predtestY,predictYear=predictYear)
    calLarge.test <- testCal$Estimate[1]
    calSlope.test <- testCal$Estimate[2]
    datModelPerformance <- data.frame(
      R2.train=R2.train,R2.test=R2.test,
      BS.train=BS.train,BS.test=BS.test,
      Cindex.train=Cindex.train,Cindex.test=Cindex.test,
      Dindex.train=Dindex.train,Dindex.test=Dindex.test,
      calSlope.train=calSlope.train,calSlope.test=calSlope.test,
      calLarge.train=calLarge.train,calLarge.test=calLarge.test,
      GNDChisq.train=GNDChisq.train,GNDChisq.test=GNDChisq.test,
      GNDPval.train=GNDPval.train,GNDPval.test=GNDPval.test,
      GNDdf=GNDrestest$df
    )
    datModelPredict <- data.frame(riskGroup=gNamestrain,
                                  predY.train=predtrainYDecile,
                                  obsY.train=obstrainYDecile,
                                  predY.test=predtestYDecile,
                                  obsY.test=obstestYDecile)
    list(datModel=datModel,datModelHR=datHR,
         datModelPerformance=datModelPerformance,
         datModelPredict=datModelPredict)
  }

  set.seed(123)
  requireNamespace(parallel)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {
    requireNamespace(survival)
    requireNamespace(dplyr)
    requireNamespace(riskRegression)
    requireNamespace(pec)
    requireNamespace(rms)
  })
  clusterExport(cl,list("dat","Rsq","brier_score","curtail","getCal",
                        "singleModel","GNDCalibrate",
                        "time","status","Xs","predictYear",
                        "KMS0","riskThreshold","calGroup","nboot"),
                envir = environment())
  datBootstrapping <- setNames(parLapply(
    cl,1:nboot,
    fun=function(i) {
      n <- nrow(dat)
      idx <- sample(1:n,size=n,replace=TRUE)
      datBoot <- dat[idx,]
      datBoot$ID <- 1:nrow(datBoot)
      datMP <- tryCatch({
        singleModel(train=as.data.frame(datBoot),
                    test=as.data.frame(dat),
                    time=time,status=status,Xs=Xs,
                    predictYear=predictYear,KMS0=KMS0,
                    riskThreshold=riskThreshold,
                    calGroup=calGroup)
      }, error=function(e) {
        NULL
      })
      if (!is.null(datMP)) {
        datMP$datModel$nboot <- i
        datMP$datModelHR$nboot <- i
        datMP$datModelPerformance$nboot <- i
        datMP$datModelPredict$nboot <- i
      }
      datMP
    }),1:nboot)
  stopCluster(cl)

  datModel <- do.call(
    "rbind",
    lapply(1:nboot, FUN=function(i) datBootstrapping[[i]]$datModel))
  datModelHR <- do.call(
    "rbind",
    lapply(1:nboot, FUN=function(i) datBootstrapping[[i]]$datModelHR))
  datMP <- do.call(
    "rbind",
    lapply(1:nboot, FUN=function(i) datBootstrapping[[i]]$datModelPerformance))
  datModelPredict <- do.call(
    "rbind",
    lapply(1:nboot, FUN=function(i) datBootstrapping[[i]]$datModelPredict))
  datMPWeighted <- data.frame(
    R2 = datMP$R2.train - with(datMP, mean(R2.train - R2.test, na.rm=TRUE)),
    BS = datMP$BS.train - with(datMP, mean(BS.train - BS.test, na.rm=TRUE)),
    Cindex = datMP$Cindex.train - with(datMP, mean(Cindex.train - Cindex.test, na.rm=TRUE)),
    Dindex = datMP$Dindex.train - with(datMP, mean(Dindex.train - Dindex.train, na.rm=TRUE)),
    calSlope = datMP$calSlope.train - with(datMP, mean(calSlope.train - calSlope.test, na.rm=TRUE)),
    calLarge = datMP$calLarge.train - with(datMP, mean(calLarge.train - calLarge.test, na.rm=TRUE)),
    # GNDChisq = datMP$GNDChisq.train - with(datMP, mean(GNDChisq.train - GNDChisq.test, na.rm=TRUE)),
    GNDChisq = datMP$GNDChisq.train,
    GNDdf=datMP$GNDdf)
  datMPWeighted$GNDPval <- 1-with(datMPWeighted,pchisq(GNDChisq,GNDdf))
  d.f <- datMPWeighted$GNDdf[1]
  idx.boot <- datModel$nboot
  best.model <- which.max(with(datMPWeighted,Cindex))
  worst.model <- which.min(with(datMPWeighted,Cindex))
  n.col <- ncol(datMPWeighted)
  datModelPerformance <- data.frame(
    Measure = c("R2",
                "Brier Score",
                "C-statistic",
                "D-statistic",
                "Calibration Slope",
                "Calibration-in-the-large",
                "Greenwood-Nam-D'agostino Chisq",
                "Greenwood-Nam-D'agostino df"),
    Estimate = c(apply(datMPWeighted[,1:(n.col-1)],2,FUN=function(x) round(quantile(x,probs=0.5,na.rm=TRUE),3))),
    Lower.95 = c(apply(datMPWeighted[,1:(n.col-1)],2,FUN=function(x) round(quantile(x,probs=0.025,na.rm=TRUE),3))),
    Upper.95 = c(apply(datMPWeighted[,1:(n.col-1)],2,FUN=function(x) round(quantile(x,probs=0.975,na.rm=TRUE),3)))
  )
  datMPpval <- data.frame(
    Measure = c("Greenwood-Nam-D'agostino P-value"),
    Estimate = round(1-with(datModelPerformance,pchisq(Estimate[grepl("Chisq",Measure)],df=d.f)),3),
    Lower.95 = round(1-with(datModelPerformance,pchisq(Lower.95[grepl("Chisq",Measure)],df=d.f)),3),
    Upper.95 = round(1-with(datModelPerformance,pchisq(Upper.95[grepl("Chisq",Measure)],df=d.f)),3))
  datModelPerformance <- rbind(datModelPerformance,datMPpval)
  row.names(datModelPerformance) <- NULL
  res <- list()
  ## Overall performance
  res$bootstrapping <- datModelPerformance
  ## Best model
  res$bestModel <- list(
    datModel = datModel[datModel$nboot==idx.boot[best.model],],
    datModelHR = datModelHR[datModelHR$nboot==idx.boot[best.model],],
    datModelPerformance = datMPWeighted[best.model,],
    datModelPredict = datModelPredict[datModelPredict$nboot==idx.boot[best.model],]
  )
  ## Worst model
  res$worstModel <- list(
    datModel = datModel[datModel$nboot==idx.boot[worst.model],],
    datModelHR = datModelHR[datModelHR$nboot==idx.boot[worst.model],],
    datModelPerformance = datMPWeighted[worst.model,],
    datModelPredict = datModelPredict[datModelPredict$nboot==idx.boot[worst.model],]
  )
  res
}
