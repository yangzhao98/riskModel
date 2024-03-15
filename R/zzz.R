.onLoad <- function(libname=find.package("drugTargetScreen"),pkgname="drugTargetScreen"){

  # CRAN Note avoidance
  if (getRversion() >= "3.1.0")
    utils::globalVariables(c("%>%", ".", "stats", "as.formula", "coef", "lm", "pchisq", "qnorm",
                             "quantile", "vcov", "Surv", "case_when",
                             "coxph", "left_join mutate", "predictSurvProb",
                             "predmort","predsurv","predsurv_null",
                             "select", "surv", "survfit", "weight",
                             "Decile", "Method", "Risk", "aes",
                             "clusterEvalQ", "clusterExport", "confint", "cph",
                             "detectCores", "element_blank", "geom_bar", "ggplot",
                             "improveProb", "makeCluster", "parLapply", "pnorm",
                             "rcorr.cens", "scale_fill_manual", "scale_x_continuous",
                             "scale_y_continuous", "sd", "setNames", "stopCluster",
                             "survest", "theme", "theme_classic","stats",
                             "confint", "pnorm", "sd", "setNames",
                             "ID", "Risk_groups", "Score", "abline", "axis",
                             "capture.output", "complete.cases",
                             "datTrainZ", "g", "gaussian", "group", "group_by",
                             "labs", "legend", "lines", "lowess", "par",
                             "plotCalibration", "polygon", "predict",
                             "qt", "risk", "rmsPredict", "segments",
                             "smooth_pseudo", "summarise","dplyr",
                             "parallel", "pec", "riskRegression", "rms", "survival"))

  invisible()

}
