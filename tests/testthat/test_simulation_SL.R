##############################
# Test main medltmle functions
##############################
library(here)
library(SuperLearner)
library(matrixStats)
library(parallel)
library(speedglm)
library(Matrix)
library(pracma)
library(reshape)
library(caret)
library(hal)

context("Overall Test for medltmle with simple SuperLearner")

#Load all scripts
setwd(here("R"))
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#Set seed:
set.seed(2)

#Generate some simulated data:
data<-GenerateData(n=400, end.time=2, abar=NULL,abar.prime=NULL)

#Generate appropriate models:
spec<-make.sim.spec(2)

#Some parameters:
end.time=2

#######################
# Super Learner
#######################

#Wrappers to tune algorithms with caret package
m=5

conf = list(
  # Number of hyperparameter configurations to try per model
  tune_length = 5,
  caret_parallel = T,
  # Internal CV folds.
  caret_cv_folds = 3,
  # Improve stability of CV.
  caret_cv_repeats = 3,
  caret_metric = "ROC"
)

# General train control settings to re-use in each caret wrapper.
caret_train_control = trainControl(method = "repeatedcv", verboseIter = F,
                                   number = conf$caret_cv_folds,
                                   # Add repeats to improve stability.
                                   repeats = conf$caret_cv_repeats,
                                   allowParallel = conf$caret_parallel,
                                   summaryFunction = twoClassSummary,
                                   classProbs = T)


SL.glmnet.caret1 <- function(...,method="glmnet",tuneLength=m, trControl=trainControl(method="cv",number=10,verboseIter=FALSE)){
  SL.caret(...,method=method,tuneLength=tuneLength,trControl=trControl)
}


#Simple SuperLearner libraries:
lib <- list(Q=c("SL.glm", "SL.mean","SL.step", "SL.nnet"),g=c("SL.glm","SL.mean","SL.step","SL.glmnet.caret1"))

result_10 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      Dnodes=NULL,
                      W2nodes=NULL,
                      survivalOutcome = T,
                      QLform=spec$QL.c,
                      QZform=spec$QZ.c,
                      gform=spec$g.c,
                      qzform=spec$qz.c,
                      qLform=spec$qL.c,
                      abar=rep(1,end.time),
                      abar.prime=rep(0,end.time),
                      gbounds=c(.01,.99),
                      deterministic.g.function = NULL,
                      stratify=FALSE,
                      SL.library=lib,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      estimand="NE",
                      time.end=end.time
)

result_00 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      Dnodes=NULL,
                      W2nodes=NULL,
                      survivalOutcome = T,
                      QLform=spec$QL.c,
                      QZform=spec$QZ.c,
                      gform=spec$g.c,
                      qzform=spec$qz.c,
                      qLform=spec$qL.c,
                      abar=rep(0,end.time),
                      abar.prime=rep(0,end.time),
                      gbounds=c(.01,.99),
                      deterministic.g.function = NULL,
                      stratify=FALSE,
                      SL.library=lib,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      estimand="NE",
                      time.end=end.time
)

result_11 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      Dnodes=NULL,
                      W2nodes=NULL,
                      survivalOutcome = T,
                      QLform=spec$QL.c,
                      QZform=spec$QZ.c,
                      gform=spec$g.c,
                      qzform=spec$qz.c,
                      qLform=spec$qL.c,
                      abar=rep(1,end.time),
                      abar.prime=rep(1,end.time),
                      gbounds=c(.01,.99),
                      deterministic.g.function = NULL,
                      stratify=FALSE,
                      SL.library=lib,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      estimand="NE",
                      time.end=end.time
)

result_01 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      Dnodes=NULL,
                      W2nodes=NULL,
                      survivalOutcome = T,
                      QLform=spec$QL.c,
                      QZform=spec$QZ.c,
                      gform=spec$g.c,
                      qzform=spec$qz.c,
                      qLform=spec$qL.c,
                      abar=rep(0,end.time),
                      abar.prime=rep(1,end.time),
                      gbounds=c(.01,.99),
                      deterministic.g.function = NULL,
                      stratify=FALSE,
                      SL.library=lib,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      estimand="NE",
                      time.end=end.time
)

#Natural Indirect Effect:
NIE<-result_11$estimates[1]-result_10$estimates[1]

#Natural Direct Effect:
NDE<-result_10$estimates[1]-result_00$estimates[1]

#Overall Natural Effect:
NE<-NIE+NDE

#Test TMLE NIE
test_that("TMLE estimate of NIE for the simulation 1 with SL matches expected", expect_equal(NIE[[1]], -0.003045925, tolerance = 0.01))

#Test TMLE NDE
test_that("TMLE estimate of NDE for the simulation 1 with SL matches expected", expect_equal(NDE[[1]], 0.04662457, tolerance = 0.01))
