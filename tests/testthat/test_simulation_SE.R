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

context("Overall Test for medltmle: Stochastic Effect")

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
                      SL.library=NULL,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      CSE=FALSE,
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
                      SL.library=NULL,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      CSE=FALSE,
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
                      SL.library=NULL,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      CSE=FALSE,
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
                      SL.library=NULL,
                      estimate.time=FALSE,
                      deterministic.Q.function=NULL,
                      rule=NULL,
                      Yrange=NULL,
                      gcomp=FALSE,
                      iptw.only=FALSE,
                      IC.variance.only=FALSE,
                      observation.weights=NULL,
                      CSE=FALSE,
                      time.end=end.time
)

#Natural Indirect Effect:
SIE<-result_11$estimates[1]-result_10$estimates[1]

#Natural Direct Effect:
SDE<-result_10$estimates[1]-result_00$estimates[1]

#Overall Natural Effect:
SE<-SIE+SDE

#Test TMLE NIE
#test_that("TMLE estimate of SIE for the simulation 1 matches expected", expect_equal(SIE[[1]], -0.02720144, tolerance = 0.01))

#Test TMLE NDE
#test_that("TMLE estimate of SDE for the simulation 1 matches expected", expect_equal(SDE[[1]], 0.04102935, tolerance = 0.01))

