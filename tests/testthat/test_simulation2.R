#########################################################
# Test main medltmle functions on a different simulation
#########################################################
library(here)
library(SuperLearner)
library(matrixStats)
library(parallel)
library(speedglm)
library(Matrix)
library(pracma)

context("Overall Test for medltmle for simulation with a single exposure")

#Load all scripts
#setwd(here("R"))
#file.sources = list.files(pattern="*.R")
#sapply(file.sources,source,.GlobalEnv)

#Load all scripts to generate data
#setwd(here("simulation"))
#file.sources = list.files(pattern="*.R")
#sapply(file.sources,source,.GlobalEnv)

#Set seed:
set.seed(2)

#Generate some simulated data:
data<-GenerateData_SingA_TimeOrdL(n=400, end.time=2, abar=NULL, abar.prime=NULL)

#Generate appropriate models:
spec<-make.sim.spec(2)

#Some parameters:
end.time=2
abar <- 1
abar.prime <- 0

result.c <- medltmle(data=data,
                     Anodes=names(data)[grep('^A',names(data))],
                     Cnodes=names(data)[grep('^C',names(data))],
                     Znodes=names(data)[grep('^Z',names(data))],
                     Lnodes=names(data)[grep('^L',names(data))],
                     Ynodes=names(data)[grep('^Y',names(data))],
                     survivalOutcome = T,
                     QLform=spec$QL.c,
                     QZform=spec$QZ.c,
                     gform=spec$g.c,
                     qzform=spec$qz.c,
                     qLform=spec$qL.c,
                     abar=rep(abar,end.time),
                     abar.prime=rep(abar.prime,end.time),
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
                     observation.weights=NULL
)

#Test TMLE
test_that("TMLE estimate for the simulation matches expected", expect_equal(result.c$estimates[1],
                                                                            0.9611661, tolerance = 0.01))

#Test IPW
test_that("IPW estimate for the simulation matches expected", expect_equal(result.c$estimates[2],
                                                                           0.9566909, tolerance = 0.01))
