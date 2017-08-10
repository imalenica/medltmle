#########################################################
# Test main medltmle functions on a different simulation
#########################################################
#library(here)
#library(SuperLearner)
#library(matrixStats)
#library(parallel)
#library(speedglm)
#library(Matrix)
#library(pracma)
#library(stremr)

context("Overall Test for medltmle for simulation with a single exposure")
#More specifically, we want to test the code when there is a single exposure and outcome, and baseline is time-dependent.
#W2 will need to be fluctuated as well.

#Load all scripts
#setwd(here("R"))
#file.sources = list.files(pattern="*.R")
#sapply(file.sources,source,.GlobalEnv)

#Set seed:
set.seed(2)
end.time=2

#Generate some simulated data:
data<-GenerateData_SingA_TimeOrdL(n=400, end.time=end.time, abar=NULL, abar.prime=NULL)

#Turn LZ into binary (if smaller that 1st quartile)
data$LZ_1<-ifelse(data$LZ_1<11,1,0)
data$LZ_2<-ifelse(data$LZ_2<11.466,1,0)

#Some parameters:
end.time=2

result_11 <- medltmle(data=data,
                     Anodes=names(data)[grep('^A',names(data))],
                     Cnodes=names(data)[grep('^C',names(data))],
                     Dnodes=names(data)[grep('^D',names(data))],
                     Znodes=names(data)[grep('^Z',names(data))],
                     Lnodes=names(data)[grep('^L',names(data))],
                     Ynodes=names(data)[grep('^Y',names(data))],
                     W2nodes=names(data)[grep('W2',names(data))],
                     survivalOutcome = F,
                     QLform=NULL,
                     QZform=NULL,
                     gform=NULL,
                     qzform=c("Z_1~LA_1+A+B6.W2+B5.2.W2+B5.1.W2+B3.1.W1", "Z_2~LA_2+A"),
                     qLform=c("LA_1~A+B6.W2+B5.2.W2+B5.1.W2","LZ_1~A+LA_1+Z_1","Y_1~Z_1+LA_1+A+LZ_1","LA_2~LA_1+A+Z_2","LZ_2~LA_2+A+Z_2","Y_2~Z_2+LA_2+A+LZ_2"),
                     abar=1,
                     abar.prime=1,
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
                     IC.variance.only=TRUE,
                     observation.weights=NULL,
                     CSE=TRUE,
                     time.end=end.time
)

result_10 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Dnodes=names(data)[grep('^D',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      W2nodes=names(data)[grep('W2',names(data))],
                      survivalOutcome = F,
                      QLform=NULL,
                      QZform=NULL,
                      gform=NULL,
                      qzform=c("Z_1~LA_1+A+B6.W2+B5.2.W2+B5.1.W2+B3.1.W1", "Z_2~LA_2+A"),
                      qLform=c("LA_1~A+B6.W2+B5.2.W2+B5.1.W2","LZ_1~A+LA_1+Z_1","Y_1~Z_1+LA_1+A+LZ_1","LA_2~LA_1+A+Z_2","LZ_2~LA_2+A+Z_2","Y_2~Z_2+LA_2+A+LZ_2"),
                      abar=1,
                      abar.prime=0,
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
                      IC.variance.only=TRUE,
                      observation.weights=NULL,
                      CSE=TRUE,
                      time.end=end.time
)

result_01 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Dnodes=names(data)[grep('^D',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      W2nodes=names(data)[grep('W2',names(data))],
                      survivalOutcome = F,
                      QLform=NULL,
                      QZform=NULL,
                      gform=NULL,
                      qzform=c("Z_1~LA_1+A+B6.W2+B5.2.W2+B5.1.W2+B3.1.W1", "Z_2~LA_2+A"),
                      qLform=c("LA_1~A+B6.W2+B5.2.W2+B5.1.W2","LZ_1~A+LA_1+Z_1","Y_1~Z_1+LA_1+A+LZ_1","LA_2~LA_1+A+Z_2","LZ_2~LA_2+A+Z_2","Y_2~Z_2+LA_2+A+LZ_2"),
                      abar=0,
                      abar.prime=1,
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
                      IC.variance.only=TRUE,
                      observation.weights=NULL,
                      CSE=TRUE,
                      time.end=end.time
)

result_00 <- medltmle(data=data,
                      Anodes=names(data)[grep('^A',names(data))],
                      Cnodes=names(data)[grep('^C',names(data))],
                      Dnodes=names(data)[grep('^D',names(data))],
                      Znodes=names(data)[grep('^Z',names(data))],
                      Lnodes=names(data)[grep('^L',names(data))],
                      Ynodes=names(data)[grep('^Y',names(data))],
                      W2nodes=names(data)[grep('W2',names(data))],
                      survivalOutcome = F,
                      QLform=NULL,
                      QZform=NULL,
                      gform=NULL,
                      qzform=c("Z_1~LA_1+A+B6.W2+B5.2.W2+B5.1.W2+B3.1.W1", "Z_2~LA_2+A"),
                      qLform=c("LA_1~A+B6.W2+B5.2.W2+B5.1.W2","LZ_1~A+LA_1+Z_1","Y_1~Z_1+LA_1+A+LZ_1","LA_2~LA_1+A+Z_2","LZ_2~LA_2+A+Z_2","Y_2~Z_2+LA_2+A+LZ_2"),
                      abar=0,
                      abar.prime=0,
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
                      IC.variance.only=TRUE,
                      observation.weights=NULL,
                      CSE=TRUE,
                      time.end=end.time
)

#Natural Indirect Effect:
NIE<-result_11$estimates[1]-result_10$estimates[1]

#Natural Direct Effect:
NDE<-result_10$estimates[1]-result_00$estimates[1]

#Overall Natural Effect:
NE<-NIE+NDE

#Test TMLE NIE
test_that("TMLE estimate of NIE for the simulation 2 matches expected", expect_equal(NIE[[1]],  0.00291363, tolerance = 0.01))

#Test TMLE NDE
test_that("TMLE estimate of NDE for the simulation 2 matches expected", expect_equal(NDE[[1]], -0.0301333, tolerance = 0.01))

#Check summary_medltmle function:
res<-summary_medltmle(nie1=result_11,nie2=result_10,nde1=result_10,nde2=result_00)

res_NDE<-res$NDE
res_NIE<-res$NIE
res_NE<-res$NE

#Test TMLE NIE SE
test_that("TMLE SE of NIE for the simulation 1 matches expected", expect_equal(res_NIE[1,2], 0.01146235, tolerance = 0.01))

#Test TMLE NDE SE
test_that("TMLE SE of NDE for the simulation 1 matches expected", expect_equal(res_NDE[1,2], 0.0557253, tolerance = 0.01))


