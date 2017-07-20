#' GenerateData_SingA_TimeOrdL
#'
#' Simulate another longitudinal data set, with uncertain time ordering of baseline covariates and single intervention.
#' Note that the data is of the following format: O={L1,L2,L3,L4,L5,L6,C1,D1,C1.2,LA1,Z1,LZ1,C2,D2,LA2,Z2,LZ2 ... Y }
#' L1,L2,L3,L4,L5 are baseline covariates accroding to time ordering, L6 is some randomized covariate (possibly intervention)
#' LA and LZ time-varying covariates, C and D are censoring and death indicators, Z is a mediator and Y is the final outcome based on LZ.
#'
#' @param n Sample size
#' @param end.time Last time point
#' @param abar Set intervention
#' @param abar.prime Control intervention
#'
#' @return Returns a dataframe with simulated covariates as described in the description.
#'
#' @export GenerateData_SingA_TimeOrdL

GenerateData_SingA_TimeOrdL <- function(n, end.time, abar=NULL,abar.prime=NULL) {

  if(is.null(abar) != is.null(abar.prime)){stop('abar and abar.prime either both given or both not given')}
  if((!is.null(abar) & length(abar)!=end.time) | (!is.null(abar.prime) & length(abar.prime)!=end.time)){stop('abar and abar.prime both need to be length of end.time')}

  #Define some functions:
  rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

  ###########################
  #Data generating process:
  ###########################

  #Note, baseline covariates are now not independent.

  #Belong to the first group of baseline covariates:
  B1.1 <- rbinom(n,1,0.4)
  B1.2 <- rbinom(n,1,0.6)

  #Belong to the second group of baseline covariates:
  B2.1 <- rexpit(1.2-0.7*B1.1)
  B2.2 <- rexpit(0.4-0.7*B1.2)

  #Belong to the third group of baseline covariates:
  B3.1 <- rexpit(1.2-0.7*B1.1+0.1*B2.2)

  #Belong to the fourth group of baseline covariates:
  B4.1 <- rexpit(0.3-0.3*B1.2+0.4*B3.1+0.5*B2.1)

  #Belong to the fifth group of baseline covariates:
  B5.1 <- rexpit(4.1+0.3*B2.2+B4.1-4.5*B1.1+0.2*B3.1)
  B5.2 <- rexpit(2.1-0.3*B2.2+0.3*B4.1-2*B1.1+2*B3.1)

  #Possibly an exposure from a randomized trial.
  B6 <- rbinom(n,1,.5)

  #Group all covariates.
  covs<-cbind.data.frame(B1.1,B1.2,B2.1,B2.2,B3.1,B4.1,B5.1,B5.2,B6)

  #Set A and order the baseline covariates.
  covs_order<-timeOrder_baseline(data=covs,A=6)

  LA <- LZ <- C <- Z <- D <- Y  <- matrix(NA, nrow=n, ncol=end.time)

  uncensored.alive <- rep(TRUE, n)

  est.psi0 <- !is.null(abar) & !is.null(abar.prime)

  #Function for Z:
  CalcZ <- function(W,LA_n,A) {
    lin <- -0.5+0.8*W+0.8*A + 1*LA_n
    Z <- rexpit(lin)
    return(Z)
  }

  for (t in 1:(end.time)) {

    ### C
    if(est.psi0){
      C[uncensored.alive,t] <- 1
    }else{
      if(t==1) C[uncensored.alive,t] <- rexpit(4-0.2*covs_order$B1.1.W1[uncensored.alive]-0.3*covs_order$B2.1.W1[uncensored.alive])
      else C[uncensored.alive,t] <-rexpit(4-0.2*covs_order$B1.1.W1[uncensored.alive]+0.2*covs_order$A[uncensored.alive]-0.01*LZ[uncensored.alive,t-1])
    }

    ### D
    if(est.psi0){
      D[uncensored.alive,t] <- 1
    }else{
      if(t==1) D[uncensored.alive,t] <- 1-rexpit(4-0.2*covs_order$B1.2.W1[uncensored.alive]-0.3*covs_order$B3.1.W1[uncensored.alive])
      else D[uncensored.alive,t] <- 1-rexpit(4-0.2*covs_order$B1.2.W1[uncensored.alive]+0.1*covs_order$A[uncensored.alive]-0.08*LZ[uncensored.alive,t-1])
    }

    #Update who died (note more might have died from outcome Y1 to now):
    uncensored.alive <- uncensored.alive & C[,t] & !D[,t]

    ## LA
    if(t==1){
      LA[uncensored.alive,t] <- rexpit(-.8+.1*covs_order$B3.1.W1[uncensored.alive]+.3*covs_order$B5.2.W2[uncensored.alive]+covs_order$A[uncensored.alive])
    }else{
      LA[uncensored.alive, t] <- rexpit(-.8+.1*LZ[uncensored.alive,t-1]+.3*LA[uncensored.alive,t-1]+covs_order$A[uncensored.alive])
    }

    ## Z
    if (est.psi0) {
      useA <- abar.prime[t]
    }else{
      useA  <- covs_order$A[uncensored.alive]
    }

    Z[uncensored.alive, t] <- CalcZ(W=covs_order$B3.1.W1[uncensored.alive],LA_n=LA[uncensored.alive,t],A = useA)

    ## LZ: should be continous.
    if(t==1) {
      LZ[uncensored.alive,t] <- 10+3*covs_order$B2.2.W1[uncensored.alive]+covs_order$A[uncensored.alive]+0.7*Z[uncensored.alive,t]
    }else {
      LZ[uncensored.alive,t] <- 10+3*covs_order$B2.2.W1[uncensored.alive]+covs_order$A[uncensored.alive]+0.7*Z[uncensored.alive,t] -0.02*LZ[uncensored.alive,t-1]
    }

    Y[uncensored.alive, t]<-ifelse(LZ[uncensored.alive,t]<12,1,0)

    ## Y: deterministic function of LZ.
    #if(t==end.time){

      #Y[uncensored.alive]<-ifelse(LZ[uncensored.alive,t]<8.5,1,0)

    #}

    }

  .BinaryToCensoring <- function(is.censored, is.uncensored) {
    if (! xor(missing(is.censored), missing(is.uncensored))) stop("exactly one of is.censored and is.uncensored must be passed")
    calling.name <- names(sys.call(0))[2]
    if (length(calling.name) == 0 || ! calling.name %in% c("is.censored", "is.uncensored")) stop("the argument to BinaryToCensoring must be completely named - see ?BinaryToCensoring")
    if (missing(is.uncensored)) {
      is.uncensored <- ! is.censored
    }
    if (! all(is.uncensored %in% c(0, 1, NA))) stop("the argument to BinaryToCensoring should be binary (0, 1, or NA) or logical")
    y <- character(length(is.uncensored))
    y[is.uncensored == 0] <- "censored"
    y[is.uncensored == 1] <- "uncensored"
    y[is.na(is.uncensored)] <- NA
    return(factor(y))
  }

  .CreateDataFrame <- function(baseline, C, D, Z, LA, LZ, Y, end.time) {

    d=baseline

    for (t in 1:(end.time)) {

      d <- data.frame(d, .BinaryToCensoring(is.uncensored=C[,t]),D[,t],LA[,t],Z[,t],LZ[,t],Y[,t])
      names(d)[ncol(d) - 5] <- paste0("C_", t)
      names(d)[ncol(d) - 4] <- paste0("D_", t)
      names(d)[ncol(d) - 3] <- paste0("LA_", t)
      names(d)[ncol(d) - 2] <- paste0("Z_", t)
      names(d)[ncol(d) - 1] <- paste0("LZ_", t)
      names(d)[ncol(d)] <- paste0("Y_", t)
    }

      #if(t!=end.time){
       # d <- data.frame(d, .BinaryToCensoring(is.uncensored=C[, t]), D[, t], LA[,t],Z[,t],LZ[,t])
        #names(d)[ncol(d) - 4] <- paste0("C_", t)
        #names(d)[ncol(d) - 3] <- paste0("D_", t)
        #names(d)[ncol(d) - 2] <- paste0("LA_", t)
        #names(d)[ncol(d) - 1] <- paste0("Z_", t)
        #names(d)[ncol(d)] <- paste0("LZ_", t)
      #}else{

        #d <- data.frame(d, .BinaryToCensoring(is.uncensored=C[, t]), D[, t], LA[,t],Z[,t],LZ[,t],Y)
        #names(d)[ncol(d) - 5] <- paste0("C_", t)
        #names(d)[ncol(d) - 4] <- paste0("D_", t)
        #names(d)[ncol(d) - 3] <- paste0("LA_", t)
        #names(d)[ncol(d) - 2] <- paste0("Z_", t)
        #names(d)[ncol(d) - 1] <- paste0("LZ_", t)
        #names(d)[ncol(d)] <- paste0("Y_", t)
      #}
  #}
    return(d)
  }

  return(.CreateDataFrame(baseline=covs_order, C = C,D = D,Z = Z,LA = LA,LZ = LZ, Y =  Y,end.time =  end.time))

}
