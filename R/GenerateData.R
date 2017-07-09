#' GenerateData
#'
#' Simulate longitudinal data as specified in the description. Note that the data is of the following format:
#' O={W1 W2 C1 A1 LA1 Z1 LZ1 Y1 C2 A2 LA2 Z2 LZ2 Y2 ... }
#'
#' @param n Sample size
#' @param end.time Last time point
#' @param abar Set intervention
#' @param abar.prime Control intervention
#'
#' @return Returns a dataframe with simulated covariates as described in the description.
#'
#' @export GenerateData

GenerateData <- function(n, end.time, abar=NULL,abar.prime=NULL) {

  if(is.null(abar) != is.null(abar.prime)){stop('abar and abar.prime either both given or both not given')}
  if((!is.null(abar) & length(abar)!=end.time) | (!is.null(abar.prime) & length(abar.prime)!=end.time)){stop('abar and abar.prime both need to be length of end.time')}

  #Define some functions:
  rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

  CalcY <- function(W1,W2,A,LA,Z,LZ,prevLA=NULL) {
    lin <- .2 +1.5*W2+1*LA + .2*LZ -.3*A-.3*Z  - .2*A*Z
    if(!is.null(prevLA)) lin <- lin - .1*prevLA

    Y <- rexpit(lin)
    return(Y)
  }

  CalcZ <- function(W2, LA, A) {
    lin <- -.5+.8*W2+.8*A + 1*LA
    Z <- rexpit(lin)
    return(Z)
  }

  est.psi0 <- !is.null(abar) & !is.null(abar.prime)

  #Data generating process:

  W1 <- rbinom(n,1,.4)
  W2 <- rbinom(n,1,.6)

  LA <- LZ <- A <- C <- Z <- Y <- matrix(NA, nrow=n, ncol=end.time)

  uncensored.alive <- rep(TRUE, n)

  for (t in 1:(end.time)) {

    ### C
    if(est.psi0){
      C[uncensored.alive,t] <- 1
    }else{
      if(t==1) C[uncensored.alive,t] <- rexpit(1.5-.8*W2[uncensored.alive]-.4*W1[uncensored.alive])
      else C[uncensored.alive,t] <- rexpit(1.5-.8*W2[uncensored.alive]+.5*A[uncensored.alive,t-1]-.4*LZ[uncensored.alive,t-1])
    }

    #Update who died (note more might have died from outcome Y1 to now):
    uncensored.alive <- uncensored.alive & C[,t]

    ##A
    if (est.psi0) {
      A[uncensored.alive, t] <- abar[t]
    } else {
      if(t==1) A[uncensored.alive, t] <- rexpit(-.1 + (0.7*W1 + 1.2*W2)[uncensored.alive])
      if(t>1)A[uncensored.alive, t] <- rexpit(-.1 + 1.2*W2[uncensored.alive]+0.7*LZ[uncensored.alive, t-1] - .1*A[uncensored.alive, t-1])
    }

    ## LA
    if(t==1){
      LA[uncensored.alive,t] <- rexpit(-.8+.1*W1[uncensored.alive]+.3*W2[uncensored.alive]+A[uncensored.alive,t])
    }else{
      LA[uncensored.alive, t] <- rexpit(-.8+.1*LZ[uncensored.alive,t-1]+.3*LA[uncensored.alive,t-1]+A[uncensored.alive,t])
    }

    ## Z
    if (est.psi0) {
      useA <- abar.prime[t]
    }else{
      useA  <- A[uncensored.alive,t]
    }

    Z[uncensored.alive, t] <- CalcZ(W2 = W2[uncensored.alive],LA=LA[uncensored.alive,t],A = useA)

    ## LZ
    if(t==1) {
      LZ[uncensored.alive,t] <- rexpit(-1+.3*W2[uncensored.alive]+A[uncensored.alive,t]+.7*Z[uncensored.alive,t])
    }else {
      LZ[uncensored.alive,t] <- rexpit(-1+.3*W2[uncensored.alive]+A[uncensored.alive,t]+.7*Z[uncensored.alive,t] -.2*LZ[uncensored.alive,t-1])
    }

    ## Y
    if(t==1) {
      Y[uncensored.alive,t] <- CalcY(W1=W1[uncensored.alive],W2 = W2[uncensored.alive],A = A[uncensored.alive,t],LA = LA[uncensored.alive,t],Z=Z[uncensored.alive,t],LZ=LZ[uncensored.alive,t])
    }else{
      Y[uncensored.alive,t] <- CalcY(W1=W1[uncensored.alive],W2 = W2[uncensored.alive],A = A[uncensored.alive,t],LA = LA[uncensored.alive,t],Z=Z[uncensored.alive,t],LZ=LZ[uncensored.alive,t],prevLA=LA[uncensored.alive,t-1])
      Y[Y[, t-1]==1, t] <- 1 #if Y was 1 last period, it stays 1 (this doesn't explicitly condition on uncensored.alive, but if you died, you must have been uncensored; we set Y=1 even after death )
    }

    uncensored.alive <- uncensored.alive & (!Y[,t])
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

  .CreateDataFrame <- function(W1, W2, C, A, Z,LA, LZ, Y, end.time) {

    d <- data.frame(W1, W2)
    for (t in 1:(end.time)) {
      d <- data.frame(d, .BinaryToCensoring(is.uncensored=C[, t]), A[, t], LA[,t],Z[,t],LZ[,t],Y[,t])
      names(d)[ncol(d) - 5] <- paste0("C_", t)
      names(d)[ncol(d) - 4] <- paste0("A_", t)
      names(d)[ncol(d) - 3] <- paste0("LA_", t)
      names(d)[ncol(d) - 2] <- paste0("Z_", t)
      names(d)[ncol(d) - 1] <- paste0("LZ_", t)
      names(d)[ncol(d)] <- paste0("Y_", t)
    }
    return(d)
  }

  return(.CreateDataFrame(W1, W2, C = C,A = A,Z = Z,LA = LA,LZ = LZ,Y =  Y,end.time =  end.time))

}
