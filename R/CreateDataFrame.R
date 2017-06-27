#' .CreateDataFrame
#'
#' Helper function for naming variables. Note the time ordering: W1,W2,C,A,Z,LA,LZ,Y
#'
#' @param W1 Baseline covariates.
#' @param W2 Baseline covariates.
#' @param C Censoring variable.
#' @param A Exposure/treatment covariate.
#' @param Z Mediator covariate.
#' @param LA L covariate, influenced by A
#' @param LZ L covariate, influenced by Z
#' @param Y Y outcome
#' @param end.time final time point.
#'

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
