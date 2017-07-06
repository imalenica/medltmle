#' CreateMedMSMInputs
#'
#' Create Mediation MSM Inputs for ltmle.
#'
#' @param data Dataframe containing the data in a wide format.
#' @param abar Intervention.
#' @param abar.prime Control for the exposure.
#' @param rule Specify rule for the intervention.
#' @param gform g form for intervention (if there is a censoring variable, include C as well).
#'
#' @return Returns regimes for abar and abar.prime in proper format, working MSM and MSM weights, summary measures, gform and final Y nodes.
#'
#' @export CreateMedMSMInputs
#'

CreateMedMSMInputs <- function(data, abar, abar.prime, rule, gform) {

  #Options for specified rule. Not yet implemented.
  if ((!missing(abar) && is.list(abar)) || is.list(rule)) {

    if (is.list(rule)) {

      if (length(rule) != 2) stop("If rule is a list, it must be of length 2")

      regimes1 <- RegimesFromAbar(data, abar, rule[[1]])
      regimes0 <- RegimesFromAbar(data, abar, rule[[2]])

    } else {
      if (length(abar) != 2) stop("If abar is a list, it must be of length 2")

      regimes1 <- RegimesFromAbar(data, abar[[1]], rule)
      regimes0 <- RegimesFromAbar(data, abar[[2]], rule)
      regimes1.prime <- RegimesFromAbar(data, abar.prime[[1]], rule)
      regimes0.prime <- RegimesFromAbar(data, abar.prime[[2]], rule)

    }

    if (ncol(regimes1) != ncol(regimes0)) stop("If abar or rule is a list, both elements must give a matrix with the same number of columns")
    if (nrow(regimes1) != nrow(regimes0)) stop("If abar or rule is a list, both elements must give a matrix with the same number of rows")
    if (ncol(regimes1.prime) != ncol(regimes0.prime)) stop("If abar or rule is a list, both elements must give a matrix with the same number of columns")
    if (nrow(regimes1.prime) != nrow(regimes0.prime)) stop("If abar or rule is a list, both elements must give a matrix with the same number of rows")
    regimes <- c(regimes1, regimes0)
    regimes.prime <- c(regimes1.prime, regimes0.prime)
    dim(regimes) <- c(nrow(regimes1), ncol(regimes1), 2)
    dim(regimes.prime) <- c(nrow(regimes1.prime), ncol(regimes1.prime), 2)
    summary.measures <- array(1:0, dim=c(2, 1, 1))
    colnames(summary.measures) <- "A"
    working.msm <- "Y ~ A"
    msm.weights <- matrix(1, nrow=2, ncol=1)

  } else {

    regimes <- RegimesFromAbar(data, abar, rule)
    regimes.prime <- RegimesFromAbar(data, abar.prime, rule)
    working.msm <- "Y ~ 1"
    msm.weights <- matrix(1, nrow=1, ncol=1)
    summary.measures <- array(dim=c(1, 0, 1))

  }

  if (is.numeric(gform)) {

    stopifnot(is.matrix(gform))
    dim(gform) <- c(nrow(gform), ncol(gform), 1)

  }

  msm.inputs <- list(regimes=regimes, regimes.prime=regimes.prime, working.msm=working.msm, summary.measures=summary.measures, gform=gform, final.Ynodes=NULL, msm.weights=msm.weights)
  return(msm.inputs)

}
