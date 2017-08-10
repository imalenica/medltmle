################################
# ltmleMediation
################################

#' ltmleMediation
#'
#' Prepares output obtained from \code{LtmleMediationMSMFromInputs} and \code{MainCalcsMediation} for final report by
#' \code{medltmle}.
#'
#' @param inputs Output of \code{CreateMedInputs}
#'
#' @return Returns $psi^{IPW}_n$, $psi^{TMLE}_n$ and if specified $psi^{NTsub}_n$ estimates,
#' corresponding influence curve, cumulative conditional density with abar and abar.prime rule for g, QL and QZ (bounded and unbounded),
#' fits for Q and Qstar and corresponding formulas, final Qstar estimates, logical indicating binary and/or transformed outcome.
#'
#' @export ltmleMediation

ltmleMediation <- function(inputs) {

  #Return main TMLE calculations and add few results from CreateMedInputs.
  msm.result <- LtmleMediationMSMFromInputs(inputs)
  #Number of regimes
  num.regimes <- dim(inputs$regimes)[3]
  stopifnot(num.regimes %in% 1:2)

  if (num.regimes == 2) {
    class(msm.result) <- "medltmleEffectMeasures"
    return(msm.result)
  }

  #IPTW
  names(msm.result$beta.iptw) <- names(msm.result$beta) <- NULL
  iptw <- plogis(msm.result$beta.iptw)
  iptw.list <- list(iptw.estimate=iptw, iptw.IC=iptw*(1-iptw)*msm.result$IC.iptw[, 1])

  r <- list()

  #TMLE
  if (inputs$iptw.only) {
    tmle <- NA
    tmle.IC <- rep(NA, nrow(inputs$data))
  } else {
    tmle <- plogis(msm.result$beta)
    tmle.IC <- msm.result$IC[, 1] #only one regime
  }

  #Estimates
  r$estimates <- c(tmle=tmle, iptw=iptw.list$iptw.estimate)
  #Influence curve
  r$IC <- list(tmle=tmle.IC * tmle * (1 - tmle), iptw=iptw.list$iptw.IC)

  if (!is.null(msm.result$variance.estimate)) {

    stopifnot(length(msm.result$variance.estimate) == 1)
    r$variance.estimate <- msm.result$variance.estimate[1] * (tmle * (1 - tmle))^2

  }

  #If gcomp=TRUE, report untargeted
  if (inputs$gcomp) {
    names(r$estimates)[1] <- names(r$IC)[1] <- "gcomp"
  }

  #Report bounded and unbounded cumulative conditional densities of g,QL and QZ
  r$cum.gq <- list(cum.gq.abar=lapply(msm.result$cum.gq.abar,function(tt)AsMatrix(tt[, , 1])),
                   cum.gq.abar.prime=lapply(msm.result$cum.gq.abar.prime,function(tt)AsMatrix(tt[, , 1])))
  r$cum.gq.unbounded <- list(cum.gq.abar=lapply(msm.result$cum.gq.abar.unbounded,function(tt)AsMatrix(tt[, , 1])),
                             cum.gq.abar.prime=lapply(msm.result$cum.gq.abar.prime.unbounded,function(tt)AsMatrix(tt[, , 1])))

  #only one regime
  r$gcomp <- inputs$gcomp
  r$fit <- msm.result$fit
  r$fit$g <- r$fit$g[[1]]  #only one regime
  r$fit$Q <- r$fit$Q[[1]]  #only one regime
  r$Qstar <- msm.result$Qstar[, 1, 1] #1 regime, 1 final.Ynode

  r$formulas <- msm.result$formulas
  r$binaryOutcome <- msm.result$binaryOutcome
  r$transformOutcome <- msm.result$transformOutcome==TRUE #Want to store transformOutcome flag without attributes

  if (msm.result$transformOutcome) {

    Yrange <- attr(msm.result$transformOutcome, "Yrange")
    #back transform estimate and IC
    r$estimates <- r$estimates*diff(Yrange) + min(Yrange)
    r$IC <- lapply(r$IC, function (IC) IC * diff(Yrange))
    r$variance.estimate <- r$variance.estimate * (diff(Yrange))^2

  }

  class(r) <- "medltmle"
  return(r)

}

################################
# LtmleMediationMSMFromInputs
################################

#' LtmleMediationMSMFromInputs
#'
#' Calls \code{MainCalcsMediation} and edits the results obtained.
#'
#' @param inputs Output of \code{CreateMedInputs}.
#'
#' @return Returns results from \code{MainCalcsMediation}, and adds Gcomp, formulas,
#' biinary outcome, transformed outcome and survival outcome logical to the results from \code{CreateMedInputs}.
#'

LtmleMediationMSMFromInputs <- function(inputs) {

  #main TMLe calculations from inputs.
  result <- MainCalcsMediation(inputs)

  #Add to results from inputs:

  result$gcomp <- inputs$gcomp
  result$formulas <- list(QLform=inputs$QLform, QZform=inputs$QZform, gform=inputs$gform, qzform=inputs$qzform,qLform=inputs$qLform)
  result$binaryOutcome <- inputs$binaryOutcome
  result$transformOutcome <- inputs$transformOutcome
  result$survivalOutcome <- inputs$survivalOutcome

  return(result)

}
