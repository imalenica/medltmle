#' CreateMedInputs
#'
#' Create Mediation Inputs for ltmle.
#'
#' @param data Dataframe containing the data in a wide format.
#' @param Anodes names of columns containing A covariates.
#' @param Znodes names of columns containing Z covariates.
#' @param Cnodes names of columns containing C covariates.
#' @param Lnodes names of columns containing L covariates.
#' @param Ynodes names of columns containing Y covariates.
#' @param survivalOutcome logical variable specifying if the outcome is survival.
#' @param QLform Q forms for L covariates.
#' @param QZform Q forms for Z covariates.
#' @param gform g form for intervention (if there is a censoring variable, include C as well).
#' @param qzform g form for Z covariates.
#' @param qLform g form for L covariates.
#' @param gbounds Bounds for the propensity score.
#' @param deterministic.g.function Logical specifying if g is a deterministic function.
#' @param stratify Logical enabling stratified outcome.
#' @param SL.library SuperLearner library for estimation.
#' @param estimate.time Measure time to fun function.
#' @param deterministic.Q.function Logical specifying if Q is a deterministic function.
#' @param gcomp
#' @param iptw.only Use IPTW estimator only
#' @param IC.variance.only Only estimate variance through the influence curve
#' @param observation.weights Provide weight for observations
#' @param Yrange Specify range for the outcome.
#' @param regimes Regimes for abar as output by the RegimesFromAbar function.
#' @param regimes.prime Regimes for abar.prime as output by the RegimesFromAbar function.
#' @param working.msm Working MSM as output by the GetMediationMSMInputsForLtmle function.
#' @param summary.measures Summary measures as output by the GetMediationMSMInputsForLtmle function.
#' @param final.Ynodes Final Y node(s).
#' @param msm.weights MSM weights.
#'
#' @return Returns output ready for ltmleMed.
#'

CreateMedInputs <- function(data, Anodes, Cnodes, Lnodes, Ynodes, Znodes, survivalOutcome, QLform, QZform, gform, qzform,qLform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, regimes.prime, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights) {

  if (is.list(regimes)) {

    if (!all(do.call(c, lapply(regimes, is.function)))) stop("If 'regimes' is a list, then all elements should be functions.")
    regimes <- aperm(simplify2array(lapply(regimes, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3))

  }

  if (is.list(regimes.prime)) {

    if (!all(do.call(c, lapply(regimes.prime, is.function)))) stop("If 'regimes.prme' is a list, then all elements should be functions.")
    regimes.prime <- aperm(simplify2array(lapply(regimes.prime, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3))

    }

  if (!(is.null(regimes) || dim(regimes) != 3)) {

    stop("regimes must be an array with 3 dimensions (unless Anodes is NULL, in which case regimes can be NULL)")

    }

  if (!(is.null(regimes.prime) || dim(regimes.prime) != 3)) {

    stop("regimes.prime must be an array with 3 dimensions (unless Anodes is NULL, in which case regimes can be NULL)")

    }

  if (is.null(regimes) || dim(regimes)[3]==0) {

    if (length(Anodes) != 0) {
      stop("regimes must not be NULL (or have dim(regimes)[3]==0) unless Anodes is also NULL")

    }

    regimes <- array(dim=c(nrow(data), 0, 1))
  }

  if (is.null(regimes.prime) || dim(regimes.prime)[3]==0) {

    if (length(Anodes) != 0) {
      stop("regimes.prime must not be NULL (or have dim(regimes)[3]==0) unless Anodes is also NULL")

    }

    regimes.prime <- array(dim=c(nrow(data), 0, 1))
  }

  if (is.logical(regimes)) {

    regimes <- regimes * 1
    message("abar or regimes was passed as logical and was converted to numeric")

    }

  if (is.logical(regimes.prime)) {

    regimes.prime <- regimes.prime * 1
    message("abar.prime or regimes.prime was passed as logical and was converted to numeric")

  }

  all.nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes,Znodes = Znodes)
  QLform <- CreateLYNodes(data, all.nodes, check.Qform=TRUE, Qform=QLform)$Qform  ### remove blocks
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions=!is.null(deterministic.g.function) && is.null(deterministic.Q.function))

  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(all.nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }

  #Using get to avoid the "no visible binding for global variable" note in R CMD check
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.Library")
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")

  if (is.null(summary.measures)) {
    summary.measures <- matrix(nrow=dim(regimes)[3], ncol=0)
  }

  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }

  if (is.null(observation.weights)) observation.weights <- rep(1, nrow(data))

  #error checking (also get value for survivalOutcome if NULL)
  #Need to fix later:
  #check.results <- CheckMediationInputs(data, all.nodes, survivalOutcome, QLform=QLform, QZform = QZform,gform=gform,qLform=qLform,qzform=qzform, gbounds, Yrange, deterministic.g.function, SL.library, regimes=regimes,regimes.prime = regimes.prime, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp, IC.variance.only)

  #survivalOutcome <- check.results$survivalOutcome

  if (!isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) {
    data <- CleanData(data, all.nodes, deterministic.Q.function, survivalOutcome)
  }

  transform.list <- TransformOutcomes(data, all.nodes, Yrange)
  data <- transform.list$data
  transformOutcome <- transform.list$transformOutcome
  #binaryOutcome <- check.results$binaryOutcome
  binaryOutcome <- all(unlist(data[, all.nodes$Y]) %in% c(0, 1, NA))

  if (is.null(QLform)) QLform <- GetDefaultFormMediation(data, all.nodes, is.Qform=TRUE, is.QLform = TRUE,is.qzform=FALSE, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(QZform)) QZform <- GetDefaultFormMediation(data, all.nodes, is.Qform=TRUE, is.QLform = FALSE,is.qzform=FALSE, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(qzform)) qzform <- GetDefaultFormMediation(data, all.nodes, is.Qform=FALSE, is.QLform = FALSE,is.qzform=TRUE, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(gform)) gform <- GetDefaultFormMediation(data, all.nodes, is.Qform=FALSE, is.QLform = FALSE,is.qzform=FALSE, stratify, survivalOutcome, showMessage=TRUE)

  # Several functions in the pooled version are only written to accept main terms MSM
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  main.terms <- ConvertToMainTerms(data, working.msm, summary.measures, all.nodes)

  intervention.match <- CalcInterventionMatchArray(data, regimes, all.nodes$A)
  intervention.match.prime <- CalcInterventionMatchArray(data, regimes.prime  , all.nodes$A)

  uncensored.array <- CalcUncensoredMatrix(data, all.nodes$C)

  #inputs <- list(data=data, all.nodes=all.nodes, survivalOutcome=survivalOutcome, QLform=QLform, QZform=QZform, gform=gform, qzform=qzform, qLform=qLform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, regimes.prime=regimes.prime,working.msm=main.terms$msm, combined.summary.measures=main.terms$summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, IC.variance.only=IC.variance.only, observation.weights=observation.weights, baseline.column.names=main.terms$baseline.column.names, beta.names=main.terms$beta.names, uncensored=check.results$uncensored, intervention.match=intervention.match, intervention.match.prime=intervention.match.prime)
  inputs <- list(data=data, all.nodes=all.nodes, survivalOutcome=survivalOutcome, QLform=QLform, QZform=QZform, gform=gform, qzform=qzform, qLform=qLform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, regimes.prime=regimes.prime,working.msm=main.terms$msm, combined.summary.measures=main.terms$summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, IC.variance.only=IC.variance.only, observation.weights=observation.weights, baseline.column.names=main.terms$baseline.column.names, beta.names=main.terms$beta.names, uncensored=uncensored.array, intervention.match=intervention.match, intervention.match.prime=intervention.match.prime)
  class(inputs) <- "ltmleInputs"
  return(inputs)
}
