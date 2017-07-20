################################
# medltmle
################################

#' medltmle
#'
#' Estimates parameters for longitudinal mediation analysis with time-varying mediators.
#'
#' @param data Dataframe containing the data in a wide format.
#' @param Anodes names of columns containing A covariates (exposure) (character).
#' @param Znodes names of columns containing Z covariates (mediator) (character).
#' @param Cnodes names of columns containing C covariates (censoring) (character).
#' @param Lnodes names of columns containing L covariates (covariate) (character).
#' @param Ynodes names of columns containing Y covariates (outcome) (character).
#' @param W2nodes names of columns containing W2 covariates (baseline covariates in need of fluctuation) (character).
#' @param Dnodes names of columns containing D covariates (death indicator) (character).
#' @param survivalOutcome logical variable specifying if the outcome is survival.
#' @param QLform Q forms for L covariates.
#' @param QZform Q forms for Z covariates.
#' @param gform g form for intervention (if there is a censoring variable, include C as well).
#' @param qzform g form for Z covariates.
#' @param qLform g form for L covariates.
#' @param abar Intervention. Dimension should correspond to the number of intervention nodes available.
#' @param abar.prime Control for the exposure. Dimension should correspond to the number of intervention nodes available.
#' @param gbounds Bounds for the propensity score.
#' @param deterministic.g.function Logical specifying if g is a deterministic function.
#' @param stratify Logical enabling stratified outcome.
#' @param SL.library SuperLearner library for estimation.
#' @param estimate.time Measure time to fun function.
#' @param deterministic.Q.function Logical specifying if Q is a deterministic function.
#' @param gcomp Logical indicating whether to use Gcomp instead (no updating if TRUE).
#' @param iptw.only Use IPTW estimator only
#' @param IC.variance.only Only estimate variance through the influence curve
#' @param observation.weights Provide weight for observations
#' @param rule Specify rule for the intervention.
#' @param Yrange Specify range for the outcome.
#' @param estimand Specifies which estimand to estimate. Options are: natural effect (NE), stochastic effect (SE), or controlled effect (CE).
#' @param time.end How many time points in the longitudinal data?
#'
#' @return Returns estimate of \eqn{E[Y_{\tau}(a, \overline{\Gamma}^{a^'})]}
#'
#' @export medltmle

medltmle <- function(data, Anodes, Znodes, Cnodes=NULL, Lnodes=NULL, Ynodes, W2nodes=NULL,Dnodes=NULL,
                           survivalOutcome=NULL,
                           QLform=NULL, QZform=NULL,gform=NULL, qzform=NULL, qLform=NULL,
                           abar, abar.prime,  rule=NULL, gbounds=c(0.01, 1), Yrange=NULL,
                           deterministic.g.function=NULL, deterministic.Q.function=NULL,
                           stratify=FALSE, SL.library=NULL,
                           estimate.time=TRUE, gcomp=FALSE,
                           iptw.only=FALSE, IC.variance.only=FALSE, observation.weights=NULL, estimand="SE", time.end, past) {

  #Implement rule and deterministic g function option. TO DO.
  if(!is.null(rule))stop('rule option not implemented yet')
  if(!is.null(deterministic.g.function)) stop('deterministic.g.function option not implemented yet')

  #Note that MSM is not implemented, yet. TO DO.
  msm.inputs <- CreateMedMSMInputs(data, abar=abar, abar.prime = abar.prime, rule = rule, gform=gform)

  inputs <- CreateMedInputs(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, Znodes=Znodes, Dnodes=Dnodes,W2nodes=W2nodes,
                                  QLform=QLform, QZform=QZform, gform=msm.inputs$gform, qLform=qLform, qzform=qzform,
                                  Yrange=Yrange, gbounds=gbounds, SL.library=SL.library, stratify=stratify,
                                  regimes=msm.inputs$regimes, regimes.prime=msm.inputs$regimes.prime,
                                  working.msm=msm.inputs$working.msm, summary.measures=msm.inputs$summary.measures,
                                  final.Ynodes=msm.inputs$final.Ynodes, msm.weights=msm.inputs$msm.weights,
                                  estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only,
                                  deterministic.Q.function=deterministic.Q.function, deterministic.g.function=deterministic.g.function,
                                  IC.variance.only=IC.variance.only,
                                  observation.weights=observation.weights, survivalOutcome=survivalOutcome, estimand=estimand, past=time.end, time.end=time.end)
  #fixme
  print(tracemem(inputs))
  result <- ltmleMediation(inputs)
  result$call <- match.call()
  return(result)

}



