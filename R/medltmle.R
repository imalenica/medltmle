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
#' @param Inodes names of columns containing I covariates (instrument) (character).
#' @param Ynodes names of columns containing Y covariates (outcome) (character).
#' @param W2nodes names of columns containing W2 covariates (baseline covariates in need of fluctuation) (character).
#' @param Dnodes names of columns containing D covariates (death indicator) (character).
#' @param survivalOutcome If TRUE, then Y nodes are indicators of an event, and if Y at some time point is 1, then all following should be 1. Required to be TRUE or FALSE if outcomes are binary and there are multiple Ynodes.
#' @param QLform character vector of regression formulas for Q corresponding to L covariates.
#' @param QZform character vector of regression formulas for Q corresponding to Z covariates.
#' @param gform character vector of regression formulas for g or a matrix/array of prob(A=1).
#' @param qzform g form for Z covariates.
#' @param qLform g form for L covariates.
#' @param abar binary vector (numAnodes x 1) or matrix (n x numAnodes) of counterfactual treatment or a list of length 2.
#' @param abar.prime binary vector (numAnodes x 1) or matrix (n x numAnodes) of counterfactual treatment or a list of length 2.
#' @param gbounds lower and upper bounds on estimated cumulative probabilities for g-factors. Vector of length 2, order unimportant.
#' @param deterministic.g.function optional information on A and C nodes that are given deterministically. Default NULL indicates no deterministic links.
#' @param stratify if TRUE stratify on following abar when estimating Q and g. If FALSE, pool over abar.
#' @param SL.library optional character vector of libraries to pass to SuperLearner. NULL indicates glm should be called instead of SuperLearner. 'default' indicates a standard set of libraries. May be separately specified for Q and g.
#' @param estimate.time if TRUE, run an initial estimate using only 50 observations and use this to print a very rough estimate of the total time to completion. No action if there are fewer than 50 observations.
#' @param deterministic.Q.function optional information on Q given deterministically. See 'Details'. Default NULL indicates no deterministic links.
#' @param gcomp if TRUE, run the maximum likelihood based G-computation estimate instead of TMLE.
#' @param iptw.only by default (iptw.only = FALSE), both TMLE and IPTW are run in ltmle and ltmleMSM. If iptw.only = TRUE, only IPTW is run, which is faster.
#' @param IC.variance.only Only estimate variance through the influence curve
#' @param observation.weights observation (sampling) weights. Vector of length n. If NULL, assumed to be all 1.
#' @param rule a function to be applied to each row (a named vector) of data that returns a numeric vector of length numAnodes.
#' @param Yrange NULL or a numerical vector where the min and max of Yrange specify the range of all Y nodes.
#' @param CSE Logical specifying if the estimand is estimated by fully conditioning on the past (TRUE), or with the data-dependent estimate (FALSE).
#' @param time.end How many time points in the longitudinal data?
#' @param YisL Logical indicating whether Y is a function of time-varying covariate.
#' @param past Number indicating Markov order for the conditional densities.
#'
#' @return Returns estimate of \eqn{E[Y_{\tau}(a, \overline{\Gamma}^{a^'})]}
#'
#' @export medltmle

medltmle <- function(data, Anodes, Znodes, Cnodes=NULL, Lnodes=NULL, Ynodes, Inodes=NULL, W2nodes=NULL,Dnodes=NULL,
                           survivalOutcome=NULL,
                           QLform=NULL, QZform=NULL,gform=NULL, qzform=NULL, qLform=NULL,
                           abar, abar.prime,  rule=NULL, gbounds=c(0.01, 1), Yrange=NULL,
                           deterministic.g.function=NULL, deterministic.Q.function=NULL,
                           stratify=FALSE, SL.library=NULL,
                           estimate.time=TRUE, gcomp=FALSE,
                           iptw.only=FALSE, IC.variance.only=FALSE, observation.weights=NULL, CSE, time.end, past=1, YisL=TRUE) {

  #Implement rule and deterministic g function option. TO DO.
  if(!is.null(rule))stop('rule option not implemented yet')
  if(!is.null(deterministic.g.function)) stop('deterministic.g.function option not implemented yet')

  #Note that MSM is not implemented, yet. TO DO.
  msm.inputs <- CreateMedMSMInputs(data, abar=abar, abar.prime = abar.prime, rule = rule, gform=gform)

  inputs <- CreateMedInputs(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, Inodes=Inodes, Znodes=Znodes, Dnodes=Dnodes,W2nodes=W2nodes,
                                  QLform=QLform, QZform=QZform, gform=msm.inputs$gform, qLform=qLform, qzform=qzform,
                                  Yrange=Yrange, gbounds=gbounds, SL.library=SL.library, stratify=stratify,
                                  regimes=msm.inputs$regimes, regimes.prime=msm.inputs$regimes.prime,
                                  working.msm=msm.inputs$working.msm, summary.measures=msm.inputs$summary.measures,
                                  final.Ynodes=msm.inputs$final.Ynodes, msm.weights=msm.inputs$msm.weights,
                                  estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only,
                                  deterministic.Q.function=deterministic.Q.function, deterministic.g.function=deterministic.g.function,
                                  IC.variance.only=IC.variance.only,
                                  observation.weights=observation.weights, survivalOutcome=survivalOutcome, CSE=CSE, past=past, time.end=time.end, YisL=YisL)
  #fixme
  #Rprofmem(filename = "Rprofmem.out", append = FALSE, threshold = 0)
  #print(tracemem(inputs))
  result <- ltmleMediation(inputs)
  result$call <- match.call()
  return(result)

}



