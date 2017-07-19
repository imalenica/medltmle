#' CreateMedInputs
#'
#' Create Mediation Inputs for ltmleMediation.
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
#' @param gbounds Bounds for the propensity score.
#' @param deterministic.g.function Logical specifying if g is a deterministic function.
#' @param stratify Logical enabling stratified outcome.
#' @param SL.library SuperLearner library for estimation.
#' @param estimate.time Measure time to fun function.
#' @param deterministic.Q.function Logical specifying if Q is a deterministic function.
#' @param gcomp gcomp formula.
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
#' @param estimand Specifies which estimand to estimate. Options are: natural effect (NE), stochastic effect (SE), or controlled effect (CE).
#' @param past Number indicating Markov order for the conditional densities.
#' @param time.end Total number of time points.
#'
#'
#' @return Returns output ready for ltmleMediation.
#'
#' @export CreateMedInputs
#'

CreateMedInputs <- function(data, Anodes, Cnodes, Lnodes, Ynodes, Znodes, Dnodes, W2nodes, survivalOutcome, QLform, QZform, gform, qzform, qLform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, regimes.prime, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights, estimand, past, time.end) {

  if (is.list(regimes)) {

    if (!all(do.call(c, lapply(regimes, is.function)))) stop("If 'regimes' is a list, then all elements should be functions.")
    regimes <- aperm(simplify2array(lapply(regimes, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3))

  }

  if (is.list(regimes.prime)) {

    if (!all(do.call(c, lapply(regimes.prime, is.function)))) stop("If 'regimes.prime' is a list, then all elements should be functions.")
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

  #Sort nodes
  all.nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes, Znodes, Dnodes, W2nodes)
  #remove blocks
  QLform <- CreateLYNodes(data, all.nodes, check.Qform=TRUE, Qform=QLform)$Qform
  #Convert censoring nodes into factors
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions=!is.null(deterministic.g.function) && is.null(deterministic.Q.function))

  #Set final Y node, if not specified:
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(all.nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }

  #Using get to avoid the "no visible binding for global variable" note in R CMD check
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.Library")
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")

  #Set summary.measures if not already set:
  if (is.null(summary.measures)) {
    summary.measures <- matrix(nrow=dim(regimes)[3], ncol=0)
  }

  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }

  #Set observation weights; if not specified, assign equal weight to all.
  if (is.null(observation.weights)) observation.weights <- rep(1, nrow(data))

  #error checking (also get value for survivalOutcome if NULL)
  #TO DO: fix
  #check.results <- CheckMediationInputs(data, all.nodes, survivalOutcome, QLform=QLform, QZform = QZform,gform=gform,qLform=qLform,qzform=qzform, gbounds, Yrange, deterministic.g.function, SL.library, regimes=regimes,regimes.prime = regimes.prime, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp, IC.variance.only)

  #survivalOutcome <- check.results$survivalOutcome

  if (!isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) {
    data <- CleanData(data, all.nodes, deterministic.Q.function, survivalOutcome, showMessage = FALSE)
  }

  #Transform the output to be in the 0-1 range. Get the Y range, if not specified.
  transform.list <- TransformOutcomes(data, all.nodes, Yrange)
  data <- transform.list$data
  transformOutcome <- transform.list$transformOutcome
  #binaryOutcome <- check.results$binaryOutcome
  #Check if binary, as expected.
  binaryOutcome <- all(unlist(data[, all.nodes$Y]) %in% c(0, 1, NA))

  #If QLform, QZform, qzform and gform are not specified, return default form.
  #Each formula will consist of all parent nodes except censoring (both C and D, if available) and event nodes.

  if (is.null(qLform)) qLform <- GetDefaultFormMediation(data, all.nodes, is.Qform=FALSE, is.QLform = FALSE,is.qzform=FALSE, is.qLform=TRUE, past=past, time.end=time.end, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(qzform)) qzform <- GetDefaultFormMediation(data, all.nodes, is.Qform=FALSE, is.QLform = FALSE,is.qzform=TRUE, is.qLform=FALSE, past=past, time.end=time.end, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(gform)) gform <- GetDefaultFormMediation(data, all.nodes, is.Qform=FALSE, is.QLform = FALSE,is.qzform=FALSE, is.qLform=FALSE, past=past, time.end=time.end, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(QZform)) QZform <- GetDefaultFormMediation(data, all.nodes, is.Qform=TRUE, is.QLform = FALSE,is.qzform=FALSE, is.qLform=FALSE, past=past, time.end=time.end, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(QLform)) QLform <- GetDefaultFormMediation(data, all.nodes, is.Qform=TRUE, is.QLform = TRUE,is.qzform=FALSE, is.qLform=FALSE, past=past, time.end=time.end, stratify, survivalOutcome, showMessage=TRUE)

  # Convert to main terms MSM.
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S2 + S3 + S4" where
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  main.terms <- ConvertToMainTerms(data, working.msm, summary.measures, all.nodes)

  #Does A observed in the data match the regime?
  intervention.match <- CalcInterventionMatchArray(data, regimes, all.nodes$A)
  intervention.match.prime <- CalcInterventionMatchArray(data, regimes.prime  , all.nodes$A)

  #Check if the patient is censored at one of the C nodes.
  uncensored.array <- CalcUncensoredMatrix(data, all.nodes$C)

  #If there is a D node, carry forward death.
  if(!is.null(Dnodes)){

    node.set <- all.nodes$D

    for(i in 2:length(Dnodes)){

      prev.node <- node.set[i-1]
      cur.node <- node.set[i]

      ind<-data[,prev.node]==1
      ind2<-ind %in% TRUE

      data[ind2,cur.node]<-1

    }
  }

  #inputs <- list(data=data, all.nodes=all.nodes, survivalOutcome=survivalOutcome, QLform=QLform, QZform=QZform, gform=gform, qzform=qzform, qLform=qLform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, regimes.prime=regimes.prime,working.msm=main.terms$msm, combined.summary.measures=main.terms$summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, IC.variance.only=IC.variance.only, observation.weights=observation.weights, baseline.column.names=main.terms$baseline.column.names, beta.names=main.terms$beta.names, uncensored=check.results$uncensored, intervention.match=intervention.match, intervention.match.prime=intervention.match.prime)
  inputs <- list(data=data, all.nodes=all.nodes, survivalOutcome=survivalOutcome, QLform=QLform, QZform=QZform, gform=gform, qzform=qzform, qLform=qLform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, regimes.prime=regimes.prime,working.msm=main.terms$msm, combined.summary.measures=main.terms$summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, IC.variance.only=IC.variance.only, observation.weights=observation.weights, baseline.column.names=main.terms$baseline.column.names, beta.names=main.terms$beta.names, uncensored=uncensored.array, intervention.match=intervention.match, intervention.match.prime=intervention.match.prime, estimand=estimand)
  class(inputs) <- "medltmleInputs"
  return(inputs)

}

################################
# GetDefaultFormMediation
################################

#' GetDefaultFormMediation
#'
#' If QLform, QZform, qzform and gform are not specified, return default form.
#' Each formula consists of all parent nodes except censoring and event nodes.
#' If \code{stratify}=TRUE, do not include A nodes.
#'
#' @param data Available data in a \code{Data.Frame} format.
#' @param nodes List of available nodes, as created by \code{CreateNodes}.
#' @param is.Qform Logical indicating whether to specify general Q formula.
#' @param is.QLform Logical indicating whether to specify Q formula for covariates.
#' @param is.qzform Logical indicating whether to specify general Q formula for Z.
#' @param is.qLform Logical indicating whether to specify general Q formula for L.
#' @param past Number indicating Markov order for the conditional densities.
#' @param time.end Total number of time points.
#' @param stratify Logical indicating whether to straify.
#' @param survivalOutcome Logical indicating if the outcome a survival outcome.
#' @param showMessage Logical indicating whether to show comments while executing.
#'
#' @return Returns default Q or g formula if not specified.
#'

GetDefaultFormMediation <- function(data, nodes, is.Qform, is.QLform, is.qzform, is.qLform, past, time.end, stratify, survivalOutcome, showMessage) {

  if (is.Qform) {
    if(is.QLform){
      lhs <- rep("Q.kplus1", length(nodes$L))
      node.set <- nodes$L
    }else{
      lhs <- rep("Q.kplus1", length(nodes$Z))
      node.set <- nodes$Z
    }
  } else {
    if(is.qzform){
      lhs <- names(data)[nodes$Z]
      node.set <- nodes$Z
    }else if(is.qLform){
      lhs <- names(data)[nodes$LY]
      node.set <- nodes$LY
    }else{
      lhs <- names(data)[nodes$AC]
      node.set <- nodes$AC
    }
  }

  if (stratify) {
    stratify.nodes <- c(nodes$C, nodes$A)
  } else {
    stratify.nodes <- c(nodes$C)
  }

  if (survivalOutcome) {
    stratify.nodes <- c(stratify.nodes, nodes$Y)
  }

  #Even if not survival outcome, in case there is survival-type node D, remove it.
  if(!is.null(nodes$D)){
    stratify.nodes <- c(stratify.nodes, nodes$D)
  }

  form <- NULL

  #First, separate the baseline covariates:
  base<-length(nodes$baseline)
  baseline.node.names<-names(data)[1:base]

  #Check out how many unique nodes per time:
  group<-length(grep("_1",names(data)[(base+1):ncol(data)]))

  for (i in seq_along(node.set)) {

    cur.node <- node.set[i]

    if (cur.node == 1) {
      #no parent nodes
      form[i] <- paste(lhs[i], "~ 1")
    } else if(past!=time.end){

      if((base+group) >= cur.node){

        parent.node.names <- names(data)[setdiff(1:(cur.node - 1), stratify.nodes)]

      }else{

        parent.node.names <- names(data)[setdiff((cur.node - 1):(cur.node-group), stratify.nodes)]

        #Check for A:
        if(length(grep("^A", parent.node.names))==0){

          #Pick the closest A
          Anode.index <- which(nodes$A < cur.node)
          parent.node.names<-c(names(data[nodes$A[Anode.index]]), parent.node.names)

        }
      }

      if (length(parent.node.names) == 0) {
        form[i] <- paste(lhs[i], "~ 1")
      } else {
        form[i] <- paste(lhs[i], "~", paste(parent.node.names, collapse=" + "))
      }

      #Include all covariates.
      }else if(past == time.end){
      parent.node.names <- names(data)[setdiff(1:(cur.node - 1), stratify.nodes)]

      if (length(parent.node.names) == 0) {
        form[i] <- paste(lhs[i], "~ 1")
      } else {
        form[i] <- paste(lhs[i], "~", paste(parent.node.names, collapse=" + "))
      }
      }

    names(form)[i] <- names(data)[cur.node]
  }

  if (showMessage) {
    #Prints formulas with automatic wrapping thanks to print.formula
    message(ifelse(is.Qform, "Qform", "gform"),
            " not specified, using defaults:")
    lapply(seq_along(form), function(i, names) {
      message("formula for ", names[i], ":")
      #Using print on a formula because it nicely wraps
      message(capture.output(print(as.formula(form[i]), showEnv=FALSE)))
    }, names=names(form))
    message("")
  }
  return(form)
}
