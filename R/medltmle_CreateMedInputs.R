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
#'
#' @return Returns output ready for ltmleMed.
#'
#' @export CreateMedInputs
#'

###################################################################################################
# Note to myself:
#
# CreateMedInputs contains the following functions within it:
#
# NodeToIndex()
# CreateLYNodes()
# CreateNodes()
# GetLibrary()
# CleanData()
# TransformOutcomes()
# ConvertToMainTerms()
# CalcInterventionMatchArray()
# ConvertCensoringNodes
# IsDeterministic()
# RhsVars()
# CalcUncensoredMatrix()
###################################################################################################

CreateMedInputs <- function(data, Anodes, Cnodes, Lnodes, Ynodes, Znodes, survivalOutcome, QLform, QZform, gform, qzform,qLform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, regimes.prime, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights) {

  # Calc logical matrix - n x numCnodes = is uncensored up to and including Cnode[i]
  CalcUncensoredMatrix <- function(data, Cnodes) {
    uncensored <- matrix(nrow=nrow(data), ncol=length(Cnodes))
    cum.uncensored <- rep(TRUE, nrow(data))
    for (Cnode.index in seq_along(Cnodes)) {
      cum.uncensored <- cum.uncensored & (data[, Cnodes[Cnode.index]] %in% c("uncensored", NA))
      uncensored[, Cnode.index] <- cum.uncensored
    }
    return(uncensored)
  }

  # Return the right hand side variables of formula f as a character vector
  RhsVars <- function(f) {
    f <- as.formula(f)
    return(all.vars(f[[3]]))
  }

  # Determine which patients have died or have Q set deterministically by user function before cur.node
  # return list:
  #    is.deterministic: vector of [numObservations x 1] - true if patient is already dead before cur.node or set by deterministic.Q.function
  #    Q.value: vector of [which(is.deterministic) x 1] - value of Q

  IsDeterministic <- function(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g, survivalOutcome) {
    #set Q.value to 1 if previous y node is 1
    if (survivalOutcome && any(nodes$Y < cur.node)) {
      last.Ynode <- max(nodes$Y[nodes$Y < cur.node])
      is.deterministic <- data[, last.Ynode] %in% TRUE
    } else {
      is.deterministic <- rep(FALSE, nrow(data))
    }
  }

  # Convert named nodes to indicies of nodes
  NodeToIndex <- function(data, node) {
    if (! is.data.frame(data)) stop("data must be a data frame")
    if (is.numeric(node) || is.null(node)) return(node)
    if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
    index <- match(node, names(data))
    if (anyNA(index)) {
      stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
    }
    return(index)
  }

  # Get the LY nodes but don't include "blocks" of L/Y nodes uninterrupted by A/C/Z nodes
  CreateLYNodes <- function(data, nodes, check.Qform, Qform) {
    LYnodes <- sort(c(nodes$L, nodes$Y))
    #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
    nodes.to.remove <- NULL
    if (length(LYnodes) > 1) {
      for (i in 1:(length(LYnodes) - 1)) {
        cur.node <- LYnodes[i]
        next.node <- LYnodes[i + 1]
        if (! any(cur.node:next.node %in% nodes$ACZ)) {
          nodes.to.remove <- c(nodes.to.remove, next.node)
        }
      }
    }
    new.LYnodes <- setdiff(LYnodes, nodes.to.remove)
    if (check.Qform) {
      removed.Qform.index <- NULL
      for (i in nodes.to.remove) {
        index <- which(names(Qform) == names(data)[i])
        if (length(index) > 0) {
          removed.Qform.index <- c(removed.Qform.index, index)
        }
      }
      if (! is.null(removed.Qform.index)) {
        message("L/Y nodes (after removing blocks)  : ", names(data)[new.LYnodes], "\n")
        message("Qform names                        : ", names(Qform), "\n")
        message(paste("The following nodes are not being considered as L/Y nodes because they are part of a block of L/Y nodes. They are being dropped from Qform:\n"), paste(names(Qform)[removed.Qform.index], "\n", collapse=" "))
        Qform <- Qform[-removed.Qform.index]
      }
      return(list(LYnodes=new.LYnodes, Qform=Qform))
    }
    return(new.LYnodes)
  }

  CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes,Znodes=NULL) {
    Anodes <- NodeToIndex(data, Anodes)
    Cnodes <- NodeToIndex(data, Cnodes)
    Lnodes <- NodeToIndex(data, Lnodes)
    Ynodes <- NodeToIndex(data, Ynodes)
    nodes <- list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, AC=sort(c(Anodes, Cnodes)))
    if(!is.null(Znodes)){
      Znodes <- NodeToIndex(data, Znodes)
      nodes$Z <- Znodes
    }
    nodes$ACZ <- sort(c(Anodes, Cnodes, Znodes))
    nodes$baseline <- seq(1, min(c(nodes$A, nodes$L, nodes$C, nodes$Y,nodes$Z)) - 1)
    nodes$LY <- CreateLYNodes(data, nodes, check.Qform=FALSE)
    return(nodes)
  }

  GetLibrary <- function(SL.library, estimate.type) {
    if (is.null(names(SL.library))) return(SL.library)
    if (! identical(sort(names(SL.library)), sort(c("Q", "g")))) stop("If SL.library has names, it must have two names: Q and g")
    if (! estimate.type %in% c("Q", "g")) stop("bad estimate.type")
    if (length(setdiff(names(attributes(SL.library)), c("names", "return.fit"))) > 0) stop("If SL.library has attributes, the only valid attributes are name and return.fit")
    lib <- SL.library[[estimate.type]]
    attr(lib, "return.fit") <- attr(SL.library, "return.fit", exact = TRUE)
    return(SL.library[[estimate.type]])
  }

  # Set all nodes (except Y) to NA after death or censoring; Set Y nodes to 1 after death
  CleanData <- function(data, nodes, deterministic.Q.function, survivalOutcome, showMessage=TRUE) {
    #make sure binaries have already been converted before calling this function
    is.nan.df <- function (x) {
      y <- if (length(x)) {
        do.call("cbind", lapply(x, "is.nan"))
      } else {
        matrix(FALSE, length(row.names(x)), 0)
      }
    }

    is.na.strict <- function (x) is.na(x) & !is.nan.df(x)  #only for data.frames
    changed <- FALSE
    ua <- rep(TRUE, nrow(data))  #uncensored and alive
    if (ncol(data) == 1) return(data)
    deterministic.Q.function.depends.on.called.from.estimate.g <- length(grep("called.from.estimate.g", as.character(body(deterministic.Q.function)))) > 0

    #ISSUEEEE
    for (i in 1:(ncol(data)-1)) {
      if (anyNA(data[ua, 1:i])) stop("NA values are not permitted in data except after censoring or a survival event")
      #is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome)$is.deterministic #check determinisitic including node i
      is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome) #check determinisitic including node i

      if (deterministic.Q.function.depends.on.called.from.estimate.g) {
        #is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)$is.deterministic
        is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)
        if (any(is.deterministic[ua] & !is.deterministic.Q[ua])) stop("Any row set deterministic by deterministic.Q.function(..., called.from.estimate.g=TRUE) must imply that the row is also set deterministic by deterministic.Q.function(..., called.from.estimate.g=FALSE)") #det.Q.fun(T) should imply det.Q.fun(F)
      }

      ua[ua] <- !is.deterministic[ua]
      if (anyNA(ua)) stop("internal ltmle error - ua should not be NA in CleanData")
      if (! all(is.na.strict(data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y), drop=FALSE]))) {
        data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y)] <- NA #if deterministic, set all nodes except Y to NA
        changed <- TRUE
      }

      if (i %in% nodes$C) {
        censored <- data[, i] == "censored" & ua
        if (! all(is.na.strict(data[censored, (i+1):ncol(data), drop=FALSE]))) {
          data[censored, (i+1):ncol(data)] <- NA  #if censored, set all nodes (including Y) to NA
          changed <- TRUE
        }
        ua[ua] <- !censored[ua]
        if (anyNA(ua)) stop("internal ltmle error - ua should not be NA in CleanData")
      }
    }

    if (changed && showMessage) {
      message("Note: for internal purposes, all nodes after a censoring event are set to NA and \n all nodes (except Ynodes) are set to NA after Y=1 if survivalFunction is TRUE (or if specified by deterministic.Q.function).\n Your data did not conform and has been adjusted. This may be relevant if you are \n writing your own deterministic function(s) or debugging ltmle.")
    }
    return(data)
  }

  TransformOutcomes <- function(data, nodes, Yrange) {
    all.Y <- unlist(data[, nodes$Y])
    transformOutcome <- FALSE
    if (!is.null(Yrange)) {
      #if Yrange was specified
      rng <- range(all.Y, na.rm=TRUE)
      if (min(rng) < min(Yrange) || max(rng) > max(Yrange)) {
        #Truncate if Y vals are outside Yrange
        message("Some Ynodes are not in [Yrange[1], Yrange[2]], Y values are truncated")
        data[,nodes$Y][data[,nodes$Y] < min(Yrange)]<- min(Yrange)
        data[,nodes$Y][data[,nodes$Y] > max(Yrange)] <- max(Yrange)
      }
      #Then transform
      transformOutcome <- TRUE
    } else {
      #if Yrange was not specified, get it
      Yrange <- range(all.Y, na.rm=TRUE)
      if (min(Yrange) < 0 || max(Yrange) > 1) {
        #And see if we need to transform
        transformOutcome <- TRUE
        message("Some Ynodes are not in [0, 1], and Yrange was NULL, so all Y nodes are being\ntransformed to (Y-min.of.all.Ys)/range.of.all.Ys")
      }
    }
    if (transformOutcome) {
      attr(transformOutcome, 'Yrange') <- Yrange
      data[,nodes$Y] <- (data[, nodes$Y]-min(Yrange))/diff(Yrange)
    }
    return(list(data=data, transformOutcome=transformOutcome))
  }

  # Converts a general formula to a main terms formula and combine summary measures with baseline covariates
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  ConvertToMainTerms <- function(data, msm, summary.measures, nodes) {
    baseline.column.names <- names(data)[nodes$baseline]
    summary.column.names <- colnames(summary.measures)
    rhs.vars <- RhsVars(msm)
    if (length(intersect(baseline.column.names, summary.column.names)) > 0) stop("Baseline covariate columns of data and columns of summary.measures may not have the same name")
    if (!all(rhs.vars %in% c(baseline.column.names, summary.column.names))) stop("All right hand side variables in working.msm must be either column names of summary.measures or column names of baseline covariates")
    baseline.column.names <- intersect(baseline.column.names, rhs.vars)
    baseline.data <- data[, baseline.column.names, drop=FALSE]
    num.regimes <- dim(summary.measures)[1]
    num.summary.measures <- dim(summary.measures)[2]
    num.final.Ynodes <- dim(summary.measures)[3]
    n <- nrow(data)
    for (j in 1:num.final.Ynodes) {
      for (i in 1:num.regimes) {
        combined.summary.measures <- model.matrix(as.formula(msm), data.frame(Y=1, baseline.data, matrix(summary.measures[i, , j], nrow=n, ncol=num.summary.measures, byrow=TRUE, dimnames=list(NULL, colnames(summary.measures)))))
        if (i == 1 && j == 1) {
          #initialize here now that we know how many columns there are
          main.terms.summary.measures <- array(dim=c(n, ncol(combined.summary.measures), num.regimes, num.final.Ynodes))
          beta.names <- colnames(combined.summary.measures) #this is the same for all i and j
        }
        main.terms.summary.measures[, , i, j] <- combined.summary.measures
      }
    }
    colnames(main.terms.summary.measures) <- paste("S", 1:ncol(main.terms.summary.measures), sep="") #temp names
    main.terms.msm <- paste("Y ~ -1 +", paste(colnames(main.terms.summary.measures), collapse=" + ")) #formula using temp names
    return(list(msm=main.terms.msm, summary.measures=main.terms.summary.measures, beta.names=beta.names, baseline.column.names=baseline.column.names))
  }

  # Calc logical array - n x num.Anodes x num.regimes = follows regime[j] up to and including Anode[i]
  CalcInterventionMatchArray <- function(data, regimes, Anodes) {
    num.regimes <- dim(regimes)[3]
    intervention.match <- array(dim=c(nrow(data), length(Anodes), num.regimes))
    cum.intervention.match <- matrix(TRUE, nrow(data), num.regimes)
    for (Anode.index in seq_along(Anodes)) {
      cum.intervention.match <- cum.intervention.match & ((data[, Anodes[Anode.index]] == regimes[, Anode.index, ]) %in% c(TRUE, NA)) #recycles regimes
      intervention.match[, Anode.index, ] <- cum.intervention.match
    }
    return(intervention.match)
  }

  # Convert censoring nodes stored as binaries into factors (factors are recommended but binaries are currently accepted)
  ConvertCensoringNodes <- function(data, Cnodes, has.deterministic.functions=FALSE) {
    error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored'\n See ?BinaryToCensoring \n (binary is also accepted, where 0=censored, 1=uncensored, but this is not recommended)"
    for (i in Cnodes) {
      col <- data[, i]
      if (is.numeric(col)) {
        if (! all(col %in% c(0, 1, NA))) stop(error.msg)
        data[, i] <- BinaryToCensoring(is.uncensored=col)
        if (has.deterministic.functions) warning("Censoring nodes have been converted from binaries to factors - see ?BinaryToCensoring.\n Note that if you are writing your own deterministic.g.function or deterministic.Q.function that censoring nodes are converted to factors\n before these functions are called.")
      } else if (is.factor(col)) {
        if (! all(levels(col) %in% c("censored", "uncensored"))) {
          stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)")
        }
        #no action required
      } else {
        stop(error.msg)
      }
    }
    return(data)
  }

  ###################################################################################################
  # End of helper functions.
  ###################################################################################################

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
