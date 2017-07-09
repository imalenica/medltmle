################################
# MainCalcsMediation
################################

#' MainCalcsMediation
#'
#' Main TMLE calculations.
#'
#' @param inputs Output of \code{CreateMedInputs}
#'
#' @return Returns
#'

# Loop over final Ynodes, run main calculations...
MainCalcsMediation <- function(inputs) {

  # not doing Estimate Var in this implementation. TO DO.

  #Change to Influence Curve Variance only:
  inputs$IC.variance.only <- T
  if (!exists("est.var.iptw")) est.var.iptw <<- F

  #Check how many final nodes
  num.final.Ynodes <- length(inputs$final.Ynodes)

  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes
  #note: num.measures is summary measures and baseline covariates, converted to main terms
  num.betas <- dim(inputs$combined.summary.measures)[2]

  #Sample size
  n <- nrow(inputs$data)

  #Regime dimension (n x end.time x num.regimes)
  num.regimes <- dim(inputs$regimes)[3]

  #Prep for updates:
  #n x num.regimes x num.final.Ynodes for Qstar
  Qstar <- array(dim=c(n, num.regimes, num.final.Ynodes))
  all.msm.weights <- GetMsmWeights(inputs)
  new.var.y <- array(dim=c(num.betas, num.betas, num.final.Ynodes))
  IC <- matrix(0, n, num.betas)

  #store IC for each final Ynode, compare var(IC) to sum(var(IC.ynode))
  IC.y <- array(dim=c(n, num.betas, num.final.Ynodes))

  #Consitional density estimates using GLM or SuperLearner.

  #Contains prob.A.is.1 (g); individual fits for each node estimated; cumulative g bounded and unbounded.
  #If specified, also cumulative g with mean L (bounded and unbounded)
  g.abar.list <- EstimateG(inputs, regimes.use =  'regimes')
  #Estimate p_z(z|past)
  cum.qz.abar <- EstimateMultiDens(inputs,use.regimes='regimes',use.intervention.match = 'intervention.match',is.Z.dens = T)
  #Estimate p_l(l|past)
  cum.qL.abar <- EstimateMultiDens(inputs,use.regimes='regimes',use.intervention.match = 'intervention.match',is.Z.dens = F)

  if(!setequal(inputs$regimes,inputs$regimes.prime)){

    g.abar.prime.list <- EstimateG(inputs,regimes.use =  'regimes.prime')
    cum.qz.abar.prime <- EstimateMultiDens(inputs,use.regimes='regimes.prime',use.intervention.match = 'intervention.match.prime',is.Z.dens = T)
    cum.qL.abar.prime <- EstimateMultiDens(inputs,use.regimes='regimes.prime',use.intervention.match = 'intervention.match.prime',is.Z.dens = F)

    }else{

    g.abar.prime.list <- g.abar.list
    cum.qz.abar.prime <- cum.qz.abar
    cum.qL.abar.prime <- cum.qL.abar

    }

  #Calculate IPTW estimate
  iptw <- CalcIPTWMediation(inputs, cum.g.abar =  g.abar.list$cum.g, cum.qz.abar = cum.qz.abar$cum.g,cum.qz.abar.prime = cum.qz.abar.prime$cum.g, all.msm.weights)

  # remove any nodes after final.Ynode
  SubsetNodes <- function(nodes, final.Ynode) {
    return(lapply(nodes, function (x) x[x <= final.Ynode]))
  }

  if (inputs$iptw.only) {

    beta <- rep(NA, length(iptw$beta))
    fitted.msm <- NULL
    variance.estimate <- NULL
    fixed.tmle <- NULL

  } else {

    for (j in 1:num.final.Ynodes) {

      print(Sys.time())
      #Update the initial estimate of all the nodes, up to the last one. Get IC and last Qstar. (mean of which is the TMLE).
      #If updating did not happen for some reason, will return gcomp.
      fixed.tmle <- FixedTimeTMLEMediation(inputs, nodes = SubsetNodes(inputs$all.nodes, final.Ynode=inputs$final.Ynodes[j]), msm.weights = drop3(all.msm.weights[, , j, drop=FALSE]), combined.summary.measures = dropn(inputs$combined.summary.measures[, , , j, drop=FALSE], n=4), g.abar.list = g.abar.list, g.abar.prime.list=g.abar.prime.list, cum.qz.abar=cum.qz.abar, cum.qz.abar.prime=cum.qz.abar.prime, cum.qL.abar=cum.qL.abar, cum.qL.abar.prime=cum.qL.abar.prime)
      #If there are multiple final nodes, will sum up the ICs.
      IC <- IC + fixed.tmle$IC
      #ICs for each final Y node seperately.
      IC.y[, , j] <- fixed.tmle$IC
      #Qstar for each final Y node separately.
      Qstar[, , j] <- fixed.tmle$Qstar
      #Non-IC estimated variance.
      new.var.y[, , j] <- fixed.tmle$est.var

    }

    #If user specified that variance should be estimated (non-IC), return IC with all NA
    if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) {
      return(list(IC=matrix(NA, 1, 1), msm=NULL, beta=qlogis(mean(Qstar)), cum.g=g.list$cum.g, cum.g.unbounded=g.list$cum.g.unbounded, fit=fixed.tmle$fit, variance.estimate=NULL, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar))
    }

    #Returns coefficients for the MSM model and predicted values for Qstar. Used for FinalizeIC() and NormalizeIC()
    fitted.msm <- FitPooledMSM(working.msm=inputs$working.msm, Qstar, combined.summary.measures=inputs$combined.summary.measures, msm.weights=all.msm.weights * inputs$observation.weights)
    #Dim: n x num.betas
    IC <- FinalizeIC(IC, inputs$combined.summary.measures, Qstar, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights)
    #C without using g.ratio
    C.old <- NormalizeIC(IC, inputs$combined.summary.measures, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights, g.ratio = NULL)

    if (inputs$IC.variance.only) {

      variance.estimate <- NULL

    } else {
      new.var <- matrix(NA, num.betas, num.betas)
      for (i in 1:num.betas) {
        for (j in 1:num.betas) {
          if (num.final.Ynodes > 1) {
            cov.IC <- cov(IC.y[, i, ], IC.y[, j, ])
            diag(cov.IC) <- new.var.y[i, j, ]
            new.var[i, j] <- sum(cov.IC)
          } else {
            new.var[i, j] <- new.var.y[i, j, 1]
          }
        }
      }

      g.ratio <- CalcGUnboundedToBoundedRatio(g.list, inputs$all.nodes, inputs$final.Ynodes)
      C <- NormalizeIC(IC, inputs$combined.summary.measures, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights, g.ratio)
      variance.estimate <- safe.solve(C) %*% new.var %*% t(safe.solve(C))

    }

    #IC %*% safe.solve(C)
    IC <- t(safe.solve(C.old, t(IC)))
    beta <- coef(fitted.msm$m)
    names(beta) <- inputs$beta.names

  }

  #note: only returns cum.g and fit for the last final.Ynode
  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.gq.abar=list(cum.g=g.abar.list$cum.g, cum.qL=cum.qL.abar$cum.g.block, cum.qZ=cum.qz.abar$cum.g),
              cum.gq.abar.prime=list(cum.g=g.abar.prime.list$cum.g, cum.qL=cum.qL.abar.prime$cum.g.block,cum.qZ=cum.qz.abar.prime$cum.g),
              cum.gq.abar.unbounded=list(cum.g=g.abar.list$cum.g.unbounded,cum.qL=cum.qL.abar$cum.g.unbounded.block,cum.qZ=cum.qz.abar$cum.g.unbounded),
              cum.gq.abar.prime.unbounded=list(cum.g=g.abar.prime.list$cum.g.unbounded,cum.qL=cum.qL.abar.prime$cum.g.unbounded.block,cum.qZ=cum.qz.abar.prime$cum.g.unbounded),
              fit=fixed.tmle$fit, variance.estimate=variance.estimate, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar))
}

################################
# GetMsmWeights
################################

#' GetMsmWeights
#'
#' Set weights for the Marginal Structural Model
#'
#' @param inputs Output of \code{CreateMedInputs}.
#'
#' @return Returns weights for the Marginal Structural Model.
#'

GetMsmWeights <- function(inputs) {

  #Number of samples
  n <- nrow(inputs$data)
  #Number of regimes
  num.regimes <- dim(inputs$regimes)[3]

  stopifnot(num.regimes >= 1)

  #Number of final Y nodes
  num.final.Ynodes <- length(inputs$final.Ynodes)

  #inputs$msm.weights is set by CreateMedMSMInputs()
  if (identical(inputs$msm.weights, "empirical")) {

    #default is probability of following abar given alive, uncensored;
    #conditioning on past treatment/no censoring, but not L, W; duplicates get weight 0

    msm.weights <- matrix(nrow=num.regimes, ncol=num.final.Ynodes)

    if (dim(inputs$regimes)[2] > 0) {
      is.duplicate <- duplicated(inputs$regimes, MARGIN=3)
    } else {
      is.duplicate <- c(FALSE, rep(TRUE, num.regimes - 1))  #in case there are C nodes but no A nodes before a Ynode
    }
    for (j in 1:num.final.Ynodes) {
      final.Ynode <- inputs$final.Ynodes[j]
      uncensored <- IsUncensored(inputs$uncensored, inputs$all.nodes$C, cur.node=final.Ynode)
      intervention.match <- InterventionMatch(inputs$intervention.match, inputs$all.nodes$A, cur.node=final.Ynode)
      for (i in 1:num.regimes) {
        if (is.duplicate[i]) {
          msm.weights[i, j] <- 0
        } else {
          msm.weights[i, j] <- sum(uncensored & intervention.match[, i]) / nrow(inputs$data)
        }
      }
    }
  } else if (is.null(inputs$msm.weights)) {

    msm.weights <- array(1, dim=c(n, num.regimes, num.final.Ynodes))

  } else {

    msm.weights <- inputs$msm.weights

  }

  if (identical(dim(msm.weights), c(num.regimes, num.final.Ynodes))) {

    msm.weights <- array(rep(msm.weights, each=n), dim=c(n, num.regimes, num.final.Ynodes))

  } else if (!identical(dim(msm.weights), c(n, num.regimes, num.final.Ynodes))) {

    stop("dim(msm.weights) should be c(n, num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes)")

  }

  if (anyNA(msm.weights) || any(msm.weights < 0)) stop("all msm.weights must be >= 0 and not NA")
  return(msm.weights)
}

################################
# EstimateG
################################

#' EstimateG
#'
#' Estimation of conditional densities of A nodes (g)
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param regimes.use
#'
#' @return Returns estimate of conditional density for each A nodes, bounded and unbounded cumulative g.
#'

EstimateG <- function(inputs,regimes.use) {

  #Which regime to use
  inputs$regimes <- inputs[[regimes.use]]
  #Number of samples
  n <- nrow(inputs$data)
  #Number of regimes
  num.regimes <- dim(inputs$regimes)[3]
  #All nodes
  nodes <- inputs$all.nodes

  #For each of the regimes, have g for each sample and for each C and A node.
  g <- cum.g <- cum.g.unbounded <- prob.A.is.1 <- array(NaN, dim=c(n, length(nodes$AC), num.regimes))

  #Supports the option of estimated variance later on
  if (inputs$IC.variance.only) {

    cum.g.meanL <- cum.g.meanL.unbounded <- NULL

  } else {

    g.meanL <- cum.g.meanL <- cum.g.meanL.unbounded <- array(NaN, dim=c(n, length(nodes$AC), num.regimes, length(nodes$LY)-1))

  }

  fit <- vector("list", length(nodes$AC))
  names(fit) <- names(inputs$data)[nodes$AC]

  if (!inputs$IC.variance.only && anyNA(inputs$regimes)) {

    regimes.meanL <- inputs$regimes

    #Replace NAs in input$regimes by Mode value
    for (i in nodes$A) {
      for (regime.index in 1:num.regimes) {
        regimes.meanL[is.na(regimes.meanL[, i, regime.index]), i, regime.index] <- Mode(inputs$regimes[, i, regime.index], na.rm = TRUE)
      }
    }
  } else {

    regimes.meanL <- NULL

  }

  #Estimates each node separately.
  for (i in 1:length(nodes$AC)) {

    cur.node <- nodes$AC[i]
    #Returns inputs$uncensored column that corresponds to censoring at the time point being considered. (looks at the last C node..)
    #uncensored will only depend on C and previous C. Does not take into account Y (death, for example).
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)
    #deterministic due to death or Q.function (now looks at Y)
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic

    if (is.numeric(inputs$gform)) {

      if (!inputs$IC.variance.only) stop("IC.variance.only=FALSE not currently compatible with numeric gform")
      if (!is.null(inputs$deterministic.g.function)) stop("deterministic.g.function is not compatible with numeric gform")

      prob.A.is.1[, i, ] <- inputs$gform[, i, ]
      g.est <- list(is.deterministic = deterministic.origdata) #note: this assumes that deterministic.Q.function doesn't depend on A (throw warning in CheckInputs)
      fit[[i]] <- "no fit due to numeric gform"

    } else {

      #Previously specified form of g
      form <- inputs$gform[i]
      #deterministic due to ACnode map - using original data; now considering g
      #Determines which patients have an Anode value which is deterministic. For example, might need to stay on specific treatment once assigned.
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, nodes, using.newdata=F)
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic

      #If stratify, use samples that are uncensored, match intervention, not deterministic for both data and g.
      if (inputs$stratify) {

        intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node=nodes$AC[i])
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata

      } else {

        #Otherwise estimate from uncensored samples with no deterministic data and g.
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata

      }

      #assume all regimes have positive weight for some final.Ynode
      #Will return estimated values for each sample, fit, and which samples are deterministic at the current node (no estimation there).
      g.est <- Estimate(inputs, form=form, Qstar.kplus1=NULL, subs=subs, family=quasibinomial(), type="response", nodes=nodes, called.from.estimate.g=TRUE, calc.meanL=!inputs$IC.variance.only, cur.node=cur.node, regimes.meanL=regimes.meanL, regimes.with.positive.weight=1:num.regimes)
      prob.A.is.1[, i, ] <- g.est$predicted.values
      fit[[i]] <- g.est$fit

    }

    #prob.A.is.1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    if (cur.node %in% nodes$A) {

      #Regime for current A node
      cur.abar <- AsMatrix(inputs$regimes[, nodes$A == cur.node, ])

      if (is.null(regimes.meanL)) {

        cur.abar.meanL <- cur.abar

      } else {

        cur.abar.meanL <- AsMatrix(regimes.meanL[, nodes$A == cur.node, ])

      }

    } else {

      #if this is a cnode, abar is always 1 (uncensored)
      cur.abar <- cur.abar.meanL <- matrix(1, nrow(inputs$data), num.regimes)
    }

    #Recall, no estimate for dead samples, but yes for censored.
    g[, i, ] <- CalcG(AsMatrix(prob.A.is.1[, i, ]), cur.abar, g.est$is.deterministic)

    if (!inputs$IC.variance.only) {
      for (j in sseq(1, dim(g.meanL)[4])) {
        g.meanL[, i, , j] <- CalcG(AsMatrix(g.est$prob.A.is.1.meanL[, , j]), cur.abar.meanL, g.est$is.deterministic)
      }
    }

    if (anyNA(g[uncensored, i, ])) stop("Error - NA in g. g should only be NA after censoring. If you passed numeric gform, make sure there are no NA values except after censoring. Otherwise something has gone wrong.")

  }

  for (regime.index in 1:num.regimes) {
    #Calculates cumulative, bounded g (multiply each estimate)
    cum.g.list <- CalcCumG(AsMatrix(g[, , regime.index]), inputs$gbounds)

    cum.g[, , regime.index] <- cum.g.list$bounded
    cum.g.unbounded[, , regime.index] <- cum.g.list$unbounded

    if (!inputs$IC.variance.only) {
      for (j in sseq(1, dim(g.meanL)[4])) {
        cum.g.list <- CalcCumG(AsMatrix(g.meanL[, , regime.index, j]), inputs$gbounds)
        cum.g.meanL[, , regime.index, j] <- cum.g.list$bounded
        cum.g.meanL.unbounded[, , regime.index, j] <- cum.g.list$unbounded
      }
    }
  }
  return(list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, cum.g.meanL=cum.g.meanL, fit=fit, prob.A.is.1=prob.A.is.1, cum.g.meanL.unbounded=cum.g.meanL.unbounded))
}

################################
# Estimate
################################

#' Estimate
#'
#' Run GLM or SuperLearner to obtain an estimate for the current node.
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param form Q form for current node.
#' @param subs Subset of samples for stratify option.
#' @param family TO DO
#' @param type Two options, "response" or "link".
#' @param nodes TO DO
#' @param Qstar.kplus1 TO DO
#' @param cur.node TO DO
#' @param calc.meanL Defaults to \code{!IC.variance.only}.
#' @param called.from.estimate.g TO DO
#' @param regimes.meanL TO DO
#' @param regimes.with.positive.weight Defaults to \code{1:num.regimes}.
#'
#' @return Returns predicted values for the fit.
#'

Estimate <- function(inputs, form, subs, family, type, nodes, Qstar.kplus1, cur.node, calc.meanL, called.from.estimate.g, regimes.meanL, regimes.with.positive.weight) {

  #Fit and predict using GLM or SuperLearner.
  FitAndPredict <- function() {

    #Check how many samples are left for estimation:
    if (length(Y.subset) < 2) stop("Estimation failed because there are fewer than 2 observations to fit")

    if (use.glm) {

      ##############################
      #estimate using GLM
      ##############################

      SuppressGivenWarnings({

        #Ex:
        #Regress Qstar.kplus1 (Y.subset) on observed past (X.subset) among uncensored/alive samples at time t.
        m <- speedglm.wfit(Y.subset, X.subset, family=family, maxit = 100, weights=observation.weights.subset, offset=offst, intercept=intercept)
        m$terms <- tf
        class(m) <- c("speedglm", "speedlm")

        #Evaluate the fitted function at the observed mediator and covariates histories and the intervened exposure = estimate of Q.k
        predicted.values <- predict(m, newdata=newdata, type=type)

      }, GetWarningsToSuppress())

    } else {

      ##############################
      #estimate using SuperLearner
      ##############################

      #remove aliased (linearly dependent) columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
      #rhs <- setdiff(RhsVars(form), rownames(alias(form, data=X.subset)$Complete))

      newX.list <- GetNewX(newdata)
      SetSeedIfRegressionTesting(inputs)

      try.result <- try({

        #Use other algorithms to evaluate ..E(Qstar.kplus1|X.subset)
        SuppressGivenWarnings(m <- SuperLearner::SuperLearner(Y=Y.subset, X=X.subset, SL.library=SL.library, verbose=FALSE, family=family, newX=newX.list$newX, obsWeights=observation.weights.subset), c("non-integer #successes in a binomial glm!", "prediction from a rank-deficient fit may be misleading"))
      })

      #Evaluate the fitted function at observed mediator and covariates histories at the intervened exposure => Q.k
      predicted.values <- ProcessSLPrediction(m$SL.predict, newX.list$new.subs, try.result)

    }

    return(list(m = m, predicted.values = predicted.values))

  }

  GetSLStopMsg <- function(Y) {
    ifelse(all(Y %in% c(0, 1, NA)), "", "\n Note that many SuperLeaner libraries crash when called with continuous dependent variables, as in the case of initial Q regressions with continuous Y or subsequent Q regressions even if Y is binary.")
  }

  ProcessSLPrediction <- function(pred, new.subs, try.result) {

    if (inherits(try.result, "try-error")) {
      stop(paste("\n\nError occured during call to SuperLearner:\n", form, GetSLStopMsg(Y.subset), "\n The error reported is:\n", try.result))
    }

    if (all(is.na(pred))) {
      stop(paste("\n\nSuperLearner returned all NAs during regression:\n", form, GetSLStopMsg(Y.subset)))
    }

    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- pred

    if (max(predicted.values, na.rm=T) > 1 || min(predicted.values, na.rm=T) < 0) {
      msg <- paste("SuperLearner returned predicted.values > 1 or < 0: [min, max] = [", min(predicted.values, na.rm=T), ",", max(predicted.values, na.rm=T), "]. Bounding to [0,1]")
      warning(msg)
      predicted.values <- Bound(predicted.values, bounds=c(0, 1))
    }

    return(ValuesByType(predicted.values))
  }

  PredictOnly <- function(newdata1) {
    #Prediction with GLM
    if (use.glm) {
      predict(m, newdata1, type)
    } else {
      #Prediction with SuperLearner
      newX.list <- GetNewX(newdata1)
      ProcessSLPrediction(predict(m, newX.list$newX, X.subset, Y.subset, onlySL = TRUE)$pred, newX.list$new.subs, try.result=NULL)
    }
  }

  ValuesByType <- function(x) {
    if (type == "link") {
      stopifnot(family$family %in% c("binomial", "quasibinomial"))
      qlogis(Bound(x, bounds=c(0.0001, 0.9999)))
    } else {
      x
    }
  }

  GetNewX <- function(newdata1) {
    new.mod.frame <- model.frame(f, data = newdata1, drop.unused.levels = TRUE, na.action = na.pass)
    newX.temp <- model.matrix(terms(f), new.mod.frame)
    new.subs <- !rowAnyMissings(newX.temp) #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    newX <- as.data.frame(newX.temp[new.subs, , drop=FALSE])
    if (ncol(X) == 1) { #fixme - prob not needed, intercept will be added unless -1 in form, could check for this in ProcessSLPred
      #SuperLearner crashes if there are screening algorithms and only one column - add a constant
      X.subset <<- cbind(X.subset, ltmle.added.constant=1)
      newX <- cbind(newX, ltmle.added.constant=1)
    }
    return(list(newX=newX, new.subs=new.subs))
  }

  #Predict the probability that A=1 if L and Y nodes are set to their mean (or median) values.

  #probAis1.meanL is n x num.LYnodes - 1
  #probAis1.meanL[, k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L

  #somewhat inefficient - for W A.1 L.2 A.2 L.3 A.3 Y, does P(A.1=1) setting L.3 to mean and then L.2 and L.3 to mean, but none of these can be used in P(A.1=1) because they're after A.1

  PredictProbAMeanL <- function() {

    #A is already set to abar in the data!
    probAis1.meanL <- matrix(NaN, nrow(inputs$data), length(nodes$LY) - 1)

    if (ncol(probAis1.meanL) == 0) return(probAis1.meanL)

    #not the same as nodes$LY, which removes blocks
    all.LY.nodes <- sort(union(nodes$L, nodes$Y))
    LYindex <- length(nodes$LY)

    for (i in length(all.LY.nodes):1) {
      regression.node <- all.LY.nodes[i]
      L <- data[single.subs, regression.node]
      if (is.numeric(L) && !IsBinary(L)) {
        meanL <- mean(L, na.rm = TRUE)
      } else {
        meanL <- Mode(L, na.rm = TRUE) #for factors and binaries
      }

      newdata.meanL[, regression.node] <- meanL

      if (regression.node %in% nodes$LY[1:length(nodes$LY)-1]) {
        LYindex <- LYindex - 1
        probAis1.meanL[, LYindex] <- PredictOnly(newdata = newdata.meanL)
      }
    }
    if (anyNA(probAis1.meanL[, 1])) stop("NA in probAis1.meanL[, 1]")
    return(probAis1.meanL)
  }

  stopifnot(type %in% c("link", "response"))

  num.regimes <- dim(inputs$regimes)[3]

  if (form == "IDENTITY") {

    stopifnot(is.vector(Qstar.kplus1) == 1)
    predicted.values <- ValuesByType(matrix(Qstar.kplus1, nrow = nrow(inputs$data), ncol = num.regimes))
    fit <- as.list(rep("no fit because form == IDENTITY", num.regimes))
    deterministic.list.olddata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g, inputs$survivalOutcome)
    is.deterministic <- matrix(deterministic.list.olddata$is.deterministic, nrow=nrow(inputs$data), ncol=num.regimes)
    deterministic.Q <- matrix(NA, nrow(inputs$data), num.regimes)
    deterministic.Q[is.deterministic, ] <- deterministic.list.olddata$Q
    return(list(predicted.values=predicted.values, fit=fit, is.deterministic=is.deterministic, deterministic.Q=deterministic.Q, prob.A.is.1.meanL=NULL))

  }

  #convert factors to binaries for compatability with glm and some SL libraries
  data <- ConvertCensoringNodesToBinary(inputs$data, nodes$C)
  f <- as.formula(form)

  #Set SL library
  SL.library <- if (called.from.estimate.g) inputs$SL.library.g else inputs$SL.library.Q

  #in a formula like "Y ~ 1", call glm
  use.glm <- (is.null(SL.library) || length(RhsVars(f)) == 0)

  #scale Lnodes to 0-1 to avoid numerical problems in speedglm.
  if (use.glm) {
    for (L in c(nodes$baseline, nodes$L, nodes$Y)) {
      if (is.numeric(data[, L])) {
        mx <- max(abs(data[, L]), na.rm = T)
        if (mx == 0) {
          data[, L] <- 1
        } else if (mx < 0.1 || mx > 10) {
          data[, L] <- data[, L] / mx
        }
      }
    }
  }

  first.regime <- min(regimes.with.positive.weight)

  #One "up" in iterative expectations; E(Qstar.kplus1|past)=Qstar.k
  if (is.null(Qstar.kplus1)) {

    data.with.Qstar <- data

  } else {

    if (is.matrix(Qstar.kplus1)) {

      data.with.Qstar <- cbind(data, Q.kplus1=Qstar.kplus1[, first.regime])

    } else {

      data.with.Qstar <- cbind(data, Q.kplus1=Qstar.kplus1)

    }
  }

  #Set up the model for estimation
  mod.frame <- model.frame(f, data = data.with.Qstar, drop.unused.levels = TRUE, na.action = na.pass)
  #Q.kplus1
  Y <- mod.frame[[1]]
  tf <- terms(f)
  #Intercept and covariates
  X <- model.matrix(tf, mod.frame)
  offst <- model.offset(mod.frame)
  intercept <- attributes(tf)$intercept

  #SL does not support quasibinomial(), change to binomial().
  if (!use.glm) {
    if (identical(family, quasibinomial())) family <- binomial()
    if (!is.null(offst)) stop("offset in formula not supported with SuperLearner")
    X <- as.data.frame(X)
  }

  fit <- vector("list", num.regimes)
  predicted.values <- deterministic.Q <- matrix(NA, nrow(data), num.regimes)
  is.deterministic <- matrix(FALSE, nrow(data), num.regimes)
  Qstar.index <- subs.index <- 1
  fit.and.predict <- NULL
  multiple.subs <- is.matrix(subs)
  multiple.Qstar <- is.matrix(Qstar.kplus1)

  #Probability of A given mean values of L nodes
  if (calc.meanL) {
    prob.A.is.1.meanL <- array(NaN, dim=c(nrow(inputs$data), num.regimes, length(nodes$LY)-1))
    Anode.index <- which(nodes$A < cur.node)
  } else {
    prob.A.is.1.meanL <- NULL
  }

  for (regime.index in regimes.with.positive.weight) {

    #Returns data with Anodes set to regime. If estimating A, will not change.
    newdata <- SetA(data = data.with.Qstar, regimes = inputs$regimes, Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)

    if (calc.meanL) {
      if (!is.null(regimes.meanL)) {
        newdata.meanL <- SetA(data = data.with.Qstar, regimes = regimes.meanL, Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)
      } else {
        newdata.meanL <- newdata
      }
    }

    deterministic.list.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g, inputs$survivalOutcome)

    #Deterministic g function and called from estimate g:
    if (called.from.estimate.g && !is.null(inputs$deterministic.g.function)) {

      newdata.with.current <- newdata
      stopifnot(cur.node %in% nodes$AC)

      if (cur.node %in% nodes$A) {
        #set current node to regime for consistency checking in IsDeterministicG
        newdata.with.current[, cur.node] <- inputs$regimes[, which(nodes$A == cur.node), regime.index]

        } else {
        newdata.with.current <- newdata
      }

      deterministic.g.list.newdata <- IsDeterministicG(newdata.with.current, cur.node, inputs$deterministic.g.function, nodes, using.newdata=T) #deterministic g - using data modified so A = abar

    } else {
      deterministic.g.list.newdata <- list(is.deterministic = rep(FALSE, nrow(data)), prob1 = NULL)
    }

    if (regime.index > first.regime && multiple.Qstar) {
      Y <- Qstar.kplus1[, Qstar.index]
    }

    if (regime.index == first.regime || multiple.subs) {
      #Subset the data
      single.subs <- if (multiple.subs) subs[, subs.index] else subs
      X.subset <- X[single.subs, , drop=FALSE]
      #if there is a column of all zeros, speedglm may crash - replace with column of 1s
      X.subset[, colAlls(X.subset == 0)] <- 1
      #Assign weights for subset of samples
      observation.weights.subset <- inputs$observation.weights[single.subs]
      offst.subset <- offst[single.subs]
    }

    if (regime.index == first.regime || multiple.subs || multiple.Qstar) {
      Y.subset <- Y[single.subs]
      if (anyNA(Y.subset)) stop("NA in Estimate")
    }

    if (is.numeric(form)) {

      #if gform is numeric, it's a matrix of prob.A.is.1
      predicted.values[, regime.index] <- form[, regime.index]
      m <- "gform passed as numeric, no estimation took place"

      #If all deterministic, no point in estimation...
    } else if (!all(deterministic.list.newdata$is.deterministic | deterministic.g.list.newdata$is.deterministic)) {

       if (is.null(fit.and.predict) || multiple.Qstar || multiple.subs) {

        #Here
        fit.and.predict <- FitAndPredict()
        m <- fit.and.predict$m
        predicted.values[, regime.index] <- fit.and.predict$predicted.values

      } else {
        #just predict
        predicted.values[, regime.index] <- PredictOnly(newdata)
      }

      #If calc.meanL was set to TRUE, use PredictProbAMeanL() function.
      if (calc.meanL) prob.A.is.1.meanL[, regime.index, ] <- PredictProbAMeanL()

    } else {
      m <- "all rows are deterministic, no estimation took place"
    }

    #For deterministic samples, assign the deterministic value from deterministic.g.list.newdata$prob1
    predicted.values[deterministic.g.list.newdata$is.deterministic, regime.index] <- deterministic.g.list.newdata$prob1

    if (calc.meanL) prob.A.is.1.meanL[deterministic.g.list.newdata$is.deterministic, regime.index, ] <- deterministic.g.list.newdata$prob1
    is.deterministic[, regime.index] <- deterministic.list.newdata$is.deterministic
    if (!called.from.estimate.g) deterministic.Q[deterministic.list.newdata$is.deterministic, regime.index] <- deterministic.list.newdata$Q
    if (!use.glm && !isTRUE(attr(SL.library, "return.fit", exact = TRUE))) m <- summary(m)
    fit[[regime.index]] <- m
    if (multiple.subs) subs.index <- subs.index + 1
    if (multiple.Qstar) Qstar.index <- Qstar.index + 1
  }

  if (all(is.na(predicted.values))) stop("??")

  return(list(predicted.values=predicted.values, fit=fit, is.deterministic=is.deterministic, deterministic.Q=deterministic.Q, prob.A.is.1.meanL=prob.A.is.1.meanL))
}

################################
# EstimateMultiDens
################################

#' EstimateMultiDens
#'
#' Estimates conditional densities of Z and L.
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param use.regimes TO DO
#' @param use.intervention.match TO DO
#' @param is.Z.dens TO DO
#'
#' @return Returns estimate of Z and L conditional densities given the regime.
#'

EstimateMultiDens <- function(inputs,use.regimes,use.intervention.match,is.Z.dens){

  #Specify regime
  inputs$regimes <- inputs[[use.regimes]]
  #Number of samples
  n <- nrow(inputs$data)
  #Number of regimes
  num.regimes <- dim(inputs[[use.regimes]])[3]
  #Specify nodes
  nodes <- inputs$all.nodes

  #Z of L density?
  if(is.Z.dens){
    dens.nodes <- nodes$Z
    dens.forms <- inputs$qzform
  }else{
    dens.nodes <- sort(c(nodes$L,nodes$Y))
    dens.forms <- inputs$qLform
  }

  fit <- vector('list',length(dens.nodes))
  g <- cum.g <- cum.g.unbounded <- prob.is.1 <- array(NaN, dim=c(n, length(dens.nodes), num.regimes))

  #estimate P(node|A=abar, past)
  for (i in 1:length(dens.nodes)) {

    #Current node being estimated
    cur.node <- dens.nodes[i]
    #Censoring based on the last C node
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)
    #Deterministic based on the last Y node (if survival outcome) or deterministic Q
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic

    if (!is.numeric(dens.forms)) {

      form <- dens.forms[i]
      #deterministic due to ACnode map - using original data
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, nodes, using.newdata=F)
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic

      if (inputs$stratify) {
        intervention.match <- InterventionMatch(inputs[[use.intervention.match]], nodes$A, cur.node=cur.node)
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata
      } else {
        #Subset to samples that are uncensored, not deterministic due to death/Q or g.
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
      }

      #assume all regimes have positive weight for some final.Ynode.
      #Estimates current node with A set to the regime, and form specified.
      g.est <- Estimate(inputs, form=form, Qstar.kplus1=NULL, subs=subs, family=quasibinomial(), type="response", nodes=nodes, called.from.estimate.g=TRUE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=1:num.regimes)
      prob.is.1[, i, ] <- g.est$predicted.values
      fit[[i]] <- g.est$fit

      # probZis1 <- rep(NaN, nrow(d))
      # if (any(subs)) {
      #   est <- Estimate(qzform[i], d=d, subs=subs, family="binomial", newdata=newdata[subs,], SL.library=SL.library)
      #   probZis1[subs] <- est$values
      #   stats <- cbind(stats, est$stats)
      #   colnames(stats)[i] <- names(d)[nodes$Z][i]
      # }
      #
      # probZis1[deterministic.g] <- deterministic.g.list$prob1
    } else {

      if (!inputs$IC.variance.only) stop("IC.variance.only=FALSE not currently compatible with numeric gform")
      if (!is.null(inputs$deterministic.g.function)) stop("deterministic.g.function is not compatible with numeric gform")

      prob.is.1[, i, ] <- dens.forms[, i, ]
      g.est <- list(is.deterministic = deterministic.origdata) #note: this assumes that deterministic.Q.function doesn't depend on A (throw warning in CheckInputs)
      fit[[i]] <- "no fit due to numeric gform"

    }

    #prob.Z.is.1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    cur.vals <- AsMatrix(inputs$data[,cur.node])
    g[, i, ] <- CalcG(AsMatrix(prob.is.1[, i, ]), cur.vals, g.est$is.deterministic)
  }

  for (regime.index in 1:num.regimes) {
    #Cumulative g: bounded and unbounded
    cum.g.list <- CalcCumG(AsMatrix(g[, , regime.index]), inputs$gbounds)
    cum.g[, , regime.index] <- cum.g.list$bounded
    cum.g.unbounded[, , regime.index] <- cum.g.list$unbounded
  }

  ret.list <- list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded,prob.is.1=prob.is.1)

  ## if L node, qL for blocks of L node.
  #For this ex, it will return LA1,Y1,LA2,Y2.
  if(!is.Z.dens){

    max.Lnodes.col <- vector('numeric',length(nodes$LY))

    for (i in 1:length(nodes$LY)) {

      cat(date(),'EstimateLYnodes node ',i,'\n')

      if(i==length(nodes$LY)){
        max.Lnodes.col[i] <- which(dens.nodes==max(dens.nodes))
      }
      else{
        ## the desn.node that is smaller than this one
        max.Lnodes.col[i] <- which(dens.nodes==nodes$LY[i+1])-1
      }
    }

    ret.list$cum.g.block <- cum.g[,max.Lnodes.col,,drop=FALSE]
    ret.list$cum.g.unbounded.block <- cum.g.unbounded[,max.Lnodes.col,,drop=FALSE]
  }
  return(ret.list)
}

################################
# CalcIPTWMediation
################################

#' CalcIPTWMediation
#'
#' Calculates IPTW estimate.
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param cum.g.abar Cumulative g obtained from \code{CalcCumG}  with regime abar.
#' @param cum.qz.abar Cumulative Q for Z obtained from \code{CalcCumG} with regime abar.
#' @param cum.qz.abar.prime Cumulative Q for Z obtained from \code{CalcCumG} with regime abar.prime.
#' @param msm.weights Marginal Structural Model weights for each sample.
#'
#' @return Returns IPTW estimate and normalized influence curve for IPTW.
#'

CalcIPTWMediation <- function(inputs, cum.g.abar, cum.qz.abar, cum.qz.abar.prime,msm.weights) {

  if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) {
    return(list(beta=NA, IC=matrix(NA, 1, 1)))
  }

  #Specify nodes
  nodes <- inputs$all.nodes
  #Number of samples
  n <- nrow(inputs$data)
  #Number of regimes
  num.regimes <- dim(inputs$regimes)[3]
  #NUmber of final Y nodes.
  num.final.Ynodes <- length(inputs$final.Ynodes)

  #Vectors to save all outputs in case of multiple outcomes.
  Y.vec <- X.mat <- weight.vec <- NULL
  save.xy <- list()

  #Estimate weights for each final Y node and each regime.
  for (j in 1:num.final.Ynodes) {

    #Which node is a final Y node.
    final.Ynode <- inputs$final.Ynodes[j]
    #Returns the intervention match for the closest A node.
    intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node=final.Ynode)
    #Check censoring for the closest C node (to Y)
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, final.Ynode)

    #For each regime:
    for (i in 1:num.regimes) {

      #Samples that are uncensored and match intervention for the regime i:
      index <- uncensored & intervention.match[, i]
      #Look for the AC node closest to final Y
      col.index <- which.max(nodes$AC[nodes$AC < final.Ynode])
      #Look for the Z node closest to final Y
      col.z.index <- which.max(nodes$Z[nodes$Z < final.Ynode])

      #Limit to samples matching index (uncensored and intervention match)
      Y <- inputs$data[index, final.Ynode]

      #Return cum.g.bar (p_A(a|past)) for closest node to Y for regime i (for sample ok by index).
      g <- cum.g.abar[index, col.index, i]

      #p_z(z|abar.prime,past)/p_z(z|abar,past)
      qz.ratio <- cum.qz.abar.prime[index, col.z.index, i]/cum.qz.abar[index, col.z.index, i]

      #Dimensions are n (here, samples that are u and fr) x num.betas x num.regimes x num.final.Ynodes
      X <- inputs$combined.summary.measures[index, , i, j]

      #if only one summary.measure or sum(index)==1, X is dropped to vector
      if (is.vector(X)) {
        dim(X) <- c(sum(index), ncol(inputs$combined.summary.measures))
      }

      #Marginal Structural Model weights
      #msm.weights is n x num.regimes x num.final.Ynodes
      weight <- msm.weights[index, i, j] * inputs$observation.weights[index] * qz.ratio / g
      #avoid problems where weight and g are both 0
      weight[msm.weights[index, i, j] == 0 | inputs$observation.weights[index] == 0] <- 0

      #Saves output for each regime and final Y node.
      save.xy[[length(save.xy) + 1]] <- list(X=X, Y=Y, weight=weight, index=index)
      Y.vec <- c(Y.vec, Y)
      X.mat <- rbind(X.mat, X)
      weight.vec <- c(weight.vec, weight)
    }
  }

  colnames(X.mat) <- colnames(inputs$combined.summary.measures)

  #this happens if there are no rows uncensored and intervention.match (no samples to estimate IPTW from)
  if (nrow(X.mat) == 0) {

    warning("no rows uncensored and matching regimes/abar - IPTW returns NA")
    num.beta <- ncol(inputs$combined.summary.measures)
    return(list(beta=rep(NA, num.beta), IC=matrix(nrow=n, ncol=num.beta)))

  }

  #working.msm: Estimate coefficient for S1 (beta) for the simple example.
  #Scale weights- large weights might cause convergence problems
  m.glm <- speedglm(formula(inputs$working.msm), family=quasibinomial(), data=data.frame(Y=Y.vec, X.mat, weight.vec), weights=as.vector(scale(weight.vec, center=FALSE)))
  beta <- coef(m.glm)

  #n x num.betas (number of estimated parameters)
  IC <- matrix(0, nrow=n, ncol=length(beta))
  #n x num.regimes x num.final.Ynodes
  m.beta <- array(dim=c(n, num.regimes, num.final.Ynodes))
  cnt <- 0

  for (j in 1:num.final.Ynodes) {

    #Set the final node.
    final.Ynode <- inputs$final.Ynodes[j]

    for (i in 1:num.regimes) {

      #All samples now
      newdata <- data.frame(inputs$combined.summary.measures[, , i, j])
      colnames(newdata) <- colnames(inputs$combined.summary.measures) #needed if only one summary measure
      SuppressGivenWarnings(m.beta[, i, j] <- predict(m.glm, newdata=newdata, type="response"), "prediction from a rank-deficient fit may be misleading")

      cnt <- cnt + 1
      XY.list <- save.xy[[cnt]]

      #recycles weight, Y, m.beta
      #IC for each uncensored and intervention.match sample
      #Penalize the difference between the Y and its estimate, m.beta
      IC[XY.list$index, ] <- IC[XY.list$index, ] + XY.list$weight * XY.list$X * (XY.list$Y - m.beta[XY.list$index, i, j])

      }
  }

  C <- NormalizeIC(IC, inputs$combined.summary.measures, m.beta, msm.weights, observation.weights=inputs$observation.weights, g.ratio=NULL)
  normalized.IC <- t(safe.solve(C, t(IC)))
  names(beta) <- inputs$beta.names
  return(list(beta=beta, IC=normalized.IC))
}

################################
# FixedTimeTMLEMediation
################################

#' FixedTimeTMLEMediation
#'
#' TMLE Mediation
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param nodes All nodes in the data.
#' @param msm.weights Marginal structural model weights.
#' @param combined.summary.measures TO DO
#' @param g.abar.list g estimate from \code{EstimateG} with regime abar.
#' @param g.abar.prime.list g estimate from \code{EstimateG} with regime abar.prime.
#' @param cum.qz.abar Cumulative estimate for Z nodes under abar regime.
#' @param cum.qz.abar.prime Cumulative estimate for Z nodes under abar.prime regime.
#' @param cum.qL.abar Cumulative estimate for L nodes under abar regime.
#' @param cum.qL.abar.prime Cumulative estimate for L nodes under abar.prime regime.
#'
#' @return Return...
#'

FixedTimeTMLEMediation <- function(inputs, nodes, msm.weights, combined.summary.measures, g.abar.list , g.abar.prime.list, cum.qz.abar, cum.qz.abar.prime, cum.qL.abar, cum.qL.abar.prime){

  #combined.summary.measures: n x num.measures x num.regimes
  #(num.measures=num.summary.measures + num.baseline.covariates)
  LYZnodes <- sort(c(nodes$LY, nodes$Z))
  data <- inputs$data

  #Number of regimes
  num.regimes <- dim(inputs$regimes)[3]
  #Number of samples
  n <- nrow(data)
  #Betas
  num.betas <- ncol(combined.summary.measures)

  #Prep for TMLE and IC
  tmle <- rep(NA, num.regimes)
  IC <- matrix(0, nrow=n, ncol=num.betas)

  est.var <- matrix(0, num.betas, num.betas)
  regimes.with.positive.weight <- which(apply(msm.weights > 0, 2, any))

  if (length(regimes.with.positive.weight) == 0) stop("All regimes have weight 0 (one possible reason is that msm.weights='emipirical' and no data rows match any of the regimes and are uncensored)")

  #Prep for updates
  fit.Qstar <- vector("list", length(LYZnodes))
  names(fit.Qstar) <- names(data)[LYZnodes]
  fit.Q <-  vector("list", length(regimes.with.positive.weight))

  #Set up fits for each LYZ node
  for (i in regimes.with.positive.weight) fit.Q[[i]] <- fit.Qstar

  #Initiate Qstar.kplus1=Y.
  Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow=n, ncol=num.regimes)
  mean.summary.measures <- apply(abs(combined.summary.measures), 2, mean)

  #at and after cur.node
  for (i in length(LYZnodes):1) {

    #Current node
    cur.node <- LYZnodes[i]
    #Looks at the most recent Y
    deterministic.list.origdata <- IsDeterministic(data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
    #Looks at the most recent C
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)

    ## if at L node: UPDAYE Q_L!
    if(cur.node %in% nodes$LY){

      #Only relevant for stratify: samples that match the rule.
      intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node)

      if (inputs$stratify) {

        subs <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic #n x num.regimes

        } else {

        #Uncensored and non-deterministic samples
        subs <- uncensored & !deterministic.list.origdata$is.deterministic

        }

      #Initial estimate of Q.k for the current node.
      #Obtained by estimating E(Qstar.kplus1|past) by either SL or regression.
      #If this is the last node, only pass the first column as a vector
      Q.est <- Estimate(inputs, form = inputs$QLform[which(nodes$LY==cur.node)], Qstar.kplus1=if (i == length(LYZnodes)) Qstar.kplus1[, 1] else Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight)
      #Initial estimate of Q.k for the current node
      logitQ <- Q.est$predicted.values
      #Fit
      fit.Q[[i]] <- Q.est$fit

      #Closest AC and Z node to the current node
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      Znode.index <- which.max(nodes$Z[nodes$Z < cur.node])

      #Update Q.k to get Qstar.k
      update.list <- UpdateQMediation(Qstar.kplus1, logitQ, combined.summary.measures, cum.g = if(length(ACnode.index)==0) 1 else g.abar.list$cum.g[, ACnode.index, ],cum.q.ratio=if(length(Znode.index)==0) 1 else cum.qz.abar.prime$cum.g[,Znode.index,]/cum.qz.abar$cum.g[,Znode.index,], working.msm=inputs$working.msm, uncensored, intervention.match, is.deterministic=deterministic.list.origdata$is.deterministic, msm.weights, gcomp=inputs$gcomp, observation.weights=inputs$observation.weights)
      Qstar <- update.list$Qstar
      #Update Qstar for samples that are deterministic.
      Qstar[Q.est$is.deterministic] <- Q.est$deterministic.Q[Q.est$is.deterministic]

      #Calculate the influence curve and relative error
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) / mean.summary.measures

      #If any IC relative error is large (sum of the IC is not ~0) and we don't want gcomp estimate, fix score equation.
      if (any(curIC.relative.error > 0.001) && !inputs$gcomp) {

        cat("fixing: ", curIC.relative.error, "\n")
        SetSeedIfRegressionTesting()
        #Sometimes GLM does not converge and the updating step of the TMLE does not solve the score equation.
        fix.score.list <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, Q.est$is.deterministic, Q.est$deterministic.Q, update.list$off, update.list$X, regimes.with.positive.weight)
        Qstar <- fix.score.list$Qstar
        #Recalculate the IC with a new Qstar
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
        update.list$fit <- fix.score.list$fit

      }

      Qstar.kplus1 <- Qstar
      fit.Qstar[[i]] <- update.list$fit
    }## done with updating QL!

    #If current node is Z, update Q_Z!
    if(cur.node %in% nodes$Z){

      #Samples for which intervention matches abar.prime. Important for stratify=TRUE
      intervention.match <- InterventionMatch(inputs$intervention.match.prime, nodes$A, cur.node)

      if (inputs$stratify) {

        #Subset to uncensored, alive and exposure matching the rule samples.
        subs <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic #n x num.regimes

      } else {

        #Subset to uncensored and alive samples.
        subs <- uncensored & !deterministic.list.origdata$is.deterministic

      }

      #Q.est <- Estimate(inputs = set(inputs,'regimes',inputs$regimes.prime), form = inputs$QZform[which(nodes$Z==cur.node)], Qstar.kplus1=Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight)

      #Get initial estimate of Q from SL or regressing Qstar.kplus1 (estimate from the previous step) on past.
      #Evaluate the fitted function at the observed mediatior and covariates and the intervened exposure.
      Q.est <- Estimate(inputs, form = inputs$QZform[which(nodes$Z==cur.node)], Qstar.kplus1=Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight)
      logitQ <- Q.est$predicted.values
      fit.Q[[i]] <- Q.est$fit

      #Get the closest AC nodes to the current node
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      #Get the closest L node to the current node.
      Lnode.index <- which.max(nodes$LY[nodes$LY < cur.node])

      #Update QMediation. Clever covariate for Z will need g and cumulative conditional density of L.
      update.list <- UpdateQMediation(Qstar.kplus1, logitQ, combined.summary.measures, cum.g = if(length(ACnode.index)==0) 1 else g.abar.prime.list$cum.g[, ACnode.index, ],cum.q.ratio=if(length(Lnode.index)==0) 1 else cum.qL.abar$cum.g.block[,Lnode.index,]/cum.qL.abar.prime$cum.g.block[,Lnode.index,], inputs$working.msm, uncensored, intervention.match, deterministic.list.origdata$is.deterministic, msm.weights, inputs$gcomp, inputs$observation.weights)
      Qstar <- update.list$Qstar
      #Update Qstar so that deterministic samples are included
      Qstar[Q.est$is.deterministic] <- Q.est$deterministic.Q[Q.est$is.deterministic]

      #Calculate the influence curve and its relative error
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) / mean.summary.measures

      #If TMLE does not solve the score equation, try manually to fix it and update IC calculation.
      if (any(curIC.relative.error > 0.001) && !inputs$gcomp) {

        cat("fixing: ", curIC.relative.error, "\n")
        SetSeedIfRegressionTesting(inputs)
        fix.score.list <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, Q.est$is.deterministic, Q.est$deterministic.Q, update.list$off, update.list$X, regimes.with.positive.weight)
        Qstar <- fix.score.list$Qstar
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
        update.list$fit <- fix.score.list$fit

      }

      Qstar.kplus1 <- Qstar
      fit.Qstar[[i]] <- update.list$fit
    }#done with updating Q_Z.

    ## if at W2 node: UPDAYE Q_W2!




    #TO DO: NON_IC VARIANCE.
    #est.var <- est.var + EstimateVariance(inputs, nodes, combined.summary.measures, regimes.with.positive.weight, uncensored, alive=!deterministic.list.origdata$is.deterministic, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, g.list$cum.g, g.list$prob.A.is.1, g.list$cum.g.meanL, g.list$cum.g.unbounded, g.list$cum.g.meanL.unbounded, inputs$observation.weights, is.last.LYnode=(LYnode.index==length(nodes$LY)), intervention.match) #fixme - remove intervention.match if not using est.var.iptw

    #Final IC will be sum of all separate ICs for each node.
    #Naturally, if we choose to estimate conditional density of Z and W and flunctuate it, this variance estimate based on IC will be larger.
    IC <- IC + curIC

  }


  return(list(IC=IC, Qstar=Qstar, est.var=NA,
              fit=list(g=NULL, Q=fit.Q, Qstar=fit.Qstar)))

}

################################
# CalcG
################################

#' CalcG
#'
#' Calculate G
#'
#' @param prob.A.is.1 TO DO
#' @param cur.abar TO DO
#' @param deterministic.newdata TO DO
#'
#' @return Returns G.
#'

CalcG <- function(prob.A.is.1, cur.abar, deterministic.newdata) {
  g <- matrix(NA_real_, nrow(prob.A.is.1), ncol(prob.A.is.1))
  g[!is.na(cur.abar) & cur.abar == 1] <- prob.A.is.1[!is.na(cur.abar) & cur.abar == 1] #matrix indexing
  g[!is.na(cur.abar) & cur.abar == 0] <- 1 - prob.A.is.1[!is.na(cur.abar) & cur.abar == 0] #matrix indexing
  g[deterministic.newdata] <- 1  #a=abar deterministically after death or other deterministic Q (matrix indexing)
  return(g)
}

################################
# CalcCumG
################################

#' CalcCumG
#'
#' Calculate bounded cumulative G.
#'
#' @param g Estimate of \code{g} derived from \code{EstimateG()}.
#' @param gbounds Lower and upper bounds for g.
#'
#' @return Returns G.
#'

CalcCumG <- function(g, gbounds) {
  cum.g <- rowCumprods(g)
  return(list(unbounded=cum.g, bounded=Bound(cum.g, gbounds)))
}

################################
# NormalizeIC
################################

#' NormalizeIC
#'
#' Normalize the influence curve matrix.
#'
#' @param IC Estimate of the influence curve.
#' @param combined.summary.measures Combined summary measure of (n x num.measures x num.regimes x num.final.Ynodes) dimension. Note that num.measures is equal to num.summary.measures + num.baseline.covariates.
#' @param m.beta Estimate of betas.
#' @param msm.weights Marginal structural model weights.
#' @param observation.weights Observation weights.
#' @param g.ratio \code{g.unbounded} / \code{g.bounded} ratio of dimensions: n x num.regimes x num.final.Ynodes. If \code{IC.variance.only}, g.ratio should be NULL.
#'
#' @return Returns normalized influence curve.
#'

NormalizeIC <- function(IC, combined.summary.measures, m.beta, msm.weights, observation.weights, g.ratio) {

  #Number of samples
  n <- nrow(IC)
  #Number of betas
  num.betas <- ncol(IC)
  #Number of regimes
  num.regimes <- dim(combined.summary.measures)[3]
  #Number of final Y nodes.
  num.final.Ynodes <- dim(combined.summary.measures)[4]

  if (is.null(g.ratio)) {
    #if IC.variance.only, g.ratio should be NULL
    g.ratio <- array(1, dim=c(n, num.regimes, num.final.Ynodes))
  }

  C <- array(0, dim=c(num.betas, num.betas))

  #For each final Y node and each regime:
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {

      tempC <- crossprod(combined.summary.measures[, , i, j] * g.ratio[, i, j], combined.summary.measures[, , i, j] * g.ratio[, i, j] * msm.weights[, i, j] * m.beta[, i, j] * (1 - m.beta[, i, j]) * observation.weights)
      if (anyNA(tempC)) stop("NA in tempC")
      C <- C + tempC

    }
  }

  C <- C / n

  if (F) { #fixme

    C2 <- array(0, dim=c(num.betas, num.betas, n))

    for (j in 1:num.final.Ynodes) {
      for (i in 1:num.regimes) {
        positive.msm.weights <- which(msm.weights[, i, j] > 0)
        tempC <- array(0, dim=c(num.betas, num.betas, n))
        for (k in positive.msm.weights) {
          if (max(m.beta[, i, j]) > (1 - 1e-6) || min(m.beta[, i, j]) < 1e-6) {
            warning("All predicted probabilities are all very close to 0 or 1. Unable to compute standard errors.")
            return(matrix(NA, nrow=num.betas, ncol=num.betas))
          }
          m.beta.temp <- m.beta[k, i, j]
          h <- matrix(combined.summary.measures[k, , i, j], ncol=1) * msm.weights[k, i, j] * g.ratio[k, i, j]

          tempC[, , k] <- h %*% t(h) * m.beta.temp * (1 - m.beta.temp) * observation.weights[k] / msm.weights[k, i, j]
        }
        if (anyNA(tempC)) stop("NA in tempC")
        C2 <- C2 + tempC
      }
    }
    C2 <- apply(C2, c(1, 2), mean)
    if (max(abs(C-C2)) > 0.00001) stop("C and C2 do not match")
  }

  if (rcond(C) < 1e-12) {
    C <- matrix(NA, nrow=num.betas, ncol=num.betas)
    warning("rcond(C) near 0, standard errors not available")
  } else {
    normalized.IC <- t(safe.solve(C, t(IC))) #IC %*% safe.solve(C)
    if (!any(abs(colSums(IC)) > 0.001) && !anyNA(normalized.IC) && any(abs(colSums(normalized.IC)) > 0.001)) {
      msg <- capture.output({
        cat("normalized IC problem", colSums(normalized.IC), "\n")
        cat("inv(C) = \n")
        print(safe.solve(C))
      })
      warning(paste(msg, collapse="\n"))
    }
  }
  return(C)
}

################################
# UpdateQMediation
################################

#' UpdateQMediation
#'
#' Update the initial fit of Q using clever covariates.
#'
#' @param Qstar.kplus1 Estimate of the expectation with respect to the distribution of the node one ahead of the current node given the past (dimension: n x num.regimes).
#' @param logitQ logit predicted values of Q for the current node, as estimated by \code{Estimate} (dimension: n x num.regimes).
#' @param combined.summary.measures TO DO (dimension: n x num.measures x num.regimes, where num.measures=num.summary.measures + num.baseline.covariates).
#' @param cum.g Cumulative g estimate up to the most recent AC node. (dimension: n x num.regimes x num.measures)
#' @param cum.q.ratio Cumulative Q estimate ratio of following abar.prime regime and following abar regime. (dimension: n x num.regimes x num.measures)
#' @param working.msm Working marginal structural model.
#' @param uncensored Uncensored samples.
#' @param intervention.match Samples with exposure that matches intervention (dimension: n x num.regimes).
#' @param is.deterministic Logical variable indicating samples deterministic due to death or deterministic Q.
#' @param msm.weights Marginal structural model weights (dimension: n x num.regimes).
#' @param gcomp Logical indicating whether to use Gcomp instead (no updating if TRUE).
#' @param observation.weights Sample weights.
#'
#' @return Returns updated estimate of Q, unless gcomp=TRUE or no samples left for estimation purposes.
#' The output also includes the offset (initial estimate of Q for the current node), fluctuation model, and the clever covariate.
#'

#Note:
#summary.measures: num.regimes x num.summary.measures
#baseline.covariates: names/indicies: num.baseline.covariates x 1
#combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)

UpdateQMediation <- function(Qstar.kplus1, logitQ, combined.summary.measures, cum.g, cum.q.ratio, working.msm, uncensored, intervention.match, is.deterministic, msm.weights, gcomp, observation.weights) {

  #Number of samples
  n <- nrow(logitQ)
  #Number of regimes
  num.regimes <- ncol(logitQ)
  #Offset
  off <- as.vector(logitQ)
  Y <- as.vector(Qstar.kplus1)

  #stacked.summary.measures: (n*num.regimes) x num.measures
  stacked.summary.measures <- apply(combined.summary.measures, 2, rbind)

  #recycles uncensored and is.deterministic, adds intervention.match
  subs.vec <- uncensored & !is.deterministic & as.vector(intervention.match)
  weight.vec <- numeric(n * num.regimes)
  #Calculate sample weight for subset of samples as specified above (subsetting avoids problems with NA in cum.g)
  #Part of clever covariate:
  weight.vec[subs.vec] <- (observation.weights * as.vector(msm.weights) * as.vector(cum.q.ratio)/ as.vector(cum.g))[subs.vec]

  if (anyNA(weight.vec)) stop("NA in weight.vec")

  #For the first update:
  #\epsilon is the coefficient of a weighted logistic regression of Qstar.kplus1 onto the
  #intercept model (S1) with an offset logitQ (off) and weights weight.vec.
  f <- as.formula(paste(working.msm, "+ offset(off)"))

  #Contains output, intercept S1 and offset (logitQ):
  data.temp <- data.frame(Y, stacked.summary.measures, off)

  if (gcomp) {

    Qstar <- plogis(logitQ)
    m <- "no Qstar fit because gcomp=TRUE (so no updating step)"

  } else {
    if (any(weight.vec > 0)) {

      #m is the model:
      #weighted logistic regression of Y onto the intercept model with an offset logitQ (no need to estimate coeff for it) for the subset of samples.
      #this should include the indicators
      SuppressGivenWarnings(m <- speedglm(f, data=data.temp[weight.vec > 0, ], family=quasibinomial(), weights=as.vector(scale(weight.vec[weight.vec > 0], center=FALSE)), maxit=100), GetWarningsToSuppress(TRUE))

      #UPDATING STEP:
      #this should NOT include the indicators. Estimate for all, with a different intercept now (\epsilon) and offset as before
      SuppressGivenWarnings(Qstar <- matrix(predict(m, newdata=data.temp, type="response"), nrow=nrow(logitQ)), GetWarningsToSuppress(TRUE))

      } else {
      Qstar <- plogis(logitQ)
      m <- "no Qstar fit because no subjects alive, uncensored, following intervention. Returns gcomp estimate."
    }

  }

  #I(A=rule and uncensored) * observation.weights
  #Note: followed rule at some point, if closest A is NA
  indicator <- matrix(uncensored * observation.weights, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * matrix(intervention.match, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures))
  # I() * h * observation.weights * cum.q.ratio / g
  h.g.ratio <- stacked.summary.measures * matrix(cum.q.ratio/cum.g, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * indicator
  dim(h.g.ratio) <- c(n, num.regimes, ncol(h.g.ratio))

  for (i in 1:num.regimes) {

    #Add msm.weights
    h.g.ratio[, i, ] <- h.g.ratio[, i, ] * msm.weights[, i]
    weight.zero.index <- msm.weights[, i] == 0
    #cum.g is 0 so X is NA so h.g.ratio is NA when weight is 0
    h.g.ratio[weight.zero.index, i, ] <- 0

  }

  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=stacked.summary.measures, off=off, fit=m))
}

################################
# CalcIC
################################

#' CalcIC
#'
#' Calculate the TMLE influence curve for one node Q estimate.
#'
#' @param Qstar.kplus1 Estimate of the expectation with respect to the distribution of the node one ahead of the current node given the past (dimension: n x num.regimes).
#' @param Qstar Possibly targeted estimate of Q for the current node.
#' @param h.g.ratio I() * h * observation.weights * cum.q.ratio / g
#' @param uncensored Uncensored samples
#' @param intervention.match Samples that match the rule.
#' @param regimes.with.positive.weight Regimes with positive weight.
#'
#' @return Returns IC for one node Q estimate.
#'

CalcIC <- function(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight) {

  #Number of samples
  n <- nrow(Qstar)
  #Number of regimes
  num.regimes <- ncol(Qstar)
  #Betas: h.g.ratio: n x num.regimes x num.betas
  num.betas <- dim(h.g.ratio)[3]

  #IC for each beta
  IC <- matrix(0, nrow=n, ncol=num.betas)

  for (i in regimes.with.positive.weight) {

    #Samples that are uncensored and receive exposure matching the rule i
    #No IC for samples that are censored or do not match rule (since Qstar is for A that matches the rule- diff between Qstar.kplus1 and Qstar does not make sense otherwise) )
    index <- uncensored & intervention.match[, i]

    #If all h.g.ratios are 0, IC will be 0 anyways.
    #As per theory, h.g.ratio will depend on which node we are estimating (clevel covariate).
    if (any(h.g.ratio[index, i, ] != 0)) {

      IC[index, ] <- IC[index, ] + (Qstar.kplus1[index, i] - Qstar[index, i]) * h.g.ratio[index, i, ]

      }
  }
  return(IC)
}

################################
# FitPooledMSM
################################

#' FitPooledMSM
#'
#' Fit pooled MSM.
#'
#' @param working.msm Working marginal structural model.
#' @param Qstar Final output of the iterated conditional expectations, as created by \code{FixedTimeTMLEMediation}.
#' @param combined.summary.measures Betas for coefficient estimation. (TO DO)
#' @param msm.weights Weights for the marginal structural model.
#'
#' @return Returns coefficients obtained from the marginal structural model and predicted values for Qstar.
#'

#Some notes:
#Qstar: n x num.regimes x num.final.Ynodes
#combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
#msm.weights: n x num.regimes x num.final.Ynodes

FitPooledMSM <- function(working.msm, Qstar, combined.summary.measures, msm.weights) {

  #Number of samples
  n <- dim(Qstar)[1]
  #Number of regimes
  num.regimes <- dim(Qstar)[2]
  #Number of final Y nodes
  num.final.Ynodes <- dim(Qstar)[3]
  #Betas
  num.summary.measures <- dim(combined.summary.measures)[2]

  #For simplest MSM, just intercept.
  X <- apply(combined.summary.measures, 2, rbind)
  Y <- as.vector(Qstar)
  weight.vec <- as.vector(msm.weights)
  data.pooled <- data.frame(Y, X)
  #speedglm crashes if Y is NA even if weight is 0
  positive.weight <- weight.vec > 0

  #Estimate coefficients for betas using the working MSM as formula and samples with positive weight.
  m <- speedglm(formula(working.msm), data=data.pooled[positive.weight, ], family=quasibinomial(), weights=weight.vec[positive.weight], maxit=100)
  #Predict on all samples (even if they do not have positive weight).
  SuppressGivenWarnings(m.beta <- predict(m, newdata=data.pooled, type="response"), "prediction from a rank-deficient fit may be misleading")
  dim(m.beta) <- dim(Qstar)

  return(list(m=m, m.beta=m.beta))

}

################################
# FinalizeIC
################################

#' FinalizeIC
#'
#' Final step in calculating TMLE influence curve
#'
#' @param IC TO DO
#' @param combined.summary.measures TO DO
#' @param Qstar TO DO
#' @param m.beta TO DO
#' @param msm.weights TO DO
#' @param observation.weights
#'
#' @return Returns final Influence Curve.
#'

#Some notes:
#mBeta, Qstar: n x num.regimes x num.final.Ynodes
#combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)

#summary.measures: num.regimes x num.summary.measures x num.final.Ynodes
#msm.weights: n x num.regimes x num.final.Ynodes

FinalizeIC <- function(IC, combined.summary.measures, Qstar, m.beta, msm.weights, observation.weights) {

  num.betas <- ncol(IC)
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.final.Ynodes <- dim(Qstar)[3]

  stopifnot(num.betas == ncol(combined.summary.measures))

  finalIC <- matrix(0, nrow=n, ncol=num.betas)
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      if (any(msm.weights[, i, j] > 0)) {
        m1 <- matrix(Qstar[, i, j] - m.beta[, i, j], ncol=1)   #n x 1
        for (k in 1:num.betas) {
          m2 <- combined.summary.measures[, k, i, j] # n x 1
          finalIC[, k] <- finalIC[, k] + msm.weights[, i, j] * observation.weights * (m1 * m2)
        }
      }
    }
  }
  if (any(abs(colSums(finalIC)) > 0.001 )) {
    msg <- capture.output(cat("final IC problem", colSums(finalIC)))
    warning(paste(msg, collapse="\n"))
  }
  IC <- IC + finalIC
  return(IC)
}

################################
# FixScoreEquation
################################

#' FixScoreEquation
#'
#' In case the GLM did not converge and the TMLE does not solve the score equation, attempt to solve the score equation directly.
#'
#' @param Qstar.kplus1 Estimate of the expectation with respect to the distribution of the node one ahead of the current node given the past (dimension: n x num.regimes).
#' @param h.g.ratio I() * h * observation.weights * cum.q.ratio / g
#' @param uncensored Uncensored samples
#' @param intervention.match Samples that match the rule.
#' @param is.deterministic TO DO
#' @param deterministic.Q TO DO
#' @param off offset used for TMLE.
#' @param X TO DO
#' @param regimes.with.positive.weight Regimes with positive weight.
#'
#' @return Returns solved score equation.
#'

FixScoreEquation <- function(Qstar.kplus1, h.g.ratio, uncensored, intervention.match, is.deterministic, deterministic.Q, off, X, regimes.with.positive.weight) {

  CalcScore <- function(e) {
    Qstar <- QstarFromE(e)
    ICtemp <- CalcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
    return(sum(colSums(ICtemp) ^ 2)) #each column has to add to zero
  }

  QstarFromE <- function(e) {
    Qstar <- plogis(off + X %*% e) #X: n x (num.summary.measures + num.baseline.covariates) (which should be num.beta);  e: num.beta x 1
    dim(Qstar) <- dim(Qstar.kplus1)
    Qstar[is.deterministic] <- deterministic.Q[is.deterministic] #matrix indexing
    return(Qstar)
  }

  FindMin <- function(minimizer) {
    num.tries <- 20
    init.e <- numeric(num.betas) #first try an initial estimate of epsilon=0
    for (i in 1:num.tries) {
      m <- nlminb(start=init.e, objective=CalcScore, control=list(abs.tol=max.objective, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
      e <- m$par
      obj.val <- m$objective
      if (obj.val < max.objective) {
        m$ltmle.msg <- "updating step using glm failed to solve score equation; solved using nlminb"
        return(list(e=e, solved=TRUE, m=m))
      }
      init.e <- rnorm(num.betas) #if the first try didn't work, try a random initial estimate of epsilon
    }
    return(list(e=numeric(num.betas), solved=FALSE, m="score equation not solved!")) #nocov - return Q (not updated)
  }

  max.objective <- 0.0001 ^ 2
  num.betas <- ncol(X)
  for (offset.lbound in c(1e-8, 0.0001, 0.001, 0.01)) {
    off <- Bound(off, qlogis(c(offset.lbound, 1-offset.lbound)))
    l <- FindMin("nlminb")
    if (l$solved) break
  }
  if (! l$solved) stop("minimizer failed")
  Qstar <- QstarFromE(l$e)
  return(list(Qstar=Qstar, fit=l$m))
}

################################
# CalcGUnboundedToBoundedRatio
################################

#' CalcGUnboundedToBoundedRatio
#'
#' In case the GLM did not converge and the TMLE does not solve the score equation, attempt to solve the score equation directly.
#'
#' @param g.list TO DO
#' @param nodes All available nodes in the data.
#' @param final.Ynodes Final Y nodes.
#'
#' @return Returns...
#'

CalcGUnboundedToBoundedRatio <- function(g.list, nodes, final.Ynodes) {

  CalcForFinalYNode <- function(num.AC.nodes) {
    if (num.AC.nodes == 0) return(1)
    if (! anyNA(g.list$cum.g)) return(AsMatrix(g.list$cum.g.unbounded[, num.AC.nodes, ] / g.list$cum.g[, num.AC.nodes, ]))
    #cum.g is NA after censoring - for censored observations use cum.g.meanL
    #If censored at node j, set all nodes > j to meanl.
    #[,,k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L (na.rm=T)
    g.ratio1 <- matrix(NA, n, num.regimes)
    for (i in 1:num.regimes) {
      g.ratio.temp <- cbind(g.list$cum.g.meanL.unbounded[, num.AC.nodes, i, ] / g.list$cum.g.meanL[, num.AC.nodes, i, ], g.list$cum.g.unbounded[, num.AC.nodes, i] / g.list$cum.g[, num.AC.nodes, i])
      index <- max.col(!is.na(g.ratio.temp), "last")
      g.ratio1[, i] <- g.ratio.temp[sub2ind(1:n, col = index, num.rows = n)]
    }
    return(g.ratio1)
  }

  #calc for each final.ynode - num.AC.nodes varies
  n <- dim(g.list$cum.g)[1]
  num.regimes <- dim(g.list$cum.g)[3]
  num.final.Ynodes <- length(final.Ynodes)
  g.ratio <- array(dim=c(n, num.regimes, num.final.Ynodes))
  for (j in 1:num.final.Ynodes) {
    num.AC.nodes <- sum(nodes$AC < final.Ynodes[j])
    g.ratio[, , j] <- CalcForFinalYNode(num.AC.nodes)
  }

  return(g.ratio)
}
