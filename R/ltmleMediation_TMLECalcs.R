################################
# MainCalcsMediation
################################

#' MainCalcsMediation
#'
#' Main TMLE calculations
#'
#' @param inputs Output of \code{CreateMedInputs}
#'
#' @return Returns estimate of $E[Y_{\tau}(a, \overline{\Gamma}^{a^'})]$
#'

# Loop over final Ynodes, run main calculations...
MainCalcsMediation <- function(inputs) {

  # TO DO:
  # not doing Estimate Var in this implementation.

  #Change to Influence Curve Variance only:
  inputs$IC.variance.only <- T
  if (!exists("est.var.iptw")) est.var.iptw <<- F #fixme!

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
  Qstar <- array(dim=c(n, num.regimes, num.final.Ynodes))
  #n x num.regimes x num.final.Ynodes
  all.msm.weights <- GetMsmWeights(inputs)
  new.var.y <- array(dim=c(num.betas, num.betas, num.final.Ynodes))
  IC <- matrix(0, n, num.betas)

  #store IC for each final Ynode, compare var(IC) to sum(var(IC.ynode))
  IC.y <- array(dim=c(n, num.betas, num.final.Ynodes))

  #TO DO: Check EstimateMultiDens function
  #Estimate for all A and C nodes.
  g.abar.list <- EstimateG(inputs,regimes.use =  'regimes') ### -Wen: exisiting ltmle implmentation will use inputs$regimes to set newdata
  #Estimate p_z(z|past)
  cum.qz.abar <- EstimateMultiDens(inputs,use.regimes='regimes',use.intervention.match = 'intervention.match',is.Z.dens = T)
  #Estimate p_l(l|past)
  cum.qL.abar <- EstimateMultiDens(inputs,use.regimes='regimes',use.intervention.match = 'intervention.match',is.Z.dens = F)

  if(!setequal(inputs$regimes,inputs$regimes.prime)){

    g.abar.prime.list <- EstimateG(inputs,regimes.use =  'regimes.prime') ### -Wen: exisiting ltmle implmentation will use inputs$regimes to set newdata
    cum.qz.abar.prime <- EstimateMultiDens(inputs,use.regimes='regimes.prime',use.intervention.match = 'intervention.match.prime',is.Z.dens = T)
    cum.qL.abar.prime <- EstimateMultiDens(inputs,use.regimes='regimes.prime',use.intervention.match = 'intervention.match.prime',is.Z.dens = F)

    }else{

    g.abar.prime.list <- g.abar.list
    cum.qz.abar.prime <- cum.qz.abar
    cum.qL.abar.prime <- cum.qL.abar

    }

  #safe.solve function?
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
      fixed.tmle <- FixedTimeTMLEMediation(inputs, nodes = SubsetNodes(inputs$all.nodes, final.Ynode=inputs$final.Ynodes[j]), msm.weights = drop3(all.msm.weights[, , j, drop=FALSE]), combined.summary.measures = dropn(inputs$combined.summary.measures[, , , j, drop=FALSE], n=4), g.abar.list = g.abar.list, g.abar.prime.list=g.abar.prime.list, cum.qz.abar=cum.qz.abar, cum.qz.abar.prime=cum.qz.abar.prime, cum.qL.abar=cum.qL.abar, cum.qL.abar.prime=cum.qL.abar.prime)
      IC <- IC + fixed.tmle$IC
      IC.y[, , j] <- fixed.tmle$IC
      Qstar[, , j] <- fixed.tmle$Qstar # n x num.regimes
      new.var.y[, , j] <- fixed.tmle$est.var

    }

    if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) {
      return(list(IC=matrix(NA, 1, 1), msm=NULL, beta=qlogis(mean(Qstar)), cum.g=g.list$cum.g, cum.g.unbounded=g.list$cum.g.unbounded, fit=fixed.tmle$fit, variance.estimate=NULL, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar))
    }

    fitted.msm <- FitPooledMSM(inputs$working.msm, Qstar, inputs$combined.summary.measures, all.msm.weights * inputs$observation.weights)
    IC <- FinalizeIC(IC, inputs$combined.summary.measures, Qstar, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights) #n x num.betas
    C.old <- NormalizeIC(IC, inputs$combined.summary.measures, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights, g.ratio = NULL) #C without using g.ratio

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

    IC <- t(safe.solve(C.old, t(IC))) #IC %*% safe.solve(C)
    beta <- coef(fitted.msm$m)
    names(beta) <- inputs$beta.names

  }

  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.gq.abar=list(cum.g=g.abar.list$cum.g,cum.qL=cum.qL.abar$cum.g.block,cum.qZ=cum.qz.abar$cum.g),
              cum.gq.abar.prime=list(cum.g=g.abar.prime.list$cum.g,cum.qL=cum.qL.abar.prime$cum.g.block,cum.qZ=cum.qz.abar.prime$cum.g),
              cum.gq.abar.unbounded=list(cum.g=g.abar.list$cum.g.unbounded,cum.qL=cum.qL.abar$cum.g.unbounded.block,cum.qZ=cum.qz.abar$cum.g.unbounded),
              cum.gq.abar.prime.unbounded=list(cum.g=g.abar.prime.list$cum.g.unbounded,cum.qL=cum.qL.abar.prime$cum.g.unbounded.block,cum.qZ=cum.qz.abar.prime$cum.g.unbounded),
              fit=fixed.tmle$fit, variance.estimate=variance.estimate, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar)) #note: only returns cum.g and fit for the last final.Ynode
}

################################
# GetMsmWeights
################################

#' GetMsmWeights
#'
#' Set weights for the Marginal Structural Model
#'
#' @param inputs Output of \code{CreateMedInputs}
#'
#' @return Returns weights for the Marginal Structural Model.
#'

GetMsmWeights <- function(inputs) {

  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]

  stopifnot(num.regimes >= 1)

  num.final.Ynodes <- length(inputs$final.Ynodes)

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
#' Parametric estimation of each g-factor.
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param regimes.use
#'
#' @return Returns
#'

EstimateG <- function(inputs,regimes.use) {

  inputs$regimes <- inputs[[regimes.use]]
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
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

    for (i in nodes$A) {
      for (regime.index in 1:num.regimes) {
        regimes.meanL[is.na(regimes.meanL[, i, regime.index]), i, regime.index] <- Mode(inputs$regimes[, i, regime.index], na.rm = TRUE)
      }
    }
  } else {

    regimes.meanL <- NULL

  }

  for (i in 1:length(nodes$AC)) {

    cur.node <- nodes$AC[i]
    #Returns inputs$uncensored column that corresponds to censoring at the time point being considered. (looks at the last C node..)
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)
    #deterministic due to death or Q.function
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic

    if (is.numeric(inputs$gform)) {

      if (!inputs$IC.variance.only) stop("IC.variance.only=FALSE not currently compatible with numeric gform")
      if (!is.null(inputs$deterministic.g.function)) stop("deterministic.g.function is not compatible with numeric gform")
      prob.A.is.1[, i, ] <- inputs$gform[, i, ]
      g.est <- list(is.deterministic = deterministic.origdata) #note: this assumes that deterministic.Q.function doesn't depend on A (throw warning in CheckInputs)
      fit[[i]] <- "no fit due to numeric gform"

    } else {

      form <- inputs$gform[i]
      #deterministic due to ACnode map - using original data; now considering g
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, nodes, using.newdata=F)
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic

      if (inputs$stratify) {

        intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node=nodes$AC[i])
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata

      } else {

        #If all uncensored and no deterministic g and q, all samples.
        #Otherwise estimate from uncensored samples.
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata

      }

      #assume all regimes have positive weight for some final.Ynode
      g.est <- Estimate(inputs, form=form, Qstar.kplus1=NULL, subs=subs, family=quasibinomial(), type="response", nodes=nodes, called.from.estimate.g=TRUE, calc.meanL=!inputs$IC.variance.only, cur.node=cur.node, regimes.meanL=regimes.meanL, regimes.with.positive.weight=1:num.regimes)
      prob.A.is.1[, i, ] <- g.est$predicted.values
      fit[[i]] <- g.est$fit

    }

    #prob.A.is.1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    if (cur.node %in% nodes$A) {

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

    g[, i, ] <- CalcG(AsMatrix(prob.A.is.1[, i, ]), cur.abar, g.est$is.deterministic)

    if (!inputs$IC.variance.only) {
      for (j in sseq(1, dim(g.meanL)[4])) {
        g.meanL[, i, , j] <- CalcG(AsMatrix(g.est$prob.A.is.1.meanL[, , j]), cur.abar.meanL, g.est$is.deterministic)
      }
    }

    if (anyNA(g[uncensored, i, ])) stop("Error - NA in g. g should only be NA after censoring. If you passed numeric gform, make sure there are no NA values except after censoring. Otherwise something has gone wrong.")

  }

  for (regime.index in 1:num.regimes) {
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
#' Run GLM or SuperLearner
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
#' @return Returns
#'

Estimate <- function(inputs, form, subs, family, type, nodes, Qstar.kplus1, cur.node, calc.meanL, called.from.estimate.g, regimes.meanL, regimes.with.positive.weight) {

  library(speedglm)

  FitAndPredict <- function() {
    if (length(Y.subset) < 2) stop("Estimation failed because there are fewer than 2 observations to fit")
    if (use.glm) {
      #estimate using GLM
      # cat(" fit: ", form, "\n")
      SuppressGivenWarnings({
        m <- speedglm.wfit(Y.subset, X.subset, family=family, maxit = 100, weights=observation.weights.subset, offset=offst, intercept=intercept)
        m$terms <- tf
        class(m) <- c("speedglm", "speedlm")
        predicted.values <- predict(m, newdata=newdata, type=type)
      }, GetWarningsToSuppress())
    } else {
      #estimate using SuperLearner

      #rhs <- setdiff(RhsVars(form), rownames(alias(form, data=X.subset)$Complete))  #remove aliased (linearly dependent) columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped ## 2/2 is this really needed?
      newX.list <- GetNewX(newdata)
      SetSeedIfRegressionTesting(inputs)
      try.result <- try({
        SuppressGivenWarnings(m <- SuperLearner::SuperLearner(Y=Y.subset, X=X.subset, SL.library=SL.library, verbose=FALSE, family=family, newX=newX.list$newX, obsWeights=observation.weights.subset), c("non-integer #successes in a binomial glm!", "prediction from a rank-deficient fit may be misleading"))
      })
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
    if (use.glm) {
      predict(m, newdata1, type)
    } else {
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

  PredictProbAMeanL <- function() {
    #Predict the probability that A=1 if L and Y nodes are set to their mean (or median) values

    #probAis1.meanL is n x num.LYnodes - 1
    #probAis1.meanL[, k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L

    #somewhat inefficient - for W A.1 L.2 A.2 L.3 A.3 Y, does P(A.1=1) setting L.3 to mean and then L.2 and L.3 to mean, but none of these can be used in P(A.1=1) because they're after A.1

    #A is already set to abar in data
    probAis1.meanL <- matrix(NaN, nrow(inputs$data), length(nodes$LY) - 1)
    if (ncol(probAis1.meanL) == 0) return(probAis1.meanL)
    all.LY.nodes <- sort(union(nodes$L, nodes$Y)) #not the same as nodes$LY, which removes blocks
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
  SL.library <- if (called.from.estimate.g) inputs$SL.library.g else inputs$SL.library.Q

  #in a formula like "Y ~ 1", call glm
  use.glm <- (is.null(SL.library) || length(RhsVars(f)) == 0)

  if (use.glm) {
    #scale Lnodes to 0-1 to avoid numerical problems in speedglm
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

  if (is.null(Qstar.kplus1)) {

    data.with.Qstar <- data

  } else {

    if (is.matrix(Qstar.kplus1)) {

      data.with.Qstar <- cbind(data, Q.kplus1=Qstar.kplus1[, first.regime])

    } else {

      data.with.Qstar <- cbind(data, Q.kplus1=Qstar.kplus1)

    }
  }

  mod.frame <- model.frame(f, data = data.with.Qstar, drop.unused.levels = TRUE, na.action = na.pass)
  Y <- mod.frame[[1]]
  tf <- terms(f)
  X <- model.matrix(tf, mod.frame)
  offst <- model.offset(mod.frame)
  intercept <- attributes(tf)$intercept

  #SL does not support quasibinomial()
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

  if (calc.meanL) {
    prob.A.is.1.meanL <- array(NaN, dim=c(nrow(inputs$data), num.regimes, length(nodes$LY)-1))
    Anode.index <- which(nodes$A < cur.node)
  } else {
    prob.A.is.1.meanL <- NULL
  }

  for (regime.index in regimes.with.positive.weight) {

    newdata <- SetA(data = data.with.Qstar, regimes = inputs$regimes, Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)

    if (calc.meanL) {
      if (!is.null(regimes.meanL)) {
        newdata.meanL <- SetA(data = data.with.Qstar, regimes = regimes.meanL, Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)
      } else {
        newdata.meanL <- newdata
      }
    }

    deterministic.list.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g, inputs$survivalOutcome)

    if (called.from.estimate.g && !is.null(inputs$deterministic.g.function)) {

      newdata.with.current <- newdata
      stopifnot(cur.node %in% nodes$AC)

      if (cur.node %in% nodes$A) {
        newdata.with.current[, cur.node] <- inputs$regimes[, which(nodes$A == cur.node), regime.index] #set current node to regime for consistency checking in IsDeterministicG
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
      single.subs <- if (multiple.subs) subs[, subs.index] else subs
      X.subset <- X[single.subs, , drop=FALSE]
      X.subset[, colAlls(X.subset == 0)] <- 1 #if there is a column of all zeros, speedglm may crash - replace with column of 1s
      observation.weights.subset <- inputs$observation.weights[single.subs]
      offst.subset <- offst[single.subs]
    }

    if (regime.index == first.regime || multiple.subs || multiple.Qstar) {
      Y.subset <- Y[single.subs]
      if (anyNA(Y.subset)) stop("NA in Estimate")
    }

    if (is.numeric(form)) {

      predicted.values[, regime.index] <- form[, regime.index]  #if gform is numeric, it's a matrix of prob.A.is.1
      m <- "gform passed as numeric, no estimation took place"

    } else if (!all(deterministic.list.newdata$is.deterministic | deterministic.g.list.newdata$is.deterministic)) {

       if (is.null(fit.and.predict) || multiple.Qstar || multiple.subs) {

        fit.and.predict <- FitAndPredict()
        m <- fit.and.predict$m
        predicted.values[, regime.index] <- fit.and.predict$predicted.values

      } else {
        #just predict
        predicted.values[, regime.index] <- PredictOnly(newdata)
      }
      if (calc.meanL) prob.A.is.1.meanL[, regime.index, ] <- PredictProbAMeanL()
    } else {
      m <- "all rows are deterministic, no estimation took place"
    }

    predicted.values[deterministic.g.list.newdata$is.deterministic, regime.index] <- deterministic.g.list.newdata$prob1
    if (calc.meanL) prob.A.is.1.meanL[deterministic.g.list.newdata$is.deterministic, regime.index, ] <- deterministic.g.list.newdata$prob1
    is.deterministic[, regime.index] <- deterministic.list.newdata$is.deterministic
    if (!called.from.estimate.g) deterministic.Q[deterministic.list.newdata$is.deterministic, regime.index] <- deterministic.list.newdata$Q
    if (!use.glm && !isTRUE(attr(SL.library, "return.fit", exact = TRUE))) m <- summary(m)
    fit[[regime.index]] <- m
    if (multiple.subs) subs.index <- subs.index + 1
    if (multiple.Qstar) Qstar.index <- Qstar.index + 1
  }

  if (all(is.na(predicted.values))) stop("??") #fixme - remove this

  return(list(predicted.values=predicted.values, fit=fit, is.deterministic=is.deterministic, deterministic.Q=deterministic.Q, prob.A.is.1.meanL=prob.A.is.1.meanL))
}

################################
# EstimateMultiDens
################################

#' EstimateMultiDens
#'
#' If not is.Z.dens, estimating L densities
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param use.regimes TO DO
#' @param use.intervention.match TO DO
#' @param is.Z.dens TO DO
#'
#' @return Returns estimate of L densities.
#'

EstimateMultiDens <- function(inputs,use.regimes,use.intervention.match,is.Z.dens){
  inputs$regimes <- inputs[[use.regimes]]
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs[[use.regimes]])[3]
  nodes <- inputs$all.nodes

  if(is.Z.dens){
    dens.nodes <- nodes$Z
    dens.forms <- inputs$qzform
  }else{
    dens.nodes <- sort(c(nodes$L,nodes$Y))
    dens.forms <- inputs$qLform
  }

  fit <- vector('list',length(dens.nodes))
  g <- cum.g <- cum.g.unbounded <- prob.is.1 <- array(NaN, dim=c(n, length(dens.nodes), num.regimes))

  for (i in 1:length(dens.nodes)) {

    cur.node <- dens.nodes[i]
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic #deterministic due to death or Q.function

    if (! is.numeric(dens.forms)) {
      form <- dens.forms[i]
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, nodes, using.newdata=F) #deterministic due to acnode map - using original data
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic

      if (inputs$stratify) {
        intervention.match <- InterventionMatch(inputs[[use.intervention.match]], nodes$A, cur.node=cur.node)
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata
      } else {
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
      }
      g.est <- Estimate(inputs, form=form, Qstar.kplus1=NULL, subs=subs, family=quasibinomial(), type="response", nodes=nodes, called.from.estimate.g=TRUE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=1:num.regimes) #assume all regimes have positive weight for some final.Ynode
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
    cum.g.list <- CalcCumG(AsMatrix(g[, , regime.index]), inputs$gbounds)
    cum.g[, , regime.index] <- cum.g.list$bounded
    cum.g.unbounded[, , regime.index] <- cum.g.list$unbounded
  }
  ret.list <- list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded,prob.is.1=prob.is.1)

  ## if L node, qL for blocks of L node.
  if(!is.Z.dens){
    max.Lnodes.col <- vector('numeric',length(nodes$LY))
    for (i in 1:length(nodes$LY)) {

      cat(date(),'EstimateLYnodes node ',i,'\n')
      if(i==length(nodes$LY)){
        max.Lnodes.col[i] <- which(dens.nodes==max(dens.nodes))
      }
      else{
        max.Lnodes.col[i] <- which(dens.nodes==nodes$LY[i+1])-1 ## the desn.node that is smaller than this one
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
#' If not is.Z.dens, estimating L densities
#'
#' @param inputs Output of \code{CreateMedInputs}
#' @param cum.g.abar TO DO
#' @param cum.qz.abar TO DO
#' @param cum.qz.abar.prime TO DO
#' @param msm.weights
#'
#' @return Returns estimate of L densities.
#'

CalcIPTWMediation <- function(inputs, cum.g.abar, cum.qz.abar, cum.qz.abar.prime,msm.weights) {

  if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) {
    return(list(beta=NA, IC=matrix(NA, 1, 1)))
  }

  nodes <- inputs$all.nodes
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  num.final.Ynodes <- length(inputs$final.Ynodes)
  Y.vec <- X.mat <- weight.vec <- NULL
  save.xy <- list()

  for (j in 1:num.final.Ynodes) {

    final.Ynode <- inputs$final.Ynodes[j]
    #Returns the intervention match for the closest A node. TO DO: check how we get intervention.match, again.
    intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node=final.Ynode)
    #Check censoring for the closest C node (to Y)
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, final.Ynode)

    for (i in 1:num.regimes) {

      index <- uncensored & intervention.match[, i]
      #Look for the AC node closest to final Y
      col.index <- which.max(nodes$AC[nodes$AC < final.Ynode])
      #Look for the Z node closest to final Y
      col.z.index <- which.max(nodes$Z[nodes$Z < final.Ynode])

      #Limit to samples matching index (uncensored and intervention match)
      Y <- inputs$data[index, final.Ynode]

      #Return cum.g.bar for closest node to Y for regime i (for sample ok by index).
      #This is p_A(a|past)
      g <- cum.g.abar[index, col.index, i]

      #Ratio part of IPTW and clever covariate.
      #This is p_z(z|abar.prime,past)/p_z(z|abar,past)
      qz.ratio <- cum.qz.abar.prime[index, col.z.index, i]/cum.qz.abar[index, col.z.index, i]

      #Dimensions are n x ? x num.regimes x num.final.Ynodes
      X <- inputs$combined.summary.measures[index, , i, j]

      if (is.vector(X)) { #if only one summary.measure or sum(index)==1, X is dropped to vector
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

  #this happens if there are no rows uncensored and intervention.match
  if (nrow(X.mat) == 0) {

    warning("no rows uncensored and matching regimes/abar - IPTW returns NA")
    num.beta <- ncol(inputs$combined.summary.measures)
    return(list(beta=rep(NA, num.beta), IC=matrix(nrow=n, ncol=num.beta)))

  }

  #working.msm: why Y ~ -1 + S1? TO DO (no intercept?)
  #working.msm: note that S1 will be X.mat, Y is Y.vec and finally we have weight.vec. Will estimate coefficient for S1 (beta)
  #note: scale weights because there were rare problems where large weights caused convergence problems
  m.glm <- speedglm(formula(inputs$working.msm), family=quasibinomial(), data=data.frame(Y=Y.vec, X.mat, weight.vec), weights=as.vector(scale(weight.vec, center=FALSE)))
  beta <- coef(m.glm)

  #n x num.betas
  IC <- matrix(0, nrow=n, ncol=length(beta))
  #n x num.regimes x num.final.Ynodes
  m.beta <- array(dim=c(n, num.regimes, num.final.Ynodes))
  cnt <- 0

  for (j in 1:num.final.Ynodes) {

    final.Ynode <- inputs$final.Ynodes[j]

    for (i in 1:num.regimes) {

      newdata <- data.frame(inputs$combined.summary.measures[, , i, j])
      colnames(newdata) <- colnames(inputs$combined.summary.measures) #needed if only one summary measure
      SuppressGivenWarnings(m.beta[, i, j] <- predict(m.glm, newdata=newdata, type="response"), "prediction from a rank-deficient fit may be misleading")

      cnt <- cnt + 1
      XY.list <- save.xy[[cnt]]

      #recycles weight, Y, m.beta
      #IC for each uncensored and intervention.match sample
      #weight is the clever covariate; penalize the difference between the Y and its estimate, m.beta
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
#' @param nodes TO DO
#' @param msm.weights TO DO
#' @param combined.summary.measures TO DO
#' @param g.abar.list TO DO
#' @param g.abar.prime.list TO DO
#' @param cum.qz.abar TO DO
#' @param cum.qz.abar.prime TO DO
#' @param cum.qL.abar TO DO
#' @param cum.qL.abar.prime TO DO
#'
#' @return
#'

FixedTimeTMLEMediation <- function(inputs, nodes, msm.weights, combined.summary.measures, g.abar.list , g.abar.prime.list, cum.qz.abar, cum.qz.abar.prime, cum.qL.abar, cum.qL.abar.prime){

  #combined.summary.measures: n x num.measures x num.regimes
  #(num.measures=num.summary.measures + num.baseline.covariates)
  LYZnodes <- sort(c(nodes$LY, nodes$Z))
  data <- inputs$data

  num.regimes <- dim(inputs$regimes)[3]
  n <- nrow(data)
  num.betas <- ncol(combined.summary.measures)
  tmle <- rep(NA, num.regimes)
  IC <- matrix(0, nrow=n, ncol=num.betas)

  est.var <- matrix(0, num.betas, num.betas)
  regimes.with.positive.weight <- which(apply(msm.weights > 0, 2, any))

  if (length(regimes.with.positive.weight) == 0) stop("All regimes have weight 0 (one possible reason is that msm.weights='emipirical' and no data rows match any of the regimes and are uncensored)")

  fit.Qstar <- vector("list", length(LYZnodes))
  names(fit.Qstar) <- names(data)[LYZnodes]
  fit.Q <-  vector("list", length(regimes.with.positive.weight))

  for (i in regimes.with.positive.weight) fit.Q[[i]] <- fit.Qstar

  Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow=n, ncol=num.regimes)
  mean.summary.measures <- apply(abs(combined.summary.measures), 2, mean)

  #at and after cur.node
  for (i in length(LYZnodes):1) {

    cur.node <- LYZnodes[i]
    #TRUE if death in previous Y or censored (for current node).
    deterministic.list.origdata <- IsDeterministic(data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
    #Looks at all C nodes, with the focus on the closest one to the cur.node
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)

    ## if at QLupdate:
    #UPDAYE Q_L!
    if(cur.node %in% nodes$LY){

      #Only relevant for stratify
      intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node)

      if (inputs$stratify) {

        subs <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic #n x num.regimes

        } else {

        #Uncensored and non-deterministic samples
        subs <- uncensored & !deterministic.list.origdata$is.deterministic #vector

        }

      #if this is the last node, only pass the first column as a vector
      Q.est <- Estimate(inputs, form = inputs$QLform[which(nodes$LY==cur.node)], Qstar.kplus1=if (i == length(LYZnodes)) Qstar.kplus1[, 1] else Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight)
      logitQ <- Q.est$predicted.values
      fit.Q[[i]] <- Q.est$fit
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      Znode.index <- which.max(nodes$Z[nodes$Z < cur.node])
      update.list <- UpdateQMediation(Qstar.kplus1, logitQ, combined.summary.measures, cum.g = if(length(ACnode.index)==0) 1 else g.abar.list$cum.g[, ACnode.index, ],cum.q.ratio=if(length(Znode.index)==0) 1 else cum.qz.abar.prime$cum.g[,Znode.index,]/cum.qz.abar$cum.g[,Znode.index,], inputs$working.msm, uncensored, intervention.match, deterministic.list.origdata$is.deterministic, msm.weights, inputs$gcomp, inputs$observation.weights)
      Qstar <- update.list$Qstar
      Qstar[Q.est$is.deterministic] <- Q.est$deterministic.Q[Q.est$is.deterministic] #matrix indexing
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) / mean.summary.measures

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
    }## done update QL

    #If current node is Z, update Q_Z!
    if(cur.node %in% nodes$Z){

      intervention.match <- InterventionMatch(inputs$intervention.match.prime, nodes$A, cur.node)

      if (inputs$stratify) {

        subs <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic #n x num.regimes

      } else {

        subs <- uncensored & !deterministic.list.origdata$is.deterministic #vector

      }

      #Issue
      #Q.est <- Estimate(inputs = set(inputs,'regimes',inputs$regimes.prime), form = inputs$QZform[which(nodes$Z==cur.node)], Qstar.kplus1=Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight)
      Q.est <- Estimate(inputs, form = inputs$QZform[which(nodes$Z==cur.node)], Qstar.kplus1=Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight)
      logitQ <- Q.est$predicted.values
      fit.Q[[i]] <- Q.est$fit
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      Lnode.index <- which.max(nodes$LY[nodes$LY < cur.node])
      update.list <- UpdateQMediation(Qstar.kplus1, logitQ, combined.summary.measures, cum.g = if(length(ACnode.index)==0) 1 else g.abar.prime.list$cum.g[, ACnode.index, ],cum.q.ratio=if(length(Lnode.index)==0) 1 else cum.qL.abar$cum.g.block[,Lnode.index,]/cum.qL.abar.prime$cum.g.block[,Lnode.index,], inputs$working.msm, uncensored, intervention.match, deterministic.list.origdata$is.deterministic, msm.weights, inputs$gcomp, inputs$observation.weights)
      Qstar <- update.list$Qstar
      Qstar[Q.est$is.deterministic] <- Q.est$deterministic.Q[Q.est$is.deterministic] #matrix indexing
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) / mean.summary.measures

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
    }

    #
    # est.var <- est.var + EstimateVariance(inputs, nodes, combined.summary.measures, regimes.with.positive.weight, uncensored, alive=!deterministic.list.origdata$is.deterministic, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, g.list$cum.g, g.list$prob.A.is.1, g.list$cum.g.meanL, g.list$cum.g.unbounded, g.list$cum.g.meanL.unbounded, inputs$observation.weights, is.last.LYnode=(LYnode.index==length(nodes$LY)), intervention.match) #fixme - remove intervention.match if not using est.var.iptw

    IC <- IC + curIC

  }
  #tmle <- colMeans(Qstar)
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
#' Calculate bounded cumulative G
#'
#' @param g TO DO
#' @param gbounds TO DO
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
#' Normalize the influence curve matrix
#'
#' @param IC TO DO
#' @param combined.summary.measures TO DO
#' @param m.beta TO DO
#' @param msm.weights TO DO
#' @param observation.weights TO DO
#' @param g.ratio TO DO
#'
#' @return Returns
#'

NormalizeIC <- function(IC, combined.summary.measures, m.beta, msm.weights, observation.weights, g.ratio) {
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  #g.ratio = g.unbounded / g.bounded : n x num.regimes x num.final.Ynodes
  n <- nrow(IC)
  num.betas <- ncol(IC)
  num.regimes <- dim(combined.summary.measures)[3]
  num.final.Ynodes <- dim(combined.summary.measures)[4]

  if (is.null(g.ratio)) {
    g.ratio <- array(1, dim=c(n, num.regimes, num.final.Ynodes)) #if IC.variance.only, g.ratio should be NULL
  }
  C <- array(0, dim=c(num.betas, num.betas))
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
#' @param Qstar.kplus1 TO DO
#' @param logitQ TO DO
#' @param combined.summary.measures TO DO
#' @param cum.g TO DO
#' @param cum.q.ratio TO DO
#' @param working.msm TO DO
#' @param uncensored TO DO
#' @param intervention.match TO DO
#' @param is.deterministic TO DO
#' @param msm.weights TO DO
#' @param gcomp TO DO
#' @param observation.weights TO DO
#'
#' @return Returns targeted Q.
#'

#Some notes:
#logitQ, Qstar.kplus1: n x num.regimes
#cum.g: n x num.regimes (already indexed for this node)
#uncensored: n x 1
#intervention.match: n x num.regimes
#is.deterministic: n x 1
#summary.measures: num.regimes x num.summary.measures
#baseline.covariates: names/indicies: num.baseline.covariates x 1
#msm.weights: n x num.regimes
#stacked.summary.measures: (n*num.regimes) x num.measures
#combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)
#h.g.ratio: n x num.regimes x num.measures
#observation.weights: n x 1

UpdateQMediation <- function(Qstar.kplus1, logitQ, combined.summary.measures, cum.g, cum.q.ratio, working.msm, uncensored, intervention.match, is.deterministic, msm.weights, gcomp, observation.weights) {

  n <- nrow(logitQ)
  num.regimes <- ncol(logitQ)
  off <- as.vector(logitQ)
  Y <- as.vector(Qstar.kplus1)

  stacked.summary.measures <- apply(combined.summary.measures, 2, rbind)

  subs.vec <- uncensored & !is.deterministic & as.vector(intervention.match) #recycles uncensored and is.deterministic
  weight.vec <- numeric(n * num.regimes)
  weight.vec[subs.vec] <- (observation.weights * as.vector(msm.weights)*as.vector(cum.q.ratio)/ as.vector(cum.g))[subs.vec] #recycles observation.weights (subsetting avoids problems with NA in cum.g)
  if (anyNA(weight.vec)) stop("NA in weight.vec")

  f <- as.formula(paste(working.msm, "+ offset(off)"))
  data.temp <- data.frame(Y, stacked.summary.measures, off)
  if (gcomp) {
    Qstar <- plogis(logitQ)
    m <- "no Qstar fit because gcomp=TRUE (so no updating step)"
  } else {
    if (any(weight.vec > 0)) {
      SuppressGivenWarnings(m <- speedglm(f, data=data.temp[weight.vec > 0, ], family=quasibinomial(), weights=as.vector(scale(weight.vec[weight.vec > 0], center=FALSE)), maxit=100), GetWarningsToSuppress(TRUE)) #this should include the indicators
      SuppressGivenWarnings(Qstar <- matrix(predict(m, newdata=data.temp, type="response"), nrow=nrow(logitQ)), GetWarningsToSuppress(TRUE))  #this should NOT include the indicators  #note: could also use plogis(off + X %*% coef(m)) [but this has problems with NAs in coef(m)?]
    } else {
      Qstar <- plogis(logitQ)
      m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
    }

  }
  indicator <- matrix(uncensored * observation.weights, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * matrix(intervention.match, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) #I(A=rule and uncensored) * observation.weights
  h.g.ratio <- stacked.summary.measures * matrix(cum.q.ratio/cum.g, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * indicator # I() * h * observation.weights * cum.q / g
  dim(h.g.ratio) <- c(n, num.regimes, ncol(h.g.ratio))
  for (i in 1:num.regimes) {
    h.g.ratio[, i, ] <- h.g.ratio[, i, ] * msm.weights[, i] #recycles msm.weights
    weight.zero.index <- msm.weights[, i] == 0
    h.g.ratio[weight.zero.index, i, ] <- 0  #cum.g is 0 so X is NA so h.g.ratio is NA when weight is 0
  }
  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=stacked.summary.measures, off=off, fit=m))
}

################################
# CalcIC
################################

#' CalcIC
#'
#' Calculate the TMLE influence curve for one node.
#'
#' @param Qstar.kplus1 TO DO
#' @param Qstar TO DO
#' @param h.g.ratio TO DO
#' @param uncensored TO DO
#' @param intervention.match TO DO
#' @param regimes.with.positive.weight TO DO
#'
#' @return Returns IC for one node.
#'

CalcIC <- function(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight) {
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.betas <- dim(h.g.ratio)[3] #h.g.ratio: n x num.regimes x num.betas

  IC <- matrix(0, nrow=n, ncol=num.betas)
  for (i in regimes.with.positive.weight) {
    index <- uncensored & intervention.match[, i]
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
#' @param working.msm TO DO
#' @param Qstar TO DO
#' @param combined.summary.measures TO DO
#' @param msm.weights TO DO
#'
#' @return Returns ...
#'

#Some notes:
#Qstar: n x num.regimes x num.final.Ynodes
#combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
#msm.weights: n x num.regimes x num.final.Ynodes

FitPooledMSM <- function(working.msm, Qstar, combined.summary.measures, msm.weights) {

  n <- dim(Qstar)[1]
  num.regimes <- dim(Qstar)[2]
  num.final.Ynodes <- dim(Qstar)[3]
  num.summary.measures <- dim(combined.summary.measures)[2]

  X <- apply(combined.summary.measures, 2, rbind)
  Y <- as.vector(Qstar)
  weight.vec <- as.vector(msm.weights)
  data.pooled <- data.frame(Y, X)
  positive.weight <- weight.vec > 0 #speedglm crashes if Y is NA even if weight is 0

  m <- speedglm(formula(working.msm), data=data.pooled[positive.weight, ], family=quasibinomial(), weights=weight.vec[positive.weight], maxit=100)
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

