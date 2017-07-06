################################
# HELPER Functions
################################

################################
# IsUncensored()
################################

#' IsUncensored
#'
#' Determine which patients are uncensored
#'
#' @param uncensored.matrix TO DO
#' @param Cnodes TO DO
#' @param cur.node TO DO
#'
#' @return Returns vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node
#'

# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
IsUncensored <- function(uncensored.matrix, Cnodes, cur.node) {
  index <- which.max(Cnodes[Cnodes < cur.node])
  if (length(index) == 0) return(rep(TRUE, nrow(uncensored.matrix)))
  return(uncensored.matrix[, index])
}

################################
# IsDeterministic()
################################

#' IsDeterministic
#'
#' Determine which patients have died or have Q set deterministically by user function before cur.node.
#'
#' @param data TO DO
#' @param cur.node TO DO
#' @param deterministic.Q.function TO DO
#' @param nodes TO DO
#' @param called.from.estimate.g TO DO
#' @param survivalOutcome TO DO
#'
#' @return Return list contains:
#' is.deterministic: vector of [numObservations x 1] - true if patient is already dead before cur.node or set by deterministic.Q.function
#' Q.value: vector of [which(is.deterministic) x 1] - value of Q

IsDeterministic <- function(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g, survivalOutcome) {
  #set Q.value to 1 if previous y node is 1
  if (survivalOutcome && any(nodes$Y < cur.node)) {
    last.Ynode <- max(nodes$Y[nodes$Y < cur.node])
    is.deterministic <- data[, last.Ynode] %in% TRUE
  } else {
    is.deterministic <- rep(FALSE, nrow(data))
  }

  #get Q values from deterministic.Q.function
  default <- list(is.deterministic=is.deterministic, Q.value=1)
  if (is.null(deterministic.Q.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.Q.function(data=data, current.node=cur.node, nodes=nodes, called.from.estimate.g=called.from.estimate.g)
  if (is.null(det.list)) return(default)
  if (called.from.estimate.g) {
    #it's ok if Q.value isn't returned if called.from.estimate.g
    if (!is.list(det.list) || !("is.deterministic" %in% names(det.list)) || !(length(det.list) %in% 1:2)) stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
  } else {
    if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "Q.value")) || length(det.list) != 2) stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
  }

  if (! length(det.list$Q.value) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'Q.value' element of deterministic.Q.function's return argument should be either 1 or length(which(det.list$is.deterministic))")

  #check that these observations where Q.value is 1 due to death (previous y is 1) aren't set to anything conflicting by deterministic.Q.function
  Q.value.from.function <- rep(NA, nrow(data))
  Q.value.from.function[det.list$is.deterministic] <- det.list$Q.value
  set.by.function.and.death <- is.deterministic & det.list$is.deterministic
  if (any(Q.value.from.function[set.by.function.and.death] != 1)) {
    stop(paste("inconsistent deterministic Q at node:", names(data)[cur.node]))
  }
  finalY <- data[, max(nodes$Y)]
  inconsistent.rows <- (det.list$Q.value %in% c(0,1)) & (det.list$Q.value != finalY[det.list$is.deterministic]) & !is.na(finalY[det.list$is.deterministic])
  if (any(inconsistent.rows)) stop(paste("At node:",names(data)[cur.node], "deterministic.Q.function is inconsistent with data - Q.value is either 0 or 1 but this does not match the final Y node value\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" ")))

  #return combined values
  Q.value <- rep(NA, nrow(data))
  Q.value[is.deterministic] <- 1
  Q.value[det.list$is.deterministic] <- det.list$Q.value
  is.deterministic <- is.deterministic | det.list$is.deterministic
  Q.value <- Q.value[is.deterministic]
  if (anyNA(is.deterministic) || anyNA(Q.value)) stop("NA in is.deterministic or Q.value")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))
}

################################
# CalcUncensoredMatrix()
################################

#' CalcUncensoredMatrix
#'
#' Determine which patients have died or have Q set deterministically by user function before cur.node.
#'
#' @param data TO DO
#' @param Cnodes TO DO
#'
#' @return Calc logical matrix - n x numCnodes = is uncensored up to and including Cnode[i]

CalcUncensoredMatrix <- function(data, Cnodes) {
  uncensored <- matrix(nrow=nrow(data), ncol=length(Cnodes))
  cum.uncensored <- rep(TRUE, nrow(data))
  for (Cnode.index in seq_along(Cnodes)) {
    cum.uncensored <- cum.uncensored & (data[, Cnodes[Cnode.index]] %in% c("uncensored", NA))
    uncensored[, Cnode.index] <- cum.uncensored
  }
  return(uncensored)
}

################################
# RhsVars()
################################

#' RhsVars
#'
#' Returns the right hand side variables of a formula.
#'
#' @param f TO DO
#'
#' @return Return the right hand side variables of formula f as a character vector

RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

################################
# NodeToIndex()
################################

#' NodeToIndex
#'
#' Convert named nodes to indicies of nodes
#'
#' @param data TO DO
#' @param node TO DO
#'
#' @return Returns index of the nodes.
#'

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

################################
# CreateLYNodes()
################################

#' CreateLYNodes
#'
#' Get the LY nodes without "blocks" of L/Y nodes uninterrupted by A/C/Z nodes
#'
#' @param data TO DO
#' @param nodes TO DO
#' @param check.Qform TO DO
#' @param Qform TO DO
#'
#' @return Returns LYnodes.
#'

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

################################
# CreateNodes()
################################

#' CreateNodes
#'
#' Get all the nodes from input data.
#'
#' @param data TO DO
#' @param Anodes TO DO
#' @param Cnodes TO DO
#' @param Lnodes TO DO
#' @param Ynodes TO DO
#' @param Znodes TO DO
#'
#' @return Returns all nodes.
#'

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

################################
# GetLibrary()
################################

#' GetLibrary
#'
#' Get Super Learner library by type (Q or g).
#'
#' @param SL.library TO DO
#' @param estimate.type TO DO
#'
#' @return Returns Super Learner library.
#'

GetLibrary <- function(SL.library, estimate.type) {
  if (is.null(names(SL.library))) return(SL.library)
  if (! identical(sort(names(SL.library)), sort(c("Q", "g")))) stop("If SL.library has names, it must have two names: Q and g")
  if (! estimate.type %in% c("Q", "g")) stop("bad estimate.type")
  if (length(setdiff(names(attributes(SL.library)), c("names", "return.fit"))) > 0) stop("If SL.library has attributes, the only valid attributes are name and return.fit")
  lib <- SL.library[[estimate.type]]
  attr(lib, "return.fit") <- attr(SL.library, "return.fit", exact = TRUE)
  return(SL.library[[estimate.type]])
}

################################
# CleanData()
################################

#' CleanData
#'
#' Set all nodes (except Y) to NA after death or censoring; Set Y nodes to 1 after death.
#'
#' @param data TO DO
#' @param nodes TO DO
#' @param deterministic.Q.function TO DO
#' @param survivalOutcome TO DO
#' @param showMessage TO DO
#'
#' @return Returns cleaned data.
#'

data <- CleanData(data, all.nodes, deterministic.Q.function, survivalOutcome)

CleanData <- function(data, nodes, deterministic.Q.function, survivalOutcome, showMessage=TRUE) {
  #make sure binaries have already been converted before calling this function
  options(warn=-1)

  is.nan.df <- function (x) {
    y <- if (length(x)) {
      do.call("cbind", lapply(x, "is.nan"))
    } else {
      matrix(FALSE, length(row.names(x)), 0)
    }
  }

  is.na.strict <- function (x) is.na(x) & !is.nan.df(x)  #only for data.frames
  changed <- FALSE

  #uncensored and alive
  ua <- rep(TRUE, nrow(data))
  if (ncol(data) == 1) return(data)

  deterministic.Q.function.depends.on.called.from.estimate.g <- length(grep("called.from.estimate.g", as.character(body(deterministic.Q.function)))) > 0

  #ISSUE sometimes.
  for (i in 1:(ncol(data)-1)) {

    if (anyNA(data[ua, 1:i])) stop("NA values are not permitted in data except after censoring or a survival event")
    is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome)$is.deterministic #check determinisitic including node i
    #is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome)

    if (deterministic.Q.function.depends.on.called.from.estimate.g) {
      is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)$is.deterministic
      #is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)
      if (any(is.deterministic[ua] & !is.deterministic.Q[ua])) stop("Any row set deterministic by deterministic.Q.function(..., called.from.estimate.g=TRUE) must imply that the row is also set deterministic by deterministic.Q.function(..., called.from.estimate.g=FALSE)") #det.Q.fun(T) should imply det.Q.fun(F)
    }

    ua[ua] <- !is.deterministic[ua]

    if (anyNA(ua)) stop("internal medltmle error - ua should not be NA in CleanData")

    if ( !all(is.na.strict( data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y), drop=FALSE])) ) {
      #if deterministic, set all nodes except Y to NA
      data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y)] <- NA
      changed <- TRUE
    }

    if (i %in% nodes$C) {
      #ua: uncensored and alive
      censored <- data[, i] == "censored" & ua

      #If not already set to NA, set to NA for all censored samples (including Y).
      if (! all(is.na.strict(data[censored, (i+1):ncol(data), drop=FALSE]))) {
        data[censored, (i+1):ncol(data)] <- NA
        changed <- TRUE
      }

      ua[ua] <- !censored[ua]

      if (anyNA(ua)) stop("internal medltmle error - ua should not be NA in CleanData")
    }
  }

  if (changed && showMessage) {
    message("Note: for internal purposes, all nodes after a censoring event are set to NA and \n all nodes (except Ynodes) are set to NA after Y=1 if survivalFunction is TRUE (or if specified by deterministic.Q.function).\n Your data did not conform and has been adjusted. This may be relevant if you are \n writing your own deterministic function(s) or debugging medltmle.")
  }
  return(data)
}

################################
# TransformOutcomes()
################################

#' TransformOutcomes
#'
#' Get the range of Y. Transform outcome to be in the 0-1 range, if necessary.
#'
#' @param data TO DO
#' @param nodes TO DO
#' @param Yrange TO DO
#'
#' @return Returns data with outcome in the 0-1 range.
#'
TransformOutcomes(data, all.nodes, Yrange)

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
    #if Yrange was not specified, get it:
    Yrange <- range(all.Y, na.rm=TRUE)
    if (min(Yrange) < 0 || max(Yrange) > 1) {
      #Will need to transform outcome if outside the 0-1 range.
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

################################
# ConvertToMainTerms()
################################

#' ConvertToMainTerms
#'
#' Converts a general formula to a main terms formula and combine summary measures with baseline covariates.
#' Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2.
#'
#' @param data TO DO
#' @param msm TO DO
#' @param summary.measures TO DO
#' @param nodes
#'
#' @return Returns main terms Marginal Structural Model, summary measures, beta names and baseline column names.
#'

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

################################
# CalcInterventionMatchArray()
################################

#' CalcInterventionMatchArray
#'
#' Calculates a logical array - n x num.Anodes x num.regimes = follows regime[j] up to and including Anode[i]
#'
#' @param data TO DO
#' @param regimes TO DO
#' @param Anodes Index of A nodes. Can be easily derived from \code{all.nodes$A}.
#'
#' @return Returns intervention that matches the array.
#'

CalcInterventionMatchArray <- function(data, regimes, Anodes) {
  num.regimes <- dim(regimes)[3]
  intervention.match <- array(dim=c(nrow(data), length(Anodes), num.regimes))
  cum.intervention.match <- matrix(TRUE, nrow(data), num.regimes)

  for (Anode.index in seq_along(Anodes)) {
    #Check if A as in data equals the regime. NA will count as a match as well.
    cum.intervention.match <- cum.intervention.match & ((data[, Anodes[Anode.index]] == regimes[, Anode.index, ]) %in% c(TRUE, NA)) #recycles regimes
    intervention.match[, Anode.index, ] <- cum.intervention.match
  }
  return(intervention.match)
}

################################
# ConvertCensoringNodes()
################################

#' ConvertCensoringNodes
#'
#' Convert censoring nodes stored as binaries into factors (factors are recommended but binaries are currently accepted)
#'
#' @param data TO DO
#' @param Cnodes TO DO
#' @param as.deterministic.functions TO DO
#'
#' @return Returns data with converted censoring nodes.
#'

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

################################
# RegimesFromAbar()
################################

#' RegimesFromAbar
#'
#' Creates regimes from abar and abar.prime
#'
#' @param data TO DO
#' @param abar TO DO
#' @param rule TO DO
#'
#' @return Returns regimes.
#'

RegimesFromAbar <- function(data, abar, rule) {
  if (!is.null(rule)) {
    if (!(missing(abar) || is.null(abar))) stop("'abar' should not be specified when using a 'rule' function")
    abar <- t(apply(data, 1, rule))
  }
  if (is.vector(abar)) {
    abar <- matrix(rep(abar, each=nrow(data)), nrow=nrow(data))
  } else if (is.null(abar)) {
    abar <- matrix(nrow=nrow(data), ncol=0)
  }
  regimes <- abar
  dim(regimes) <- c(nrow(regimes), ncol(regimes), 1)
  return(regimes)
}

################################
# IsDeterministicG()
################################

#' IsDeterministicG
#'
#' Determines which patients have an Anode value which is deterministic. For example, deterministic.g.function may be used to
#' specify that once a patient starts treatment, they stay on treatment and this should be taken into
#' consideration during estimation of G.
#'
#' @param data TO DO
#' @param cur.node TO DO
#' @param deterministic.g.function TO DO
#' @param nodes TO DO
#' @param using.newdata TO DO
#'
#' @return Returns TRUE for patients that have a deterministic Anode value.
#'

IsDeterministicG <- function(data, cur.node, deterministic.g.function, nodes, using.newdata) {
  default <- list(is.deterministic=rep(FALSE, nrow(data)), prob1=NULL)
  if (is.null(deterministic.g.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.g.function(data=data, current.node=cur.node, nodes=nodes)
  if (is.null(det.list)) return(default)
  if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "prob1")) || length(det.list) != 2) stop("deterministic.g.function should return a list with names: is.deterministic, prob1")
  if (! length(det.list$prob1) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'prob1' element of deterministic.g.function's return argument should be either 1 or length(which(det.list$is.deterministic))")

  inconsistent.rows <- (det.list$prob1 %in% c(0,1)) & (det.list$prob1 != data[det.list$is.deterministic, cur.node]) & !is.na(data[det.list$is.deterministic, cur.node])
  if (any(inconsistent.rows)) {
    err.msg <- paste("At node:",names(data)[cur.node], "deterministic.g.function is inconsistent with data - prob1 is either 0 or 1 but this does not match the node value.\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" "))
    if (using.newdata) {
      err.msg <- paste(err.msg, "\n This error occured while calling deterministic.g.function on data where Anodes are set to abar.")
      cat("deterministic.g.function is inconsistent with data.\nAfter setting Anodes to abar, the data looks like this:\n")
      print(head(data[det.list$is.deterministic[inconsistent.rows], ]))
    }
    stop(err.msg)
  }
  return(det.list)
}

################################
# InterventionMatch()
################################

#' InterventionMatch
#'
#' Determine which patients are following specified treatment regime.
#'
#' @param intervention.match.array TO DO
#' @param Anodes TO DO
#' @param cur.node TO DO
#'
#' @return Returns matrix of [numObservations x numRegimes] I(A==abar) from Anodes[1] to the Anode just before cur.node.
#'

# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
InterventionMatch <- function(intervention.match.array, Anodes, cur.node) {
  index <- which.max(Anodes[Anodes < cur.node])
  if (length(index) == 0) return(matrix(TRUE, nrow(intervention.match.array), dim(intervention.match.array)[3]))
  return(as.matrix(intervention.match.array[, index, ]))
}

################################
# ConvertCensoringNodesToBinary()
################################

#' ConvertCensoringNodesToBinary
#'
#' Before passing data to SuperLearner, convert factors to binary.
#'
#' @param data TO DO
#' @param Cnodes TO DO
#'
#' @return Returns data with censoring nodes converted to binary.
#'

ConvertCensoringNodesToBinary <- function(data, Cnodes) {
  CensoringToBinary <- function(x) {
    if (! all(levels(x) %in% c("censored", "uncensored"))) {
      stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)")
    }
    b <- rep(NA_integer_, length(x))
    b[x == "censored"] <- 0L
    b[x == "uncensored"] <- 1L
    return(b)
  }

  error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored' \n (binary is also accepted, where 0=censored, 1=uncensored, but is not recommended)"
  for (i in Cnodes) {
    col <- data[, i]
    if (is.numeric(col)) {
      if (! all(col %in% c(0, 1, NA))) stop(error.msg)
    } else if (is.factor(col)) {
      data[, i] <- CensoringToBinary(col)
    } else {
      stop(error.msg)
    }
  }
  return(data)
}

################################
# SetA()
################################

#' SetA
#'
#' Set the Anodes of data to regime[, , regime.index] up to cur.node
#'
#' @param data TO DO
#' @param regimes TO DO
#' @param Anodes TO DO
#' @param regime.index TO DO
#' @param cur.node TO DO
#'
#' @return Returns data with Anodes set to regime.
#'

SetA <- function(data, regimes, Anodes, regime.index, cur.node) {
  Anode.index <- which(Anodes < cur.node)
  data[, Anodes[Anode.index]] <- regimes[, Anode.index, regime.index]
  return(data)
}

################################
# SuppressGivenWarnings
################################

SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

################################
# safe.solve
################################

#Strange errors were reported on solaris-sparc, this attempts to avoid them
safe.solve <- function(a, b) {
  if (missing(b)) {
    try.result <- try(x <- solve(a))
  } else {
    try.result <- try(x <- solve(a, b))
  }
  if (inherits(try.result, "try-error")) {
    if (missing(b)) {
      x <- matrix(nrow = nrow(a), ncol = ncol(a))
    } else {
      x <- matrix(nrow = ncol(a), ncol = ncol(AsMatrix(b)))
    }
    warning("Error in solve(), standard errors not available")
  }
  return(x)
}

################################
# AsMatrix
################################

# If x is a matrix, keep it; if x is a vector, make it a 1 column matrix
AsMatrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  } else if (is.vector(x)) {
    dim(x) <- c(length(x), 1)
    return(x)
  } else {
    stop("AsMatrix input should be a matrix or vector") # nocov (should never occur - ignore in code coverage checks)
  }
}

################################
# drop3
################################

#if x is an array with 3 dimensions and third dimension has one level, return a matrix with it dropped; otherwise error
drop3 <- function(x) {

  return(dropn(x, 3))
}

################################
# dropn
################################

#if x is an array with n dimensions and nth dimension has one level, return a matrix with it dropped; otherwise error
dropn <- function(x, n) {

  stopifnot(length(dim(x))==n)
  stopifnot(dim(x)[n]==1)
  dn <- dimnames(x)
  dim(x) <- dim(x)[1:(n-1)]
  dimnames(x) <- dn[1:(n-1)]
  return(x)

}

################################
# Bound
################################

Bound <- function(x, bounds) {
  stopifnot(length(bounds) == 2 && !anyNA(bounds))
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}
