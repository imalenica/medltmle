#' CheckMediationInputs
#'
#' Error checking for inputs.
#'

##########################################
# TO DO: Will need to fix this function.
##########################################

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}

# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

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

# Determine which patients are uncensored
#return vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
IsUncensored <- function(uncensored.matrix, Cnodes, cur.node) {
  index <- which.max(Cnodes[Cnodes < cur.node])
  if (length(index) == 0) return(rep(TRUE, nrow(uncensored.matrix)))
  return(uncensored.matrix[, index])
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

#Error checking for inputs
CheckMediationInputs <- function(data, nodes, survivalOutcome, QLform, QZform, qLform, qzform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes,regimes.prime, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp, IC.variance.only) {

  stopifnot(length(dim(regimes)) == 3)
  num.regimes <- dim(regimes)[3]
  stopifnot(length(dim(regimes.prime)) == 3)
  num.regimes.prime <- dim(regimes.prime)[3]
  if (!all(is.null(GetLibrary(SL.library, "Q")), is.null(GetLibrary(SL.library, "g")))) {
    if (!requireNamespace("SuperLearner")) stop("SuperLearner package is required if SL.library is not NULL")
  }
  #each set of nodes should be sorted - otherwise causes confusion with gform, Qform, abar
  if (is.unsorted(nodes$A, strictly=TRUE)) stop("Anodes must be in increasing order")
  if (is.unsorted(nodes$C, strictly=TRUE)) stop("Cnodes must be in increasing order")
  if (is.unsorted(nodes$L, strictly=TRUE)) stop("Lnodes must be in increasing order")
  if (is.unsorted(nodes$Y, strictly=TRUE)) stop("Ynodes must be in increasing order")
  if (is.unsorted(nodes$Z, strictly=TRUE)) stop("Znodes must be in increasing order")
  if (is.unsorted(final.Ynodes, strictly=TRUE)) stop("final.Ynodes must be in increasing order")

  if (length(nodes$L) > 0) {
    if (min(nodes$L) < min(nodes$ACZ)) stop("Lnodes are not allowed before A/C nodes. If you want to include baseline nodes, include them in data but not in Lnodes")
    if (max(nodes$L) > max(nodes$Y)) stop("Lnodes are not allowed after the final Y node")
  }
  if (min(nodes$Y) < min(nodes$ACZ)) stop("Ynodes are not currently allowed before A/C/Z nodes.")

  all.nodes <- c(nodes$A, nodes$C,nodes$Z, nodes$L, nodes$Y)
  if (length(all.nodes) > length(unique(all.nodes))) stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes")
  if (is.null(nodes$Y)) stop("Ynodes cannot be null")
  if (is.null(nodes$Z)) stop("Znodes cannot be null")
  if (is.null(nodes$AC)) stop("Anodes and Cnodes cannot both be null")

  if (min(all.nodes) < ncol(data)) {
    if (!all((min(all.nodes):ncol(data)) %in% all.nodes)) {
      stop("All nodes after the first of A-, C-, Z- L-, or Ynodes must be in A-, C-,Z-, L-, or Ynodes")
    }
  }
  for (reserved.name in c("observation.weights", "Q.kplus1", "Qstar")) {
    #these might cause conflicts
    if (reserved.name %in% names(data)) stop(paste(reserved.name, "is reserved and may not be used as a column name of data"))
  }

  #If gform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(gform)) {
    if (is.character(gform)) {
      if (length(gform) != length(nodes$AC)) stop("length(gform) != length(c(Anodes, Cnodes))")
      for (i in 1:length(gform)) {
        if (LhsVars(gform[i]) != names(data)[nodes$AC[i]]) {
          stop("The LHS of gform[", i, "] should be the name of the ", i, "th A or C node")
        }
        parents <- if(nodes$AC[i] > 1) {
          names(data)[1:(nodes$AC[i]-1)]
        } else {
          NULL
        }
        if (!all(RhsVars(gform[i]) %in% parents)) {
          stop("Some nodes in gform[", i, "] are not parents of ", LhsVars(gform[i]))
        }
        if (any(RhsVars(gform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in gform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    } else {
      if (! is.numeric(gform)) stop("gform should be a character vector or numeric")
      if (nrow(gform) != nrow(data)) stop("if gform is numeric, it should have the same number of rows as data")
      if (ncol(gform) != length(nodes$AC)) stop("if gform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
      if (length(dim(gform)) != 3 || dim(gform)[3] != num.regimes) stop("if gform is numeric, dim[3] should be num.regimes")
      if (!is.null(deterministic.g.function)) stop("if gform is numeric, deterministic.g.function must be NULL")
      if (max(gform, na.rm=T) > 1 || min(gform, na.rm=T) < 0) stop("if gform is numeric, all values should be probabilities")
    }
  }
  #If gform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(qzform)) {
    if (is.character(qzform)) {
      if (length(qzform) != length(nodes$Z)) stop("length(qzform) != length(c(Znodes))")
      for (i in 1:length(qzform)) {
        if (LhsVars(qzform[i]) != names(data)[nodes$Z[i]]) {
          stop("The LHS of qzform[", i, "] should be the name of the ", i, "th Z node")
        }
        parents <- if(nodes$Z[i] > 1) {
          names(data)[1:(nodes$Z[i]-1)]
        } else {
          NULL
        }
        if (!all(RhsVars(qzform[i]) %in% parents)) {
          stop("Some nodes in qzform[", i, "] are not parents of ", LhsVars(qzform[i]))
        }
        if (any(RhsVars(qzform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in qzform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    } else {
      stop('qZ form should be function specifications')
      # if (! is.numeric(qzform)) stop("qzform should be a character vector or numeric")
      # if (nrow(qzform) != nrow(data)) stop("if qzform is numeric, it should have the same number of rows as data")
      # if (ncol(qzform) != length(nodes$Z)) stop("if qzform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
      # if (length(dim(qzform)) != 3 || dim(qzform)[3] != num.regimes.prime) stop("if qzform is numeric, dim[3] should be num.regimes")
      # if (!is.null(deterministic.g.function)) stop("if qzform is numeric, deterministic.g.function must be NULL")
      # if (max(qzform, na.rm=T) > 1 || min(qzform, na.rm=T) < 0) stop("if qzform is numeric, all values should be probabilities")
    }
  }
  #  If gform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(qLform)) {
    if (is.character(qLform)) {
      if (length(qLform) != length(c(nodes$L, nodes$Y))) stop("length(qLform) != length(c(Lnodes, Ynodes))")
      for (i in 1:length(qLform)) {
        if (LhsVars(qLform[i]) != names(data)[sort(c(nodes$L, nodes$Y))[i]]) {
          stop("The LHS of qLform[", i, "] should be the name of the ", i, "th L or Y node")
        }
        parents <- if(sort(c(nodes$L, nodes$Y))[i] > 1) {
          names(data)[1:(sort(c(nodes$L, nodes$Y))[i]-1)]
        } else {
          NULL
        }
        if (!all(RhsVars(qLform[i]) %in% parents)) {
          stop("Some nodes in qLform[", i, "] are not parents of ", LhsVars(qLform[i]))
        }
        if (any(RhsVars(qLform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in qLform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    } else {
      stop('qL form should be function specifications')
      # if (! is.numeric(qLform)) stop("qLform should be a character vector or numeric")
      # if (nrow(qLform) != nrow(data)) stop("if qLform is numeric, it should have the same number of rows as data")
      # if (ncol(qLform) != length(nodes$LY)) stop("if qLform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
      # if (length(dim(qLform)) != 3 || dim(qLform)[3] != num.regimes) stop("qLform cannot be numeric")
      # if (!is.null(deterministic.g.function)) stop("if qLform is numeric, deterministic.g.function must be NULL")
      # if (max(qLform, na.rm=T) > 1 || min(qLform, na.rm=T) < 0) stop("if qLform is numeric, all values should be probabilities")
    }
  }

  #If QLform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(QLform)) {
    if (! is.character(QLform)) stop("QLform should be a character vector")
    if (length(QLform) != length(nodes$LY)) {
      stop("length of Qform is not equal to number of L/Y nodes")
    }
    for (i in 1:length(QLform)) {
      if (length(names(QLform[i])) == 0) stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
      if (names(QLform[i]) != names(data)[nodes$LY[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
      if (QLform[i] != "IDENTITY") {
        #This is only meant to be used in a call by ltmle:::EstimateVariance
        if (LhsVars(QLform[i]) != "Q.kplus1") stop("LHS of each Qform should be Q.kplus1")
        parents <- names(data)[1:(nodes$LY[i]-1)]
        if (!all(RhsVars(QLform[i]) %in% parents)) {
          stop("Some nodes in QLform[", i, "] are not parents of ", names(QLform[i]))
        }
        if (any(RhsVars(QLform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in QLform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    }
  }
  #If QLform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(QZform)) {
    if (! is.character(QZform)) stop("QZform should be a character vector")
    if (length(QZform) != length(nodes$Z)) {
      stop("length of Qform is not equal to number of L/Y nodes")
    }
    for (i in 1:length(QZform)) {
      if (length(names(QZform[i])) == 0) stop("Each element of QZform must be named. The name must match the name of the corresponding L/Y node in data.")
      if (names(QZform[i]) != names(data)[nodes$Z[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
      if (QZform[i] != "IDENTITY") {
        #This is only meant to be used in a call by ltmle:::EstimateVariance
        if (LhsVars(QZform[i]) != "Q.kplus1") stop("LHS of each QZform should be Q.kplus1")
        parents <- names(data)[1:(nodes$Z[i]-1)]
        if (!all(RhsVars(QZform[i]) %in% parents)) {
          stop("Some nodes in QZform[", i, "] are not parents of ", names(QZform[i]))
        }
        if (any(RhsVars(QZform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in QZform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    }
  }

  if (length(gbounds) != 2) stop("gbounds should have length 2")
  if (! (is.null(deterministic.g.function) || is.function(deterministic.g.function))) {
    stop("deterministic.g.function should be a function or NULL")
  }

  if (! all(unlist(data[, nodes$A]) %in% c(0, 1, NA))) stop("in data, all Anodes should be binary")
  #note: Cnodes are checked in ConvertCensoringNodes

  all.Y <- unlist(data[, nodes$Y])

  binaryOutcome <- all(all.Y %in% c(0, 1, NA))

  if (binaryOutcome) {
    if (is.null(survivalOutcome)) {
      if (length(nodes$Y) == 1) {
        survivalOutcome <- FALSE #doesn't matter
      } else {
        stop("All Ynodes are 0, 1, or NA; the outcome is treated as binary. The 'survivalOutcome' argument must be specified if there are multiple Ynodes.")
      }
    }
    if (!is.null(Yrange) && !is_equivalent_to(Yrange, c(0L, 1L))) {
      stop("All Ynodes are 0, 1, or NA, but Yrange is something other than NULL or c(0, 1)")
    }
  } else {
    if (is.null(survivalOutcome)) {
      survivalOutcome <- FALSE
    }
    if (survivalOutcome) {
      stop("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA")
    }
  }

  uncensored.array <- CalcUncensoredMatrix(data, nodes$C)

  for (i in nodes$Y) {
    uncensored <- IsUncensored(uncensored.array, nodes$C, cur.node=i)
    deterministic <- IsDeterministic(data, cur.node=i, deterministic.Q.function=NULL, nodes, called.from.estimate.g=FALSE, survivalOutcome)$is.deterministic #pass deterministic.Q.function=NULL so we're only picking up deaths (if surivalOutcome=FALSE, deterministic will all be FALSE)
    if (anyNA(data[deterministic, i]) || ! all(data[deterministic, i] == 1)) stop("For survival outcomes, once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1.")
    if (anyNA(data[uncensored, i])) stop("Ynodes may not be NA except after censoring")
  }

  stopifnot(num.regimes == nrow(summary.measures))
  if (!all(regimes %in% c(0, 1, NA))) stop("all regimes should be binary")

  if (!all(regimes.prime %in% c(0, 1, NA))) stop("all regimes.prime should be binary")

  for (i in seq_along(nodes$A)) {
    cur.node <- nodes$A[i]
    uncensored <- IsUncensored(uncensored.array, nodes$C, cur.node=i)
    deterministic <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic
    if (anyNA(regimes[uncensored & !deterministic, i, ])) {
      stop("NA in regimes/abar not allowed (except after censoring/death)")
    }
  }

  if ((length(dim(summary.measures)) != 3) || !all.equal(dim(summary.measures)[c(1, 3)], c(num.regimes, length(final.Ynodes)))) stop("summary.measures should be an array with dimensions num.regimes x num.summary.measures x num.final.Ynodes")
  if (class(working.msm) != "character") stop("class(working.msm) must be 'character'")
  if (LhsVars(working.msm) != "Y") stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
  if (!is.vector(observation.weights) || length(observation.weights) != nrow(data) || anyNA(observation.weights) || any(observation.weights < 0) || max(observation.weights) == 0) stop("observation.weights must be NULL or a vector of length nrow(data) with no NAs, no negative values, and at least one positive value")

  if (!IC.variance.only) {
    if (!binaryOutcome) stop("IC.variance.only=FALSE not currently compatible with non binary outcomes")
    if (!is.null(deterministic.Q.function)) stop("IC.variance.only=FALSE not currently compatible with deterministic.Q.function")
    #if (gcomp) stop("IC.variance.only=FALSE not currently compatible with gcomp")
    if (stratify) stop("IC.variance.only=FALSE not currently compatible with stratify=TRUE")
  }
  return(list(survivalOutcome=survivalOutcome, binaryOutcome=binaryOutcome, uncensored=uncensored.array))
}

