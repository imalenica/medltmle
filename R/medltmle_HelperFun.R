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
#' @return Returns vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node.
#' Note that if the current node is first C node, it will return that all samples are uncensored.
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
  } else if (!is.null(nodes$D) && any(nodes$D < cur.node)) {

    #Add option to split C into death and censoring due to other reasons.
    last.Dnode <- max(nodes$D[nodes$D < cur.node])
    is.deterministic <- data[, last.Dnode] %in% TRUE
  } else{
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

  #TO DO: Add D option as well
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
#' @return Returns a logical matrix - n x numCnodes = is uncensored up to and including Cnode[i]

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
#' Get the LY nodes without "blocks" of L/Y nodes uninterrupted by A/C/Z nodes.
#'
#' @param data Original wide-form data.
#' @param nodes List object containing A,C,L,Y,D,W2,AC,ACD,ACZ,Z node index.
#' @param check.Qform Logical variable indicating whether the corresponding Q forms need to be edited.
#' @param Qform If \code{check.Qform} is set to \code{TRUE}, specify the Q forms.
#'
#' @return Returns LY nodes that contain A/C/Z nodes between them.
#'

CreateLYNodes <- function(data, nodes, check.Qform, Qform) {

  LYnodes <- sort(c(nodes$W2, nodes$L, nodes$Y))
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
#' @param data Available data in a wide format.
#' @param Anodes Character variable with exposure nodes.
#' @param Cnodes Character variable with censoring nodes.
#' @param Lnodes Character variable with covariate nodes.
#' @param Ynodes Character variable with outcome nodes.
#' @param Znodes Character variable with mediator nodes.
#' @param Dnodes Character variable with death as censoring nodes.
#' @param W2nodes Character variable with baseline nodes after treatment that need to be fluctuated during the TMLE procedure.
#'
#' @return Returns all nodes in a format necessary for further functions.
#'

CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes,Znodes=NULL,Dnodes=NULL,W2nodes=NULL) {

  Anodes <- NodeToIndex(data, Anodes)
  Cnodes <- NodeToIndex(data, Cnodes)
  Lnodes <- NodeToIndex(data, Lnodes)
  Ynodes <- NodeToIndex(data, Ynodes)
  Dnodes <- NodeToIndex(data, Dnodes)
  W2nodes <- NodeToIndex(data, W2nodes)

  nodes <- list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, D=Dnodes, W2=W2nodes, AC=sort(c(Anodes, Cnodes)), ACD=sort(c(Anodes, Cnodes, Dnodes)), LW2=sort(c(Lnodes, W2nodes)))

  if(!is.null(Znodes)){
    Znodes <- NodeToIndex(data, Znodes)
    nodes$Z <- Znodes
  }

  nodes$ACZ<-sort(c(Anodes, Cnodes, Znodes))
  nodes$baseline <- seq(1, min(c(nodes$A, nodes$L, nodes$C, nodes$Y,nodes$Z)) - 1)
  #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
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
#' @param SL.library List containing two elements: Super Learner algorithms for Q estimation, and Super Learner algorithms for g estimation.
#' @param estimate.type Q or g library.
#'
#' @return Returns Super Learner library for Q or g.
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

      #If not already set to NA, set to NA for all censored samples (including Y and D).
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
#' Check if A as in data equals the regime.
#'
#' @param data TO DO
#' @param regimes TO DO
#' @param Anodes Index of A nodes. Can be easily derived from \code{all.nodes$A}.
#'
#' @return Returns a logical array, n x num.Anodes x num.regimes = follows regime[j] up to and including Anode[i]. NA will count as a match as well.
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

  #If deterministic.g.function is already set to NULL, return FALSE for all.
  if (is.null(deterministic.g.function)) return(default)
  #TO DO: put this in a try-catch
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

################################
# deterministic.g.function
################################

#Note: needs better implementation.

deterministic.g.function <- function(data, current.node, nodes) {

    if (names(data)[current.node] == "A1") {
       det <- (data$L1 < -2 | data$L1 > 2) & !is.na(data$L1)
       prob1 <- ((data$L1 < -2) * 1 + (data$L1 > 2) * 0.1)[det]
     } else if (names(data)[current.node] == "A2") {
       det <- data$A1 == 1 & !is.na(data$A1)
       prob1 <- 1
      } else if (names(data[current.node]) %in% c("C1", "C2")){
       return(NULL)  #this returns the default of no deterministic links
       #note that it is not necessary to specify that prior censoring indicates future censoring
     } else {
       stop("unexpected current.node")
     }
     return(list(is.deterministic=det, prob1=prob1))
   }

################################
# SetSeedIfRegressionTesting
################################

#if comparing outputs between different versions of medltmle, we need to sync random numbers
#before calling SuperLearner or FixScoreEquation since these use random numbers

SetSeedIfRegressionTesting <- function() {

  seed <- Sys.getenv("MEDLTMLE.REGRESSION.TESTING.SEED")
  stopifnot(length(seed) == 1)

  if (seed != "") {
    seed <- as.numeric(seed)
    stopifnot(is.finite(seed))
    set.seed(seed)
  }

  invisible(NULL)
}

################################
# IsBinary
################################

IsBinary <- function(mat) {
  identical(mat, as.numeric(as.logical(mat)))
}

################################
# sseq
################################

#like seq, but returns integer(0) if from > to   (always increments by 1)
sseq <- function(from, to) {
  if (from > to) return(integer(0))
  seq(from, to)
}

################################
# timeOrder_baseline
################################

#' timeOrder_baseline
#'
#' Function to group baseline covariates according to time ordering, in case there is time dependence.
#'
#' @param data data.frame object containing all the baseline covariates and possibly an exposure.
#' @param A Numeric value indicating the column number of exposure. Must be within (1,ncol(data)) range.
#'
#' @return Returns data with set exposure A and ordered baseline covariates. All baseline covariates before A
#' are W1, whereas all after A are W2 covariates. This part is necessary later for TMLE fluctuation.
#'
#' @export timeOrder_baseline

timeOrder_baseline<-function(data, A){

  #How many unique groups of covariates:
  num_uniqe<-length(unique(unlist(lapply(strsplit(names(data),split='[.]'),`[[`, 1))))

  #First, get ordering and possibly subordering of A:
  if(!grepl('[.]', names(data)[A])){

    ord_A<-unlist(strsplit(unlist(strsplit(names(data)[A],"B"))[2],'[.]'))[1]

  }else{

    ord_A<-unlist(strsplit(unlist(strsplit(names(data)[A],"B"))[2],'[.]'))[1]
    subord_A<- unlist(strsplit(names(data)[A],'[.]'))[2]

  }

  if(ord_A==num_uniqe){

    #A belongs to the last group; assign all baseline covariates to W1, none in W2.
    names(data)[A]<-"A"
    names(data)[-A]<-paste(names(data)[-A],"W1",sep=".")

  }else{

    #A is not part of the last group; assign all baseline covariates before A to W1, and all after A to W2.
    names(data)[A]<-"A"

    #Get order for all baseline covariates:
    ord<-sapply(strsplit(sapply(strsplit(names(data[,-A]),"B"), "[[", 2),'[.]'),"[[", 1)
    data_temp<-data[,-A]

    for(i in 1:length(ord)){

      if(ord[i]<ord_A){

        names(data_temp)[i]<-paste(names(data_temp)[i],"W1",sep=".")

      }else if(ord[i]>ord_A){

        names(data_temp)[i]<-paste(names(data_temp)[i],"W2",sep=".")

      }else{
        #A and baseline covariate have the same order.Check suborder.

        subord_cov<- unlist(strsplit(names(data_temp)[i],'[.]'))[2]

        if(subord_cov<subord_A){
          names(data_temp)[i]<-paste(names(data_temp)[i],"W1",sep=".")
        }else{
          names(data_temp)[i]<-paste(names(data_temp)[i],"W2",sep=".")
        }
      }
    }
    data<-cbind.data.frame(data_temp[,1:(as.numeric(ord_A)+1)], A=data$A, data_temp[,(as.numeric(ord_A)+2):ncol(data_temp)])
  }

  return(data)

}

################################
# GetWarningsToSuppress
################################

GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                            "prediction from a rank-deficient fit may be misleading",
                            "non-integer #successes in a binomial glm!",
                            "the matrix is either rank-deficient or indefinite")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

################################
# dataConvert()
################################

#' dataConvert
#'
#' Function to convert wide to long, or long to wide dataframe.
#'
#' @param data Dataframe object in a wide format.
#' @param nBCov number of baseline covariates. They must the listed before C,A,Z,LA,LZ,Y,D nodes.
#' @param type Which format to convert the current data to, long or wide.
#'
#' @return Dataframe object of the original data in a long or wide format.
#'
#' @export dataConvert

dataConvert<-function(data,nBCov,type){

  if(type=="long"){

    #Adds ID identifier as needed by stremr
    data<-cbind.data.frame(row.names(data),data)
    names(data)[1]<-"ID"

    newdata<-reshape(data, idvar="ID", varying=(nBCov+2):length(data), sep="_", direction="long")
    newdata$ID<-as.numeric(levels(newdata$ID))[newdata$ID]
    newdata<-cbind.data.frame(newdata$ID,newdata[,grep("time", colnames(newdata))],newdata[,2:(grep("t", colnames(newdata))-1)],newdata[,(grep("t", colnames(newdata))+1):ncol(newdata)])
    names(newdata)[1:2]<-c("ID","t")

    #Remove all NA from C node:
    newdata<-newdata[!is.na(newdata$C),]

  }else if(type=="wide"){

    #Just use Oleg's convert.to.wide() function

  }else{
    print("Can convert data to either long or wide format, no other options.")
  }

  return(newdata)

}

################################
# summary_medltmle
################################

#' summary_medltmle
#'
#' Function to summarize Natural Indirect Effect and Natural Direct Effect parameters estimated by \code{medltmle}.
#'
#' @param nie1 First part of the treatment/mediator specific mean returned by \code{medltmle} for the Natural Indirect Effect. It should correspond to the output generated by \code{medltmle} with abar=1, abar.prime=1 for NIE.
#' @param nie2 Second part of the treatment/mediator specific mean returned by \code{medltmle} for the Natural Indirect Effect.  It should correspond to the output generated by \code{medltmle} with abar=1, abar.prime=0 for NIE.
#' @param nde1 First part of the treatment/mediator specific mean returned by \code{medltmle} for the Natural Direct Effect.  It should correspond to the output generated by \code{medltmle} with abar=1, abar=0 for NDE.
#' @param nde2 Second part of the treatment/mediator specific mean returned by \code{medltmle} for the Natural Direct Effect.  It should correspond to the output generated by \code{medltmle} with abar=0, abar.prime=0 for NDE.
#' @param type Options are "NE" or "SE" for non-data adaptive and data-adaptive treatment/mediator specific mean returned by \code{medltme}.
#'
#' @return Returns summary statistics corresponding to Natural Indirect Effect and Natural Direct Effect parameters.
#'
#' @export summary_medltmle

summary_medltmle<-function(nie1,nie2,nde1,nde2,type="NE"){

  #Natural Indirect Effect
  #TMLE or GComp
  param_nie_tmle<-nie1$estimates[1]-nie2$estimates[1]
  var_nie_tmle<-var(nie1$IC$tmle-nie2$IC$tmle, na.rm = TRUE)/sum(!is.na(nie1$IC$tmle))
  se_nie_tmle<-sqrt(var_nie_tmle)

  CI_nie_tmle_l<-param_nie_tmle-1.96*se_nie_tmle
  CI_nie_tmle_u<-param_nie_tmle+1.96*se_nie_tmle

  #IPTW
  param_nie_iptw<-nie1$estimates[2]-nie2$estimates[2]
  var_nie_iptw<-var(nie1$IC$iptw-nie2$IC$iptw, na.rm=TRUE)/sum(!is.na(nie1$IC$iptw))
  se_nie_iptw<-sqrt(var_nie_iptw)

  CI_nie_iptw_l<-param_nie_iptw-1.96*se_nie_iptw
  CI_nie_iptw_u<-param_nie_iptw+1.96*se_nie_iptw

  res_nie_tmle<-data.frame(param_nie_tmle,var_nie_tmle,se_nie_tmle,CI_nie_tmle_l,CI_nie_tmle_u)
  names(res_nie_tmle)<-c("NIE", "Var NIE", "SE NIE", "CI lower", "CI upper")
  res_nie_iptw<-data.frame(param_nie_iptw,var_nie_iptw,se_nie_iptw,CI_nie_iptw_l,CI_nie_iptw_u)
  names(res_nie_iptw)<-c("NIE", "Var NIE", "SE NIE", "CI lower", "CI upper")

  res_nie<-rbind.data.frame(res_nie_tmle,res_nie_iptw)

  #Natural Direct Effect
  #TMLE
  param_nde_tmle<-nde1$estimates[1]-nde2$estimates[1]
  var_nde_tmle<-var(nde1$IC$tmle-nde2$IC$tmle, na.rm=TRUE)/sum(!is.na(nde1$IC$tmle))
  se_nde_tmle<-sqrt(var_nde_tmle)

  CI_nde_tmle_l<-param_nde_tmle-1.96*se_nde_tmle
  CI_nde_tmle_u<-param_nde_tmle+1.96*se_nde_tmle

  #IPTW
  param_nde_iptw<-nde1$estimates[2]-nde2$estimates[2]
  var_nde_iptw<-var(nde1$IC$iptw-nde2$IC$iptw, na.rm=TRUE)/sum(!is.na(nde1$IC$iptw))
  se_nde_iptw<-sqrt(var_nde_iptw)

  CI_nde_iptw_l<-param_nde_iptw-1.96*se_nde_iptw
  CI_nde_iptw_u<-param_nde_iptw+1.96*se_nde_iptw

  res_nde_tmle<-data.frame(param_nde_tmle,var_nde_tmle,se_nde_tmle,CI_nde_tmle_l,CI_nde_tmle_u)
  names(res_nde_tmle)<-c("NDE", "Var NDE", "SE NDE", "CI lower", "CI upper")
  res_nde_iptw<-data.frame(param_nde_iptw,var_nde_iptw,se_nde_iptw,CI_nde_iptw_l,CI_nde_iptw_u)
  names(res_nde_iptw)<-c("NDE", "Var NDE", "SE NDE", "CI lower", "CI upper")

  res_nde<-rbind.data.frame(res_nde_tmle,res_nde_iptw)

  #Total effect
  #TMLE
  param_tmle<-param_nie_tmle+param_nde_tmle
  var_tmle<-var(nie1$IC$tmle-nde2$IC$tmle, na.rm=TRUE)/sum(!is.na(nie1$IC$tmle))
  se_tmle<-sqrt(var_tmle)

  CI_tmle_l<-param_tmle-1.96*se_tmle
  CI_tmle_u<-param_tmle+1.96*se_tmle

  #IPTW
  param_iptw<-param_nie_iptw+param_nde_iptw
  var_iptw<-var(nie1$IC$iptw-nde2$IC$iptw, na.rm=TRUE)/sum(!is.na(nie1$IC$iptw))
  se_iptw<-sqrt(var_iptw)

  CI_iptw_l<-param_iptw-1.96*se_iptw
  CI_iptw_u<-param_iptw+1.96*se_iptw

  res_tmle<-data.frame(param_tmle,var_tmle,se_tmle,CI_tmle_l,CI_tmle_u)
  names(res_tmle)<-c("NE", "Var NE", "SE NE", "CI lower", "CI upper")
  res_iptw<-data.frame(param_iptw,var_iptw,se_iptw,CI_iptw_l,CI_iptw_u)
  names(res_iptw)<-c("NE", "Var NE", "SE NE", "CI lower", "CI upper")

  res_ne<-rbind.data.frame(res_tmle,res_iptw)

  return(list(NDE=res_nde,NIE=res_nie,NE=res_ne))

}
