######################################################################
## These functions are minor modifications or directly copied from the
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.


#'
#' get coefficients or make coefficient predictions from and "FuncompCGL" object
#'
#' @description
#'          Computes the coefficients at the requested values for \code{lam} from a fitted
#'          \code{\link{FuncompCGL}} object.
#'
#' @param object fitted \code{\link{FuncompCGL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required. Default is the
#'          entire sequence used to create the model.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @details
#'          \code{s} is the new vector at which predictions are requested. If \code{s} is not in the
#'          lambda sequence used for fitting the model, the \code{coef} function will use linear
#'          interpolation to make predictions. The new values are interpolated using a fraction of
#'          coefficients from both left and right \code{lam} indices.
#'
#' @return
#' The coefficients at the requested values for \code{s}.
#'
#' @export
#'


coef.FuncompCGL <- function(object, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1)
    {
      beta = beta[, lamlist$left, drop=FALSE] *  (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    rownames(seq(s))
  }

  if(is.numeric(object$ref)) {
    # checked dimension after 'apply'
    beta <- apply(beta, 2, function(x, k, p1, ref)
                  c(x[ifelse(ref==1, 0, 1):((ref-1)*k)], -colSums(matrix(x[1:p1], byrow = TRUE, ncol = k)), x[((ref-1)*k+1):(length(x))]),
                  ref = object$ref, k = as.list(object$call)$k, p1 = ncol(object$Z) )
  }

  return(beta)
}


#'
#' get coefficients or make coefficient predictions from and "compCL" object
#'
#' @description
#'          Computes the coefficients at the requested values for \code{lam} from a fitted
#'          \code{\link{compCL}} object.
#'
#' @param object fitted \code{\link{compCL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required. Default is the
#'          entire sequence used to create the model.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @details
#'          \code{s} is the new vector at which predictions are requested. If \code{s} is not in the
#'          lambda sequence used for fitting the model, the \code{coef} function will use linear
#'          interpolation to make predictions. The new values are interpolated using a fraction of
#'          coefficients from both left and right \code{lam} indices.
#'
#' @return
#' The coefficients at the requested values for \code{s}.
#'
#' @export
#'
#'

# object <- compCL.object
# s <- 1
# s <- c(1, 0.5, 0.1)
# coef(compCL.object, s = 1)
# coef(compCL.object, s = c(1, 0.5, 0.1))
#coef(compCL.object)


# TO add baseline to coef.compCL
coef.compCL <- function(object, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1) {
      beta = beta[, lamlist$left, drop=FALSE] *  (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    colnames(beta) <- paste(seq(along = s))
  }
  return(beta)
}



#'
#' make predictions from a "FuncompCGL" object
#'
#' @description
#'  predicts fitted values and class labels from a fitted \code{\link{FuncompCGL}} object.
#'
#' @param object fitted \code{\link{FuncompCGL}} object.
#'
#' @param newx matrix of new values for x at which predictions are to be made. Must be a matrix.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the entire sequence used to create the model.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @details
#'        \code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda
#'        sequence used for fitting the model, the \code{predict} function will use linear interpolation
#'        to make predictions. The new values are interpolated using a fraction of predicted values from
#'        both left and right \code{lam} indices.
#'
#' @return
#' prediction values at the requested values for \code{s}.
#'
#' @export
#'
#' @importFrom methods cbind2
#'


predict.FuncompCGL <- function(object, newx, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1)
    {
      beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    #dimnames(beta) <- list(vnames, paste(seq(along = s)))
  }
  if (is.null(dim(newx))) newx = matrix(newx, nrow = 1)
  fitting <- cbind2(newx, 1) %*% beta #as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  return(fitting)
}


#'
#' make predictions from a "compCL" object
#'
#' @description
#'  predicts fitted values and class labels from a fitted \code{\link{compCL}} object.
#'
#' @param object fitted \code{\link{compCL}} object.
#'
#' @param Znew,Zcnew Znew is the matrix of new values for log composition and Zcnew is the matrix for confounding
#'                   matrix. Joint of Znew and Zcnew is the design matrix at which predictions are to be made.
## at which predictions are to be made. Must be a matrix.
#'                   Default value for Zcnew is \code{NULL}.
#'
## @param  matrix for new values for Zc.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the entire sequence used to create the model.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @details
#'        \code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda
#'        sequence used for fitting the model, the \code{predict} function will use linear interpolation
#'        to make predictions. The new values are interpolated using a fraction of predicted values from
#'        both left and right \code{lam} indices.
#'
#' @return
#' prediction values at the requested values for \code{s}.
#'
#' @export
#'

# object <- compCL.object
# Znew <- Z[1, ]
# Zcnew = NULL
# s <- 1
# s <- c(1, 0.5, 0.1)
# Znew <- Z[1:5,]
# predict(compCL.object, Znew = Z[1:5, ], s = 1)
# predict(compCL.object, Znew = Z[1:5, ])

predict.compCL <- function(object, Znew, Zcnew = NULL,  s = NULL, ...) {
  #print(object$beta)
  beta <- object$beta
  nvars <- dim(beta)[1] - 1
  if( is.null(dim(Znew)) ) Znew = matrix(Znew, nrow = 1)
  if( !is.null(Zcnew) && is.null(dim(Zcnew)) ) Zcnew = matrix(Zcnew, nrow = 1)
  new <- cbind(Znew, Zcnew)
  if( nvars != dim(new)[2] ) stop("numbers of variables in data is not consistant with estimated coefficients")

  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1) {
      beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    colnames(beta) <- paste(seq(along = s))
  }

  #colnames(beta) <- paste(seq(along = 1:dim(beta)[2]))
  fitting <- cbind2(new, 1) %*% beta   #as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  return(fitting)
}



#'
#' print a "FuncompCGL" object
#'
#' @description
#' Print the nonzero counts for composition varaibles at each lambda along the FuncompCGL path.
#'
#' @param x fitted \code{\link{FuncompCGL}} object.
#'
#' @param digits significant digits in printout.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @return a two-column matrix, the first columns is the number of nonzero group counts and
#'         the second column is Lambda.
#'
#' @export
#'

print.FuncompCGL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  #print(cbind(Df = x$df, Lam = signif(x$lam, digits)))
  show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
  colnames(show) <- c("DF", "Lam")
  rownames(show) <- seq(nrow(show))
  print(show)
}



#'
#' print a "compCL" object
#'
#' @description
#' Print the nonzero counts for composition varaibles at each lambda along the compCL path.
#'
#' @param x fitted \code{\link{compCL}} object.
#'
#' @param digits significant digits in printout.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @return a two-column matrix, the first columns is the number of nonzero group counts and
#'         the second column is Lambda.
#'
#' @export
#'

#print(compCL.object)

print.compCL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
  colnames(show) <- c("Df", "Lam")
  rownames(show) <- paste0("L", seq(nrow(show)))
  print(show)
  #print(cbind(Df = object$df, Lam = signif(object$lam, digits)))
}


# print.compCL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
#   cat("\nCall: ", deparse(x$call), "\n\n")
#   show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
#   colnames(show) <- c("DF", "Lam")
#   rownames(show) <- seq(nrow(show))
#   print(show)
# }

#'
#' get coefficients or make coefficient predictions from a cross-validation "cv.FuncompCGL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a cross-validated
#' \code{FuncompCGL} model, using the stored \code{"FuncompCGL.fit"} object,
#' and the optimal value chosen for \code{lam} and\code{k} (degree of freedom of basis).
#' Questionable!!!
#'
#' @param object fitted \code{\link{cv.FuncompCGL}} object.
#'
#' @param trim logical, used the trimmed result or not.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the value \code{s="lam.min"} stored on the CV \code{object}, it is the
#'          value of \code{lam} for optimal \code{k} such that gives minimum mean cross-validated error.
#'          If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param k value of df of basis at which predictions are required. Default is \code{k = "k.min"},
#'          it is the value of \code{k} that gives minimum mean cross-validated error. If \code{k}
#'          is a scaler, it is taken as the value of \code{k} to be used and it must be stored in
#'          \code{\link{cv.FuncompCGL}} object.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @return The coefficients at the requested values for \code{lam} and \code{k}.
#'
#' @export
#'


# s <- "lam.1se"
# object <- cv.m1
# trim = FALSE
coef.cv.FuncompCGL <- function(object, trim = FALSE, s = c("lam.min", "lam.1se"), k = NULL, ...) {
  trim <- ifelse(trim, "Ttrim", "Ftrim")

  if (is.numeric(s)) {
    lam_se <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object[[trim]][[s]][1]
    k_se <- object[[trim]][[s]][2]
  } else stop("Invalid form for s")

  if(is.numeric(k)){
    if(! (k %in% as.integer(names(object$Funcomp.CGL.fit)) )) stop("K is out of range")
    k_se <- k
  }
  k_se <- as.character(k_se)
  cv.fit <- object$Funcomp.CGL.fit[[k_se]]

  #??????????????????????
  # pass arugment from inner function ??????????
  cv.fit$call <- as.list(cv.fit$call)
  cv.fit$call$k <- as.numeric(k_se)
  cv.fit$call$ref <-  as.numeric(cv.fit$ref)
  cv.fit$call
  #??????????????????????
  coef( cv.fit, s = lam_se) #coef.FuncompCGL

}


##
## make predictions from a "cv.FuncompCGL" object.
## This function makes prediction from a cross-validated \code{FuncompCGL} model,
## using the stored \code{FuncompCGL.fit} object and the optimal value chose for \code{lambda}.
##
## @inheritParams
##
## @export

# predict.cv.FuncompCGL <- function(object, Znew, Zcnew = NULL,
#                                   s = c("lam.1se", "lam.min"), k = NULL,
#                                   trim = FALSE, ...) {
#
#
#   trim <- ifelse(trim, "Ttrim", "Ftrim")
#
#   if (is.numeric(s)) {
#     lam_se <- s
#   } else if (is.character(s)) {
#     s <- match.arg(s)
#     lam_se <- object[[trim]][[s]][1]
#     k_se <- object[[trim]][[s]][2]
#   } else stop("Invalid form for s")
#
#   if(is.numeric(k)){
#     if(! (k %in% as.integer(names(object$Funcomp.CGL.fit)) )) stop("K is out of range")
#     k_se <- k
#   }
#   k_se <- as.character(k_se)
#   #lam_loc = match(lam_se, object$lam)
#   param_list <- as.list(object$Funcomp.CGL.fit[[k_se]]$call)
#   param_list$outer_maxiter <- 0
#   param_list$X <- as.data.frame(Znew)
#   param_list$T.name <- colnames(Znew)[2]
#   param_list$ID.name <- colnames(Znew)[1]
#   param_list$k <- as.numeric(k_se)
#   Znew <- FuncompCGL()
#   beta <- predict(object$Funcomp.CGL.fit[[k_se]], newx = cbind(Znew, Zcnew),
#                   s = lam_se)
#
# }





#'
#' get coefficients or make coefficient predictions from a cross-validation "cv.compCL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a cross-validated
#' \code{compCL} model, using the stored \code{"compCL.fit"} object,
#' and the optimal value chosen for \code{lam}.
#'
#' @param object fitted \code{\link{cv.compCL}} object.
#'
#' @param trim logical, used the trimmed result or not. Default is FASLE.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'         Default is the value \code{s="lam.min"} stored on the CV \code{object}, it is the
#'         optimal value of \code{lam} that gives minimum.
#'         Alternatively \code{s="lambda.min"} can be used,  it is the largest value of \code{lam}
#'         such that error is within 1 standard error of the minimum.
#'         If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @return
#' The coefficients at the requested values for \code{lam}.
#'
#' @export
#'

# object <- cvm
# trim <- FALSE
# s <- 1
# coef(cvm, trim = FALSE, s = "lam.min")
#coef(cvm, trim = TRUE, s = "lam.min")

coef.cv.compCL <- function(object, trim = FALSE, s = c("lam.min", "lam.1se" ),...) {
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  object_use <- object[[trim]]
  if (is.numeric(s)) {
    lam <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    # switch(s,
    #           "lam.min" = "lambda.min",
    #           "lam.1se" = "lambda.1se")
    lam <- object_use[[s]]
  } else {
    stop("Invalid form for s")
  }
  select <- list(beta = coef(object$compCL.fit, s = lam, ...), lam = lam)
  return(select)
}



#' @title
#' plot the cross-validation curve produced by cv.FuncompCGL
#'
#' @description
#' Plots the cross-validation curve for different degree of freedom \code{k},
#' and upper and lower standard deviation curves,
#' as a function of the \code{lambda} values used.
#'
#' @param x fitted \code{"cv.FuncompCGL"}.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
#' @param trim Logic value, determining use trimmed result stored in \code{cv.FumcompCGL} or not. Default
#' @param k_list a vector, determining result of which \code{k} (degree freedom) to plot.
#'               If not provided, cross-validation curve for k that associate with \code{lambda.min} (blue)
#'               and \code{lambda.1se} (red) are plotted.
#' @param \dots Other graphical parameters to plot
#'
#' @details
#' A plot is prodoced, and nothing is returned.
#'
#' @export
#' @import graphics

plot.cv.FuncompCGL <- function(x, xlab = c("log", "lambda"),
                               trim = FALSE,
                               k_list,
                               ...) {
  #cat("1 \r\n")
  cvobj = x
  xlab <- match.arg(xlab)
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  cvobj_use <- cvobj[[trim]]
  #xlab = "log"
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(cvobj$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(cvobj$lam))
         })
  #cat(xlab, "\r\n")
  k_listall <- as.numeric(names(cvobj$Funcomp.CGL.fit))
  rule_min <- cvobj_use$lam.min
  rule_1se <- cvobj_use$lam.1se
  if(missing(k_list)) k_list <- c(rule_min['df'], rule_1se['df'])
  if(is.numeric(k_list)) k_list <- k_list[k_list %in% k_listall]
  if(is.character(k_list)) k_list <- switch(k_list,
                                            "lam.1se" = rule_1se['df'],
                                            "lam.min" = rule_min['df'])

  #cat(k_list, "\r\n")
  N_list <- apply(cvobj_use$cvm[k_list- k_listall[1] + 1, ], 1, function(x) length(x[!is.na(x)]))
  plot.args = list(x = xvalue[seq(max(N_list))], xlab = xlab,
                   y = cvobj_use$cvm[k_list[which.max(N_list)] - k_listall[1] + 1, seq(max(N_list))],
                   #y = cvobj_use$cvm[k_list[1], ],
                   ylab = "MEAN-Squared Error", ##
                   type = "n",
                   ylim = range(cvobj_use$cvup[k_list- k_listall[1] + 1, ], cvobj_use$cvlo[k_list- k_listall[1] + 1, ], na.rm = TRUE)
  )
  #cat(plot.args[[1]], "\r\n")
  new.args = list(...) #list(...)
  #cat("2\r\n")
  if(length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)


  for(l in 1:length(k_list)) {
    if(k_list[l] != rule_min['df'] && k_list[l] != rule_1se['df']) {
      # ll <- N_list[l]
      # error.bars(xvalue[seq(ll)],
      #            cvobj_use$cvup[k_list[l] - k_listall[1] + 1, seq(ll)],
      #            cvobj_use$cvlo[k_list[l] - k_listall[1] + 1, seq(ll)],
      #            width=0.01,col="lightgrey")
      #

      error.bars(xvalue,
                 cvobj_use$cvup[k_list[l] - k_listall[1] + 1, ],
                 cvobj_use$cvlo[k_list[l] - k_listall[1] + 1, ],
                 width=0.01,col="lightgrey")
    }
  }


  if(rule_min['df'] %in% k_list) {
    # ll <- N_list[which(k_list == rule_min['df'])[1]]
    # error.bars(xvalue[seq(ll)],
    #            cvobj_use$cvup[rule_min['df'] - k_listall[1] + 1, seq(ll)],
    #            cvobj_use$cvlo[rule_min['df'] - k_listall[1] + 1, seq(ll)],
    #            width=0.01,col="red")


    error.bars(xvalue,
               cvobj_use$cvup[rule_min['df'] - k_listall[1] + 1, ],
               cvobj_use$cvlo[rule_min['df'] - k_listall[1] + 1, ],
               width=0.01,col="red")
  }

  if(rule_1se['df'] %in% k_list){
    # ll <- N_list[which(k_list == rule_1se['df'])[1]]
    # error.bars(xvalue[seq(ll)],
    #            cvobj_use$cvup[rule_1se['df'] - k_listall[1] + 1, seq(ll)],
    #            cvobj_use$cvlo[rule_1se['df'] - k_listall[1] + 1, seq(ll)],
    #            width=0.01,col="blue")


    error.bars(xvalue,
               cvobj_use$cvup[rule_1se['df'] - k_listall[1] + 1, ],
               cvobj_use$cvlo[rule_1se['df'] - k_listall[1] + 1, ],
               width=0.01,col="blue")
  }



  for(l in 1:length(k_list)) {
    if(k_list[l] != rule_min['df'] && k_list[l] != rule_1se['df']) {
      points(xvalue,
             cvobj_use$cvm[k_list[l] - k_listall[1] + 1, ],
             pch=18,
             col="limegreen")
    }
  }

  if(rule_min['df'] %in% k_list) {
    points(xvalue,
           cvobj$Ftrim$cvm[rule_min['df'] - k_listall[1] + 1, ],
           pch=20,
           col="red")
    axis(side=3,at=xvalue,
         labels=paste(cvobj$Funcomp.CGL.fit[[as.character(rule_min['df'])]]$df[1:length(xvalue)]),
         tick=FALSE,
         line=0,
         pos = plot.args$ylim[2],
         col.axis = "red")

  }
  if(rule_1se['df'] %in% k_list){
    points(xvalue,
           cvobj$Ftrim$cvm[rule_1se['df'] - k_listall[1] + 1, ],
           pch=20,
           col="blue")
    axis(side=3,at=xvalue,
         labels=paste(cvobj$Funcomp.CGL.fit[[as.character(rule_1se['df'])]]$df[1:length(xvalue)]),
         tick=FALSE,
         line=0,
         col.axis = "blue")
  }


  abline(v = switch(xlab,
                    "Lambda" = rule_min["lam"],
                    "Log(Lambda)" = log(rule_min["lam"])),
         lty = 3, col = "red"
  )

  abline(v = switch(xlab,
                    "Lambda" = rule_1se["lam"],
                    "Log(Lambda)" = log(rule_1se["lam"])),
         lty = 3, col = "blue"
  )

}
#
#
# coef.cv.compCL <- function(object, trim = FALSE, s = "lam.min",...) {
#   trim <- ifelse(trim, "Ttrim", "Ftrim")
#
#   cv.fit <- object$compCGL.fit
#
#   if (is.numeric(s)) {
#     lam <- s
#   } else if (is.character(s)) {
#     lam <- object[[trim]]$lam.min
#   } else stop("Invalid form for s")
#
#   coef(cv.fit, s = lam)
#
# }







#'
#' get coefficients or make coefficient predictions from a "GIC.compCL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a regulaized fitting
#' \code{compCL} model, using the stored \code{"compCL.fit"} object,
#' and the optimal value chosen for \code{lam}.
#'
#' @param object fitted \code{\link{GIC.compCL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'         Default is the value \code{s="lam.min"} stored on the CV \code{object}, it is the
#'         optimal value of \code{lam} that gives minimum.
#'         Alternatively \code{s="lambda.min"} can be used,  it is the largest value of \code{lam}
#'         such that error is within 1 standard error of the minimum.
#'         If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @return
#' The coefficients at the requested values for \code{lam}.
#'
#' @export
#'

coef.GIC.compCL <- function(object, s = "lam.min",...) {

  if (is.numeric(s)) {
    lam <- s
  } else if (s == "lam.min") {
    lam <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  select <- list(beta = coef(object$compCL.fit, s = lam, ...), lam = lam)
  return(select)
}
#
#
# coef.cv.compCL <- function(object, trim = FALSE, s = "lam.min",...) {
#   trim <- ifelse(trim, "Ttrim", "Ftrim")
#
#   cv.fit <- object$compCGL.fit
#
#   if (is.numeric(s)) {
#     lam <- s
#   } else if (is.character(s)) {
#     lam <- object[[trim]]$lam.min
#   } else stop("Invalid form for s")
#
#   coef(cv.fit, s = lam)
#
# }








#' @title
#' plot coefficients from a \code{"FuncompCGL"} object
#'
#' @description
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{"FuncompCGL"} object.
#'
#' @param x fitted \code{"FuncompCGL"} model.
#' @param p numbers of compositional variables.
#' @param k degree of freedom of basis used.
#' @param ylab What is on the Y-axis, L1-norm or L2-norm (default) of each group of coefficients.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
#' @param \dots Other graphical parameters to plot.
#'
#' @details
#' A plot is prodoced, and nothing is returned.
#'
#' @export
#'


plot.FuncompCGL <- function(x, p, k,
                            ylab = c("L2", "L1"),
                            xlab = c("log", "lambda"),
                            ...) {
  obj <- x
  ylab = match.arg(ylab)
  xlab = match.arg(xlab)
  if( is.numeric(obj$ref) && !is.numeric(as.list(obj$call)$k) ) {
    obj$call <- as.list(obj$call)
    obj$call$k <- k
  }
  beta <- coef(obj) #obj$beta
  include_which <- sort(unique(unlist(apply(beta, 2, Nzero, p = p, k = k))))
  group <- matrix(1:(p*k), nrow = k)
  beta <- beta[group[, include_which], ]
  #cat("1\r\n")
  switch(ylab,
         "L1" = {
           yvalue = apply(beta, 2, function(x, k) {
             value = vector()
             for(l in 1:(length(x)/k)) value[l] <- sum(abs(x[(l*k-k+1):(l*k)]))
             return(value)
           }, k = k)
           ylab = "L1 norm"},
         "L2" = {
           yvalue = apply(beta, 2, function(x, k) {
             value = vector()
             for(l in 1:(length(x)/k)) value[l] <- sqrt(sum((x[(l*k-k+1):(l*k)])^2))
             return(value)
           }, k = k)
           ylab = "L2 norm"
         }
  )

  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(obj$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(obj$lam))
         })


  dotlist=list(...) #list(...)
  type=dotlist$type
  if(is.null(type))
    matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab,type="l", ylim = c(-0.1, max(yvalue)), ...
    )  else matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab
                    , ylim = c(-0.1, max(yvalue)),...
    )

  #text(x = min(xvalue), y = yvalue[,dim(yvalue)[2]], labels = paste(include_which))
  # cvraw <- (drop(y) - cv.m1$fit.preval[, , 1])^2
  # N <- length(y) - apply(is.na(cv.m1$fit.preval[, , 1]), 2, sum)
  # cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  # cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

  at.txt = apply(yvalue, 1, function(x) which(x>0)[1])
  at.index <- rbind(at.txt, sort(at.txt))
  #matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab,type="l", ylim = c(-0.4, max(yvalue)))
  text(x = xvalue[at.txt[order(at.txt)[-((1:floor(length(at.txt)/2))*2)]]],
       y = -0.05, labels = paste(include_which[order(at.txt)][-((1:floor(length(at.txt)/2))*2)]),
       col = rep(c("black", "grey"), length.out = length(include_which) - length((1:floor(length(at.txt)/2))*2)) )
  text(x = xvalue[at.txt[order(at.txt)[((1:floor(length(at.txt)/2))*2)]]],
       y = -0.1, labels = paste(include_which[order(at.txt)][(1:floor(length(at.txt)/2))*2]),
       col = rep(c("black", "grey"), length.out = length((1:floor(length(at.txt)/2))*2)) )



}





#' @title
#' plot the cross-validation curve produced by "cv.compCL" object.
#'
#' @description
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used.
#'
#' @param x fitted \code{"cv.compCL"} object.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
#' @param \dots Other graphical parameters to plot.
#'
#' @details
#' A plot is produced, and nothing is returned.
#'
#' @export
#'


plot.cv.compCL <- function(x, xlab = c("log", "lambda"), ...) {
  cvobj <- x
  xlab <- match.arg(xlab)
  #trim <- ifelse(trim, "Ttrim", "Ftrim")
  cvobj_use <- cvobj[["Ftrim"]]
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(cvobj$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(cvobj$lam))
         })


  plot.args = list(x = xvalue,y = cvobj_use$cvm,
                   ylim = range(cvobj_use$cvupper,cvobj_use$cvlo),
                   xlab=xlab,
                   ylab="MEAN-Squared Error",
                   type="n")
  new.args=list(...)
  if(length(new.args)) plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  error.bars(xvalue, cvobj_use$cvupper, cvobj_use$culo, width=0.01, col="darkgrey")
  points(xvalue,cvobj_use$cvm,pch=20,col="red")
  axis(side=3,at=xvalue,labels=paste(cvobj$compCL.fit$df),tick=FALSE,line=0)
  abline(v=switch(xlab,
                  "Lambda" = cvobj_use$lam.min,
                  "Log(Lambda)" = log(cvobj_use$lam.min)
  ), lty=3)
  abline(v = switch(xlab,
                    "Lambda" = cvobj_use$lam.1se,
                    "Log(Lambda)" = log(cvobj_use$lam.1se)
  ),lty=3)
  invisible()
}



# #' @title
# #' plot the GIC curve produced by "GIC.FuncompCGL" object.
# #'
# #' @description
# #' Plots the GIC curve, as a function of the \code{lambda} values used.
# #'
# #' @param x fitted "\code{GIC.FuncompCGL}" object.
# #' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
# #' @param \dots Other graphical parameters to plot.
# #' @param alpha scaler to penalty model selection.
# #' @details
# #' $$GIC(\eqn{\lambda}) = log(MSE) + Selection * \code{alpha}$$, for normal error.
# #' If include linear constraints, like in log contract and constraints group lasso model,
# #' selection = None zero group - 1; otherwise with consideration of linear constraints,
# #' selection = None zero group.
# #'
# #' BIC:alpha = log(n)*df, GIC:alpha=log(log(n)) /n * log(max(p*df,n)) * df
# #' @export
# #'
#
#
# plot.GIC.FuncompCGL <- function(x, xlab = c("log", "lambda"), alpha, ...) {
#   cvobj <- x
#   xlab <- match.arg(xlab)
#
#   switch(xlab,
#          "lambda" = {
#            xlab = "Lambda"
#            xvalue = drop(cvobj$lam)
#          },
#          "log" = {
#            xlab = "Log(Lambda)"
#            xvalue = log(drop(cvobj$lam))
#          })
#   GIC_val <- matrix(NA, nrow = length(cvobj$fit), length(cvobj$lam) )
#   for(i in length(cvobj$fit)) {
#
#   }
#   plot.args = list(x = xvalue,y = cvobj_use$cvm,
#                    ylim = range(cvobj_use$cvupper,cvobj_use$cvlo),
#                    xlab=xlab,
#                    ylab="GIC",
#                    type="n")
#   new.args=list(...)
#   if(length(new.args)) plot.args[names(new.args)]=new.args
#   do.call("plot",plot.args)
#   error.bars(xvalue, cvobj_use$cvupper, cvobj_use$culo, width=0.01, col="darkgrey")
#   points(xvalue,cvobj_use$cvm,pch=20,col="red")
#   axis(side=3,at=xvalue,labels=paste(cvobj$compCL.fit$df),tick=FALSE,line=0)
#   abline(v=switch(xlab,
#                   "Lambda" = cvobj_use$lam.min,
#                   "Log(Lambda)" = log(cvobj_use$lam.min)
#   ), lty=3)
#   abline(v = switch(xlab,
#                     "Lambda" = cvobj_use$lam.1se,
#                     "Log(Lambda)" = log(cvobj_use$lam.1se)
#   ),lty=3)
#   invisible()
# }
#
# #
