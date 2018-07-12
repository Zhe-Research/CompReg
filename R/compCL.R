#' @title
#' Fits regularization paths for longitudinal compositional data with lasso penalty.
#'
#' @description
#' Fits regularization paths for longitudinal compositional data with lasso penalty at a sequence of regularization parameters lambda.
#'
#'
#'
# @usage
# compCL <- function(y, Z, Zc = NULL, intercept = TRUE, pf = rep(1, times = p),
#                    lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
#                    dfmax = p, pfmax = min(dfmax * 1.5, p),
#                    u = 1,mu_ratio = 1.01, tol = 1e-10,
#                    outer_maxiter = 1e+08, outer_eps = 1e-8,
#                    inner_maxiter = 1e+3, inner_eps = 1e-6)
#
#'
#'
#' @param y a vector of response variable with length n.
#'
#' @param Z a \eqn{n*p} matrix after taking log transformation on compositional data.
#'
#' @param Zc a design matrix of other covariates considered. Default is \code{NULL}.
#'
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#'
#' @param pf penalty factor, a vector in length of p. Separate penalty weights can be applied to each coefficience
#'           \eqn{\beta} for composition variates to allow differential shrinkage. Can be 0 for some \eqn{\beta}'s,
#'           which implies no shrinkage, and results in that composition always being included in the model.
#'           Default value for each entry is the 1.
#'
#' @param u \code{u} is the inital value for penalty parameter of augmented Lanrange method.
#'                   Default value is 1.
#'
#' @inheritParams FuncompCGL
#'
#' @return An object with S3 calss \code{\link{compCL}}
#' \item{beta}{a matrix of coefficients for \code{cbind{Z, Zc, 1_n}}, with \code{nlam} rows.}
#' \item{lam}{the actual sequence of \code{lam} values used.}
#' \item{df}{the number of non-zero \eqn{\beta}'s in estimated coefficients for \code{Z} at each value of \code{lam}}
#' \item{npass}{total iteration conducted in computation.}
#' \item{error}{error message for path at each each value of \code{lam}. If 0, no error occurs.}
#' \item{call}{the call that produced this object.}
#' \item{dim}{dimension of coefficient matrix}
#'
#'
#' @examples
#'
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c( beta, rep(0, times = p - length(beta)) )
#' Comp_data = comp_simulation(n = n, p = p, rho = 0.2, sigma = 0.5, gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' Comp_data$Zc
#' m1 <- compCL(y = Comp_data$y, Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'              intercept = Comp_data$intercept,
#'              pf = rep(1, times = p),
#'              lam = NULL, nlam = 100,lambda.factor = ifelse(n < p, 0.05, 0.001),
#'              dfmax = 20, #pfmax = min(dfmax * 1.5, p),
#'              mu_ratio = 1, tol = 1e-10,
#'              outer_maxiter = 1e8, outer_eps = 1e-10,
#'              inner_maxiter = 1e3, inner_eps = 1e-6)
#' drop(m1$lam)
#' coef(m1, s = c(1, 0.5, 0.1))
#' #coef(m1)
#' Znew <- Comp_data$Z[1:5,]
#' print(predict(m1, Znew = Znew[1:5, ], s = m1$lam[15:20]))
#'
#'
#' @export
#'
#'



compCL <- function(y, Z, Zc = NULL, intercept = TRUE,
                   pf = rep(1, times = p),
                   lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
                   dfmax = p, pfmax = min(dfmax * 1.5, p),
                   u = 1,
                   mu_ratio = 1.01, tol = 1e-10,
                   outer_maxiter = 1e+08, outer_eps = 1e-8,
                   inner_maxiter = 1e+3, inner_eps = 1e-6) {
  #u <- 1
  this.call <- match.call()
  y <- drop(y)
  Z <- as.matrix(Z)

  if( any(abs(rowSums(Z) - 1) > 1e-10) ) {
    message("Z is transformed into compositional data by deviding rowSums")
    Z <- Z / rowSums(Z)
  }
  if(any(Z == 0)) stop("There is zero entry in compositional data")
  Z <- log(Z)

  n <- length(y)
  p <- dim(Z)[2]
  inter <- as.integer(intercept)
  Znames <- colnames(Z)
  if (is.null(Znames)) Znames <- paste0("Z", seq(p))

  Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }
  beta_c <- as.vector(Zc_proj %*% y)
  beta.ini <- c(rep(0, times = p), beta_c)

  m <- dim(Zc)[2]
  if(m > 1) {
    Zcnames <- colnames(Z)
    if (is.null(Zcnames)) Zcnames <- paste0("Zc", seq(m-1))
  } else {
    Zcnames <- NULL
  }

  if(is.null(lam)) {
    lam0 <- t(Z) %*% (y - Zc %*% beta_c) / n
    lam0 <- max(abs(lam0) / pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  }

  pfmax <- as.integer(pfmax)
  dfmax <- as.integer(dfmax)
  inner_maxiter <- as.integer(inner_maxiter)
  outer_maxiter <- as.integer(outer_maxiter)
  inner_eps <- as.double(inner_eps)
  outer_eps <- as.double(outer_eps)
  tol <- as.double(tol)
  u <- as.double(u)
  mu_ratio <- as.double(mu_ratio)
  #estar <- as.double(estar)

  output <- ALM_B(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini,
                  lambda = lam, pf = pf, dfmax = dfmax, pfmax = pfmax,
                  inner_eps = inner_eps, outer_eps = outer_eps,
                  inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter,
                  u_ini = u, mu_ratio = mu_ratio, tol = tol

  )

  output$call <- this.call
  class(output) <- "compCL"
  output$dim <- dim(output$beta)
  output$lam <- drop(output$lam)
  output$Z_log <- Z

  dimnames(output$beta) <- list(c(Znames, Zcnames, "Intercept"), paste0("L", seq(along = output$lam)) )
  return(output)

}


