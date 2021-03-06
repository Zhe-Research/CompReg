% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compCL.R
\name{compCL}
\alias{compCL}
\title{Fits regularization paths for longitudinal compositional data with lasso penalty.}
\usage{
compCL(y, Z, Zc = NULL, intercept = TRUE, pf = rep(1, times = p),
  lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
  dfmax = p, pfmax = min(dfmax * 1.5, p), u = 1, mu_ratio = 1.01,
  tol = 1e-10, outer_maxiter = 1e+08, outer_eps = 1e-08,
  inner_maxiter = 1000, inner_eps = 1e-06)
}
\arguments{
\item{y}{a vector of response variable with length n.}

\item{Z}{a \eqn{n*p} matrix after taking log transformation on compositional data.}

\item{Zc}{a design matrix of other covariates considered. Default is \code{NULL}.}

\item{intercept}{Whether to include intercept in the model. Default is TRUE.}

\item{pf}{penalty factor, a vector in length of p. Separate penalty weights can be applied to each coefficience
\eqn{\beta} for composition variates to allow differential shrinkage. Can be 0 for some \eqn{\beta}'s,
which implies no shrinkage, and results in that composition always being included in the model.
Default value for each entry is the 1.}

\item{lam}{a user supplied lambda sequence. Typically, by leaving this option unspecified users can have the
program compute its own \code{lam} sequence based on \code{nlam} and \code{lambda.factor}
If \code{lam} is provided but a scaler, \code{lam} sequence is also created starting from \code{lam}.
Supplying a value of lambda overrides this. It is better to supply a decreasing sequence of lambda
values, if not, the program will sort user-defined \code{lambda} sequence in decreasing order
automatically.}

\item{nlam}{the length of \code{lam} sequence - default is 100.}

\item{lambda.factor}{the factor for getting the minimal lambda in \code{lam} sequence, where
\code{min(lam)} = \code{lambda.factor} * \code{max(lam)}.
\code{max(lam)} is the smallest value of \code{lam} for which all penalized group are zero's.
The default depends on the relationship between \eqn{n}
and \eqn{p1}
If \eqn{n >= p1} the default is \code{0.001}, close to zero.
If \eqn{n < p1}, the default is \code{0.05}. A very small value of
\code{lambda.factor}
will lead to a saturated fit. It takes no effect if there is user-defined lambda sequence.}

\item{dfmax}{limit the maximum number of groups in the model. Useful for very large \eqn{p},
if a partial path is desired - default is \eqn{p}.}

\item{pfmax}{limit the maximum number of groups ever to be nonzero. For example once a group enters the
model along the path, no matter how many times it exits or re-enters model through the path,
it will be counted only once. Default is \code{min(dfmax*1.5, p)}.}

\item{u}{\code{u} is the inital value for penalty parameter of augmented Lanrange method.
Default value is 1.}

\item{mu_ratio}{\code{mu_ratio} is the increasing ratio for \code{u} - default value is 1.01.
Inital values for scaled Lagrange multipliers are set as 0's.
If \code{mu_ratio} < 1,
there is no linear constraints included. Group lasso coefficients are estimated.}

\item{tol}{tolerance for vectors beta'ss to be considered as none zero's. For example, coefficient
\eqn{\beta_j} for group j, if \eqn{max(abs(\beta_j))} < \code{tol}, set \eqn{\beta_j} as 0's.
Default value is 0.}

\item{outer_maxiter}{\code{outer_maxiter} is the maximun munber of loops allowed for Augmented Lanrange method;
and \code{outer_eps} is the convergence termination tolerance.}

\item{outer_eps}{\code{outer_maxiter} is the maximun munber of loops allowed for Augmented Lanrange method;
and \code{outer_eps} is the convergence termination tolerance.}

\item{inner_maxiter}{\code{inner_maxiter} is the maximun munber of loops allowed for blockwise-GMD;
and \code{inner_eps} is the convergence termination tolerance.}

\item{inner_eps}{\code{inner_maxiter} is the maximun munber of loops allowed for blockwise-GMD;
and \code{inner_eps} is the convergence termination tolerance.}
}
\value{
An object with S3 calss \code{\link{compCL}}
\item{beta}{a matrix of coefficients for \code{cbind{Z, Zc, 1_n}}, with \code{nlam} rows.}
\item{lam}{the actual sequence of \code{lam} values used.}
\item{df}{the number of non-zero \eqn{\beta}'s in estimated coefficients for \code{Z} at each value of \code{lam}}
\item{npass}{total iteration conducted in computation.}
\item{error}{error message for path at each each value of \code{lam}. If 0, no error occurs.}
\item{call}{the call that produced this object.}
\item{dim}{dimension of coefficient matrix}
}
\description{
Fits regularization paths for longitudinal compositional data with lasso penalty at a sequence of regularization parameters lambda.
}
\examples{

p = 30
n = 50
beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
beta = c( beta, rep(0, times = p - length(beta)) )
Comp_data = comp_simulation(n = n, p = p, rho = 0.2, sigma = 0.5, gamma  = 0.5, add.on = 1:5,
                            beta = beta, intercept = FALSE)
Comp_data$Zc
m1 <- compCL(y = Comp_data$y, Z = Comp_data$X.comp, Zc = Comp_data$Zc,
             intercept = Comp_data$intercept,
             pf = rep(1, times = p),
             lam = NULL, nlam = 100,lambda.factor = ifelse(n < p, 0.05, 0.001),
             dfmax = 20, #pfmax = min(dfmax * 1.5, p),
             mu_ratio = 1, tol = 1e-10,
             outer_maxiter = 1e8, outer_eps = 1e-10,
             inner_maxiter = 1e3, inner_eps = 1e-6)
drop(m1$lam)
coef(m1, s = c(1, 0.5, 0.1))
#coef(m1)
Znew <- Comp_data$Z[1:5,]
print(predict(m1, Znew = Znew[1:5, ], s = m1$lam[15:20]))


}
