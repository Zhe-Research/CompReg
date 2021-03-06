% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.R
\name{cv.compCL}
\alias{cv.compCL}
\title{Cross-validation for compCL}
\usage{
cv.compCL(y, Z, Zc = NULL, intercept = FALSE, lam = NULL, nfolds = 10,
  foldid, trim = 0.1, ...)
}
\arguments{
\item{y}{a vector of response variable with length n.}

\item{Z}{a \eqn{n*p} matrix after taking log transformation on compositional data.}

\item{Zc}{a design matrix of other covariates considered. Default is \code{NULL}.}

\item{intercept}{Whether to include intercept in the model. Default is TRUE.}

\item{lam}{a user supplied lambda sequence. Typically, by leaving this option unspecified users can have the
program compute its own \code{lam} sequence based on \code{nlam} and \code{lambda.factor}
If \code{lam} is provided but a scaler, \code{lam} sequence is also created starting from \code{lam}.
Supplying a value of lambda overrides this. It is better to supply a decreasing sequence of lambda
values, if not, the program will sort user-defined \code{lambda} sequence in decreasing order
automatically.}

\item{nfolds}{number of folds - default is 10. Smallest value allowable is nfolds=3.}

\item{foldid}{an optional vector of values between 1 and \code{nfolds} identifying what fold each
observation is in. If supplied, \code{nfold} can be missing.}

\item{trim}{a scaler specifying percentage to be trimmed off for prediction error - default is 0.}

\item{...}{other arguments that can be passed to compCL.}
}
\value{
an object of class \code{\link{cv.compCL}} is returned.
\item{compCL.fit}{a fitted \code{\link{compCL}} object for the full data}
\item{lam}{the values of \code{lam} used in the fits}
\item{Ftrim}{a list of cross-validation result without trimming.
               \itemize{
               \item \code{cvm} the mean cross-validated error without trimming -  a vector of  \code{length(lam)}
               \item \code{cvsd} estimate of standard error of cvm without trimming- a vector of  \code{llength(lam)}
               \item \code{cvupper} upper curve = \code{cvm+cvsd}.
               \item \code{cvlo} lower curve = \code{cvm-cvsd}.
               \item \code{lam.min} The optimal value of \code{lam} that gives minimum cross validation error \code{cvm}
               \item \code{lam.1se} The largest value of lam such that error is within 1 standard error of the minimum \code{cvm}
               }
            }

\item{Ttrim}{a list of cross-validation result with \code{trim*100\%}, if provided, of tails trimmed off for cross validation error.
               \itemize{
               \item \code{cvm} the mean cross-validated error with with \code{trim*100\%} trimmed - a vector of  \code{length(lam)}
               \item \code{cvsd} estimate of standard error of cvm with \code{trim*100\%} trimmed - a vector of  \code{length(lam)}
               \item \code{cvupper} upper curve = \code{cvm+cvsd}.
               \item \code{cvlo} lower curve = \code{cvm-cvsd}.
               \item \code{lam.min} The optimal value of \code{lam} that gives minimum cross validation error cvm with \code{trim*100\%} trimmed
               \item \code{lam.1se} The largest value of lam such that error is within 1 standard error of the minimum \code{cvm} after \code{trim*100\%} trimmed.
               }
            }

\item{foldid}{the values of \code{folidi} used in fits.}
}
\description{
Does nfolds cross-validation for compCL, return value of \code{lam}.
The function is modified based on the \code{cv} function from \code{glmnet} package
}
\examples{

p = 30
n = 50
beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
beta = c( beta, rep(0, times = p - length(beta)) )
Comp_data = comp_simulation(n = n, p = p,
                            rho = 0.2, sigma = 0.5,
                            gamma  = 0.5, add.on = 1:5,
                            beta = beta, intercept = FALSE)
Comp_data$Zc
cvm <- cv.compCL(y = Comp_data$y,
                 Z = Comp_data$X.comp, Zc = Comp_data$Zc,
                 intercept = Comp_data$intercept,
                 lam = NULL, nfolds = 10, trim = 0.05, lambda.factor = 0.0001,
                 dfmax = p, mu_ratio = 1, outer_eps = 1e-10, inner_eps = 1e-8, inner_maxiter = 1e4)

plot(cvm)
coef(cvm, s = "lam.min")
cvm$compCL.fit
#apply(cvm$compCL.fit$beta[1:p, ], 2, function(x) which(abs(x) > 0))
which(abs(coef(cvm, s = "lam.min")$beta[1:p]) > 0)
which(abs(coef(cvm, s= "lam.1se")$beta[1:p]) > 0)

}
