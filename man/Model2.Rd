% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{Model2}
\alias{Model2}
\title{Simulation model2 for longitudinal composition data}
\usage{
Model2(n, p, m = 0, intercept = TRUE, interval = c(0, 1), ns = 100,
  obs_spar = 0.6, discrete = FALSE, SNR = 1, sigma = 2, rho_X,
  Corr_X = c("CorrAR", "CorrCS"), rho_W, Corr_W = c("CorrAR", "CorrCS"),
  Nzero_group = 4, range_beta = c(0.5, 1), range_beta_c = 1, beta_C,
  theta.add = c(1, 2, 5, 6), gamma = 0.5, basis_W = c("bs", "OBasis",
  "fourier"), df_W = 5, degree_W = 3, basis_beta = c("bs", "OBasis",
  "fourier"), df_beta = 5, degree_beta = 3, insert = c("FALSE", "basis"),
  method = c("trapezoidal", "step"))
}
\arguments{
\item{n}{sample size}

\item{p}{size of compositional predictors falling in \eqn{S^p}}

\item{m}{size of time-invariant predictors. First \code{ceiling(m/2)} columns are generated by  \code{bin(1,0.5)} independently;
latter \code{(m - ceiling(m/2))} columns are generated \code{norm(0, 1)} independently.}

\item{intercept}{including intercept or not to generate response variable, default is TRUE}

\item{interval}{a length 2 vector indicating time domain.}

\item{ns}{\code{ns} is a scaler specifying length of equally spaced time sequence on domian \code{interval}.}

\item{obs_spar}{a percentage used to get sparse ovbservation. Each time point probability \code{obs_spar} to be observed. It allows different subject with
different observed time points and size.
\code{obs_spar * ns > 5} is required.}

\item{discrete}{is logical, specifying whether \eqn{X} is generated at different time points.
If \code{distrete = TRUE}, generate \eqn{X} on dense sequence created by \code{max(ns_dense = 1000, 5*ns)} and Then for
each subject, random sample \code{ns} points. recommend \code{ns < 200} when \code{distrete = TRUE}.}

\item{SNR}{signal to noise ratio.}

\item{sigma, rho_X, Corr_X, rho_W, Corr_W}{linear combination scaler \code{W}, \code{df_W*p}, is generated from from Multivariate Normal Distribution with mean 0's and
ovariance matrix = \code{simga^2 * kronecker(Sigma_X, Sigma_W)}.
\code{Corr_X} is correlation structure of \code{Sigma_X} with \eqn{\rho} = \code{rho_X},
which controls canonical-correlation between groups for W;
\code{Corr_W} is correlation structure of \code{Sigma_W} with \eqn{\rho} = \code{rho_W},
which controls correlation within groups of W.}

\item{Nzero_group}{a even scaler. First \code{Nzero_group} compositional predictors are considered having none zero effect, while others are with 0 coefficients.}

\item{range_beta}{a sorted vector of length 2 used to generate coefficient matrix \code{B} for compositional predict, which is with demension \code{p*k}.
For each column of \code{B}, generate \code{Nzero_group/2} from unifom distribution with range \code{range_beta}, and together with their negatives
are ramdom assigned to the first \code{Nzero_group} rows.}

\item{range_beta_c}{value of coefficients for beta0 and beta_c (coefficients for time-invariant predictors)}

\item{beta_C}{vectorized coefficients of coefficient matrix for compositional predictors. Could be missing.}

\item{theta.add}{logical or numerical. If numerical, indicating which ones of compositional predictors of high level mean.
If logical, \code{c(1:ceiling(Nzero_group/2), Nzero_group + (1:ceiling(Nzero_group/2)))} are set to with high level mean}

\item{gamma}{high level mean groups adding log(p * gamma) before convertint into compositional data, otherwise 0.}

\item{basis_W, df_W, degree_W}{longitudinal compositional data is generated from linear combination of basis \eqn{\Psi(t)}, take exponetial and change into compositional data.
\itemize{
   \item \code{basis_W} is the basis function for \eqn{\Psi(t)}  - default is \code{"bs"}.
          Other choise are \code{"OBasis"} and \code{"fourier"};
   \item \code{df_W} is the degree of freedom for basis \eqn{\Psi(t)} - default is 10 ;
   \item \code{degree_W} the is degree for \eqn{\Psi(t)} - default is 3.
 }}

\item{basis_beta, df_beta, degree_beta}{coefficinet curve is generate by linear combination of basis \eqn{\Phi(t)}.
\itemize{
     \item \code{basis_beta} is the basis function for \eqn{\Phi(t)} - default is \code{"bs"}.
            Other choise are \code{"OBasis"} and \code{"fourier"};
      \item \code{df_beta} is the degree of freedom for basis \eqn{\Phi(t)} - default is 5;
      \item \code{degree_beta} is the degree for \eqn{\Phi(t)} - default is 3.
 }}

\item{insert}{way to interpolation.
\itemize{
   \item \code{"FALSE"} no interpolation.
   \item \code{"basis"} compositional data is considered as step function, imposing basis on un-observed time points for each subject.
}
Default is \code{"FALSE"}.}

\item{method}{method used to approximate integral.
\itemize{
\item \code{"trapezoidal"} Sum up area under trapezoidal formulated by values of function at two adjacent observed time points. See \code{\link{ITG_trap}}.
\item \code{"step"} Sum up area under rectangle formulated by step function at observed time points. See \code{\link{ITG_step}}.
}
Default is \code{"trapezoidal"}}
}
\value{
a list
\item{data}{a list, a vector \code{y} of response variable, a data frame \code{Comp} of sparse observation of longitudinal compositional data,
                    a matrix \code{Zc} of time-invariable predictors, a logtical \code{intercept}.}
\item{beta}{a length \code{p*df_beta + m + 1} vector for coefficients}
\item{basis.info}{ matrix for basis for beta, combining the first column as time sequence.}
\item{data.raw}{ a list, \code{Z_t.full} full observation of longitudinal compositional data,
                         \code{Z_ITG} integral for full observation of longitudinal compositional data,
                         \code{Y.tru} true response without noise,
                         \code{X} longitudinal before converting into compositional data,
                         \code{W} matrix of linear combination scalers, \code{n * (df_W * p)}.}

\item{parameter}{ a list of parameters.}
}
\description{
Simulate sparse observation from longitudinal compositional data \eqn{X}.
}
\examples{

df_beta = 5
p = 20
Data <- Model2(n = 50, p = p, m = 2, intercept = TRUE, ns = 50, SNR = 1,
              rho_X = 0.1, rho_W = 0.2, df_W = 10, df_beta = df_beta, obs_spar = 0.5)
names(Data$data)
Data.test <- Model2(n = 50, p = p, m = 2, intercept = TRUE, ns = 50, SNR = 1,
                   rho_X = 0.1, rho_W = 0.2, df_W = 10, df_beta = df_beta, obs_spar = 0.5,
                   beta_C = Data$beta[1:(p*df_beta)])


}
