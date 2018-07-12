#'
#' GIC cirterion selection for compCL
#'
#' Calculate GIC for compCL, return value of \code{lam}.
#' The function follows Variable selection in regression with compositional covariates by
#' WEI LIN, PIXU SHI, RUI FENG AND HONGZHE LI
#'
#' @inheritParams compCL
#'
#' @param \dots other arguments that can be passed to compCL.
#'
#' @return an object of class \code{\link{GIC.compCL}} is returned.
#' \item{compCL.fit}{a fitted \code{\link{compCL}} object for the full data}
#' \item{lam}{the values of \code{lam} used in the fits}
#' \item{GIC}{a vector of GIC values for each \code{lam}}
#' \item{lam.min}{\code{lam} value such that minimize \eqn{GIC(\lambda)} }
#'
#' @details
#' \deqn{\textrm{GIC}(\lambda) = \log{\hat{\sigma}^2_\lambda}
#'                       + (s_\lambda - 1) \frac{\log{\log{n}}}{n} \log{max(p, n)} },
#' where \eqn{\hat{\sigma}^2_\lambda} is the MSE for fitted path.
#'
#'
# \deqn{F_{n} =\frac{\phi^{n} - \psi^{n}}{\sqrt{5}}.}
#
# deqn ASCII example
#'
#\deqn{ \sigma = \sqrt{ \frac{Z}{n} \sum
#  \left[ \textstyle\frac{1}{2}\displaystyle
#    \left( \log \frac{H_i}{L_i} \right)^2  - (2\log 2-1)
#    \left( \log \frac{C_i}{O_i} \right)^2 \right] }
#}{sqrt(N/n * runSum(0.5 * log(OHLC[,2]/OHLC[,3])^2 -
#           (2*log(2)-1) * log(OHLC[,4]/OHLC[,1])^2, n))}
#'
#' @examples
#'
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_simulation(n = n, p = p,
#'                             rho = 0.2, sigma = 0.5,
#'                             gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' Comp_data$Zc
#' GICm <- GIC.compCL(y = Comp_data$y,
#'                    Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                    intercept = Comp_data$intercept,
#'                    lam = NULL,lambda.factor = 0.0001,
#'                    dfmax = p, outer_eps = 1e-10, mu_ratio = 1)
#' coef(GICm)
#' plot(y = GICm$GIC, x = log(GICm$lam), ylab = "GIC", xlab = "Log(Lambda)" )
#' @export
#'


GIC.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
                       lam = NULL, ...) {
  digits = 5
  this.call <- match.call()
  y <- drop(y)
  n <- length(y)

  compCL.object <- compCL(y = y, Z = Z, Zc = Zc, intercept = intercept, lam = lam, ...)
  #print(compCL.object)
  #cat(class(compCL.object), "\r\n")
  lam <- compCL.object$lam
  GIC <- GIC.test2(object = compCL.object, y = y, Z =  compCL.object$Z_log, Zc = Zc, intercept = intercept)
  #print(GIC.test2(object = compCL.object, y = y, Z = Z, Zc = Zc, intercept = intercept))
  GIC1 <- round(GIC, digits = digits)
  #GIC.min <- min(GIC1)
  GIC.min <- min(GIC1[drop(compCL.object$df) > 0])

  idmin <- GIC1 <= GIC.min
  idmin[drop(compCL.object$df) < 2 ] <- FALSE
  lam.min <- max(lam[idmin])

  result <- list(compCL.fit = compCL.object,
                 lam = lam,
                 GIC = GIC,
                 lam.min = lam.min)

  class(result) <- "GIC.compCL"
  result$call <- this.call
  return(result)
}



#' @title
#' GIC cirterion selection for FuncompCGL
#'
#' @description
#' Calculate GIC for compCL, return value of \code{lam}.
#'
#### @usage
#### GIC.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL,
####                            ref = NULL,
####                            W = rep(1,times = p - length(ref)),
####                            k = 4:10, nlam = 100, outer_maxiter = 1e+6,
####                            cut_off = c("Curve","Matrix", "NULL"),
####                            lower_tri = 0.01, ...)
####
#'
#'
#' @inheritParams cv.FuncompCGL
#'
#' @param \dots other arguments that could be passed to FuncompCL.
#'
#' @return an object of class \code{\link{GIC.FuncompCGL}} is returned.
#' \item{Funcomp.CGL.fit}{a list, length of \code{k},
#'                        of fitted \code{\link{FuncompCGL}} object for the full data.
#'                        objects with S3 calss \code{\link{FuncompCGL}}}
#' \item{lam}{the values of \code{lam} used in the fits}
#' \item{MSE}{matrix of mean squared error with size \code{k} by \code{nlam} (the length of
#'            actually used \code{lambda} sequence, migth pre-stop by \code{dfmax} or
#'            \code{pfmax}). MSE is equivalent to likelihood under normal error model. \cr
#'            \strong{Could be edited for other linkage.}}
#' \item{Nzero}{a \code{k} by nlam matrix for Nzero group cut-off by \code{cut_off} and \code{lower_tri}}
#'
#' @examples
#'
#' df_beta = 5
#' p = 30 #30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[3, ] <- c(-1, 0, 0, 0, -0.5) #c(-0.5, 0, 0, 0, -0.5)
#' beta_C_true[1, ] <- c(1, 0, 1 , 0, -0.5) #c(0.5, 0, 1 , 0, -0.5)
#' beta_C_true[2, ] <- c(0, 0,  -1,  0,  1)
#'
#' nfolds = 10
#' k_list <- c(4,5,6)
#' n_train = 100 #100
#' n_test = 500
#'
#' Data <- Model2(n = n_train, p = p, m = 0, intercept = TRUE,
#'                SNR = 3, sigma = 3,
#'                rho_X = 0, rho_W = 0,
#'                Corr_X = "CorrCS", Corr_W = "CorrAR",
#'                df_W = 5, df_beta = df_beta,
#'                ns = 20, obs_spar = 1, theta.add = FALSE, #c(0,0,0),
#'                beta_C = as.vector(t(beta_C_true)))
#' y <- drop(Data$data$y)
#' n <- length(y)
#' X <- Data$data$Comp
#' Zc <- Data$data$Zc
#' intercept <- Data$data$intercept
#' m <- ifelse(is.null(Zc), 0, dim(Zc)[2]) #+ as.integer(intercept)
#' m1 <- m + as.integer(intercept)
#' sseq <- Data$basis.info[,1]
#' beta_C.true <- matrix(Data$beta[1:(p*(df_beta))],
#'                       nrow = p, ncol = df_beta, byrow = TRUE)
#' beta_curve.true <- Data$basis.info[,-1] %*% t(beta_C.true)
#' Non_zero.true <- (1:p)[apply(beta_C.true, 1, function(x) max(abs(x)) > 0)]
#' foldid <- sample(rep(seq(nfolds), length = n))
#'
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Model2, arg_list)
#' y_test <- drop(Test$data$y)
#'
#' GIC_arg <- list(basis_fun = "bs", degree = 3, sseq = Data$basis.info[, 1])
#'
#' # y = y
#' # X = X
#' # Zc = Zc
#' # intercept = intercept
#' # W = rep(1, p)
#' # k = k_list
#' # nfolds = 10
#' # trim = 0
#' # tol = 0
#' # inner_eps = 1e-6
#' # inner_maxiter = 1E3
#' # dfmax = 20
#' # lambda.factor = 1e-20
#' # mu_ratio = 1
#' # outer_eps = 1e-6
#' # keep = TRUE
#' # Trange = c(0,1)
#'
#'
#'
#' GIC_m1 <- GIC.FuncompCGL( y = y, X = X, Zc = Zc, ref = NULL,
#'                           inner_eps = 1e-8, outer_eps = 1e-8, tol = 1e-8,
#'                           k = k_list)
#'
#' temp <- get.GIC(p = p, df_list = k_list, lower_tri = 0.01,
#'                 GIC_obj = GIC_m1, GIC_arg = GIC_arg,
#'                 cut_type = "Strict", GIC_type = "GIC1",
#'                 method_type = "cgl", refit = FALSE)
#' GIC_curve <- temp$GIC_curve
#' k_opt <- temp$k_opt
#' beta_GIC <- temp$beta
#'
#' plot.args = list(x = seq(length(GIC_m1$lam)), #GIC_m1$lam, #log(GIC_m1$lam),
#'                  y = GIC_curve[1, ],
#'                  ylim = range(GIC_curve),
#'                  xlab= "lambda Index",#"lambda", #"log(lambda)",
#'                  ylab="GIC",
#'                  type="n")
#' #do.call("plot",plot.args)
#' # for(i in 1:length(k_list)) {
#' #
#' #   points(x = seq(length(GIC_m1$lam)), #GIC_m1$lam, #log(GIC_m1$lam),
#' #          y = GIC_curve[i, ], col = rainbow(length(k_list))[i])
#' #   text(length(GIC_m1$lam), #tail(log(GIC_m1$lam), 1),
#' #        GIC_curve[i, length(GIC_m1$lam)], labels=paste(k_list[i]),
#' #        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' # }
#' # axis(3, at = pretty(seq(length(GIC_m1$lam))), labels = rev(pretty(GIC_m1$lam)))
#' # loc  = which(GIC_curve == min(GIC_curve), arr.ind = TRUE)
#'
#'
#'
#'
#' beta_C <- matrix(beta_GIC[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C) , "\r\n")
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#' do.call("plot",plot.args)
#' for(i in 1:length(k_list)) {
#'   points(x = seq(length(GIC_m1$lam)), #log(GIC_m1$lam),
#'          y = GIC_curve[i, ], col = rainbow(length(k_list))[i], pch = seq(length(k_list))[i])
#'   text(length(GIC_m1$lam), #tail(log(GIC_m1$lam), 1),
#'        GIC_curve[i, length(GIC_m1$lam)], labels=paste(k_list[i]),
#'        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' }
#' #axis(3, at = pretty(seq(length(GIC_m1$lam))), labels = rev(pretty(GIC_m1$lam)))
#'
#' matplot(sseq, beta_curve.true,
#'         ylab = "coeffcients curve", xlab = "TIME", #main = "TRUE",
#'         ylim = range(Data$beta[1:(p*df_beta)]),
#'         type = "l")
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("TRUE", line = 0.5)
#' text(0, beta_curve.true[1, Non_zero.true], labels = paste(Non_zero.true))
#'
#' B <- splines::bs(Data$basis.info[,1], df = k_opt, intercept = TRUE)
#' beta_curve <- B %*% t(beta_C)
#' matplot(sseq, beta_curve,
#'         ylab = "coef", xlab = "TIME", #main = "ESTI",
#'         ylim = range(Data$beta[1:(p*df_beta)])#,
#'         #type = "l"
#' )
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("Estimate", line = 0.5)
#' text(0, beta_curve[1, Non.zero], labels = paste(Non.zero))
#' text(tail(sseq, 1), beta_curve[dim(beta_curve)[1], Non.zero], labels = paste(Non.zero))
#' plot(apply(abs(beta_C),1,sum))
#' text(seq(length(GIC_m1$lam))[which(apply(abs(beta_C),1,sum) > 0)], #tail(log(GIC_m1$lam), 1),
#'      apply(abs(beta_C),1,sum)[which(apply(abs(beta_C),1,sum) > 0)],
#'      labels=paste(seq(length(GIC_m1$lam))[which(apply(abs(beta_C),1,sum) > 0)]),
#'      cex= 1, pos= 4)
#'
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method cgl"), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#'
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'   ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#'
#' cgl_GIC <- list()
#' MSE <- temp$MSE[temp$k_opt - k_list[1] + 1, temp$lam_loc]
#' R_sqr <- 1 - MSE * length(y) / crossprod(y -  mean(y))
#'
#' obj <- FuncompCGL(y = y_test, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
#' X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
#' PE <- sum((y_test - X_test %*% beta_GIC)^2) / length(y_test)
#' cgl_GIC$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#'
#' cgl_GIC$Non.zero_cut <- Non.zero
#' cgl_GIC <- c(cgl_GIC,
#'              ERROR_fun(beta_fit = beta_GIC, beta_true = Data$beta,
#'                        basis_fit = B, basis_true = Data$basis.info[,-1],
#'                        sseq = Data$basis.info[, 1],
#'                        m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#'              k = k_opt)
#' cgl_GIC$coef <- list(beta_C = beta_C, beta_c = tail(beta_GIC, m1))
#'
#' \dontrun{
#'
#' GIC_m2 <- GIC.FuncompCGL( y = y, X = X, Zc = Zc, ref = NULL,
#'                           outer_eps = 1e-8, mu_ratio = 0, tol = 1e-8,
#'                           k = k_list)
#' temp <- get.GIC(p = p, df_list = k_list, lower_tri = 0.01,
#'                 GIC_obj = GIC_m1, GIC_arg = GIC_arg,
#'                 cut_type = "Strict", GIC_type = "GIC1",
#'                 method_type = "naive", refit = FALSE)
#' GIC_curve <- temp$GIC_curve
#' k_opt <- temp$k_opt
#' beta_GIC <- temp$beta
#'
#' plot.args = list(x = seq(length(GIC_m2$lam)), #GIC_m2$lam, #log(GIC_m2$lam),
#'                  y = GIC_curve[1, ],
#'                  ylim = range(GIC_curve),
#'                  xlab= "lambda Index",#"lambda", #"log(lambda)",
#'                  ylab="GIC",
#'                  type="n")
#' # do.call("plot",plot.args)
#' #
#' # for(i in 1:length(k_list)) {
#' #
#' #   points(x = seq(length(GIC_m2$lam)), #GIC_m2$lam, #log(GIC_m2$lam),
#' #          y = GIC_curve[i, ], col = rainbow(length(k_list))[i])
#' #   text(length(GIC_m2$lam), #tail(log(GIC_m2$lam), 1),
#' #        GIC_curve[i, length(GIC_m2$lam)], labels=paste(k_list[i]),
#' #        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' # }
#' # axis(3, at = pretty(seq(length(GIC_m2$lam))), labels = rev(pretty(GIC_m2$lam)))
#'
#' beta_C <- matrix(beta_GIC[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C), "\r\n")
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#'
#' do.call("plot",plot.args)
#' for(i in 1:length(k_list)) {
#'   points(x = seq(length(GIC_m2$lam)), #GIC_m2$lam, #log(GIC_m2$lam),
#'          y = GIC_curve[i, ], col = rainbow(length(k_list))[i], pch = seq(length(k_list))[i])
#'   text(length(GIC_m2$lam), #tail(log(GIC_m2$lam), 1),
#'        GIC_curve[i, length(GIC_m2$lam)], labels=paste(k_list[i]),
#'        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' }
#' #axis(3, at = pretty(seq(length(GIC_m2$lam))), labels = rev(pretty(GIC_m2$lam)))
#'
#' matplot(sseq, beta_curve.true,
#'         ylab = "coeffcients curve", xlab = "TIME", #main = "TRUE",
#'         ylim = range(Data$beta[1:(p*df_beta)]),
#'         type = "l")
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("TRUE", line = 0.5)
#' text(0, beta_curve.true[1, Non_zero.true], labels = paste(Non_zero.true))
#'
#' B <- splines::bs(Data$basis.info[,1], df = k_opt, intercept = TRUE)
#' beta_curve <- B %*% t(beta_C)
#' matplot(sseq, beta_curve,
#'         ylab = "coef", xlab = "TIME", #main = "ESTI",
#'         ylim = range(Data$beta[1:(p*df_beta)])#,
#'         #type = "l"
#' )
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("Estimate", line = 0.5)
#' text(0, beta_curve[1, Non.zero], labels = paste(Non.zero))
#' text(tail(sseq, 1), beta_curve[dim(beta_curve)[1], Non.zero], labels = paste(Non.zero))
#'
#'
#' plot(apply(abs(beta_C),1,sum))
#' text(seq(length(GIC_m2$lam))[which(apply(abs(beta_C),1,sum) > 0)], #tail(log(GIC_m2$lam), 1),
#'      apply(abs(beta_C),1,sum)[which(apply(abs(beta_C),1,sum) > 0)],
#'      labels=paste(seq(length(GIC_m2$lam))[which(apply(abs(beta_C),1,sum) > 0)]),
#'      cex= 1, pos= 4)
#'
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method naive"), outer=TRUE, line = -2)
#'
#' par(mfrow=c(1,1))
#'
#'
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'   ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#'
#' naive_GIC <- list()
#' MSE <- temp$MSE[temp$k_opt - k_list[1] + 1, temp$lam_loc]
#' R_sqr <- 1 - MSE * length(y) / crossprod(y -  mean(y))
#' obj <- FuncompCGL(y = y_test, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
#' X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
#' PE <- sum((y_test - X_test %*% beta_GIC)^2) / length(y_test)
#' naive_GIC$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#' naive_GIC$Non.zero_cut <- Non.zero
#' naive_GIC <- c(naive_GIC,
#'                ERROR_fun(beta_fit = beta_GIC, beta_true = Data$beta,
#'                          basis_fit = B, basis_true = Data$basis.info[,-1],
#'                          sseq = Data$basis.info[, 1],
#'                          m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#'                k = k_opt)
#' naive_GIC$coef <- list(beta_C = beta_C, beta_c = tail(beta_GIC, m1))
#'
#'
#'
#' GIC_m3 <- GIC.FuncompCGL(y = y, X = X, Zc = Zc, ref = sample(4:p, 1), #sample(1:p, 1),
#'                          outer_eps = 1e-8, mu_ratio = 0, tol = 1e-8,
#'                          k = k_list)
#' temp <- get.GIC(p = p, df_list = k_list, lower_tri = 0.01,
#'                 GIC_obj = GIC_m1, GIC_arg = GIC_arg,
#'                 cut_type = "Strict", GIC_type = "GIC1",
#'                 method_type = "base", refit = FALSE)
#' GIC_curve <- temp$GIC_curve
#' k_opt <- temp$k_opt
#' beta_GIC <- temp$beta
#' ref = GIC_m3$Funcomp.CGL.fit[[1]]$ref
#' beta_GIC <-  c(beta_GIC1[ifelse(ref==1, 0, 1):((ref-1)*k_opt)],
#'                -colSums(matrix(beta_GIC1[1:((p-1)*k_opt)], byrow = TRUE, ncol = k_opt)),
#'                beta_GIC1[((ref-1)*k_opt+1):(length(beta_GIC1))])
#'
#' plot.args = list(x = seq(length(GIC_m3$lam)), #GIC_m3$lam, #log(GIC_m3$lam),
#'                  y = GIC_curve[1, ],
#'                  ylim = range(GIC_curve),
#'                  xlab= "lambda index",#"lambda", #"log(lambda)",
#'                  ylab="GIC",
#'                  type="n")
#' # do.call("plot",plot.args)
#' #
#' # for(i in 1:length(k_list)) {
#' #
#' #   points(x = seq(length(GIC_m3$lam)), #GIC_m3$lam, #log(GIC_m3$lam),
#' #          y = GIC_curve[i, ], col = rainbow(length(k_list))[i])
#' #   text(length(GIC_m3$lam), #tail(log(GIC_m3$lam), 1),
#' #        GIC_curve[i, length(GIC_m3$lam)], labels=paste(k_list[i]),
#' #        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' #
#' # }
#' # axis(3, at = seq(1, length(GIC_m3$lam), length.out = 10),
#' #      labels = rev(GIC_m3$lam[seq(1, length(GIC_m3$lam), length.out = 10)]) )
#' # axis(3, at = pretty(seq(length(GIC_m3$lam))),
#' #         labels = rev(pretty(GIC_m3$lam, n = length(pretty(seq(length(GIC_m3$lam))))))
#' #      )
#'
#'
#'
#' beta_C <- matrix(beta_GIC[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C), "\r\n")
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#' do.call("plot",plot.args)
#' for(i in 1:length(k_list)) {
#'   points(x = seq(length(GIC_m3$lam)), #log(GIC_m3$lam),
#'          y = GIC_curve[i, ], col = rainbow(length(k_list))[i], pch = seq(length(k_list))[i])
#'   text(length(GIC_m3$lam), #tail(log(GIC_m3$lam), 1),
#'        GIC_curve[i, length(GIC_m3$lam)], labels=paste(k_list[i]),
#'        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' }
#' #axis(3, at = pretty(seq(length(GIC_m3$lam))), labels = rev(pretty(GIC_m3$lam)))
#' matplot(sseq, beta_curve.true,
#'         ylab = "coeffcients curve", xlab = "TIME", #main = "TRUE",
#'         ylim = range(Data$beta[1:(p*df_beta)]),
#'         type = "l")
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("TRUE", line = 0.5)
#' text(0, beta_curve.true[1, Non_zero.true], labels = paste(Non_zero.true))
#'
#' B <- splines::bs(Data$basis.info[,1], df = k_opt, intercept = TRUE)
#' beta_curve <- B %*% t(beta_C)
#' matplot(sseq, beta_curve,
#'         ylab = "coef", xlab = "TIME", #main = "ESTI",
#'         ylim = range(Data$beta[1:(p*df_beta)])
#'         #type = "l"
#' )
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("Estimate", line = 0.5)
#' text(0, beta_curve[1, Non.zero], labels = paste(Non.zero))
#' text(tail(sseq, 1), beta_curve[dim(beta_curve)[1], Non.zero], labels = paste(Non.zero))
#' plot(apply(abs(beta_C),1,sum))
#' text(seq(length(GIC_m3$lam))[which(apply(abs(beta_C),1,sum) > 0)], #tail(log(GIC_m3$lam), 1),
#'      apply(abs(beta_C),1,sum)[which(apply(abs(beta_C),1,sum) > 0)],
#'      labels=paste(seq(length(GIC_m3$lam))[which(apply(abs(beta_C),1,sum) > 0)]),
#'      cex= 1, pos= 4)
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method base",  ", ref=", ref), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'   ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#'
#'
#' base_GIC <- list()
#' MSE <- temp$MSE[temp$k_opt - k_list[1] + 1, temp$lam_loc]
#' R_sqr <- 1 - MSE * length(y) / crossprod(y -  mean(y))
#' obj <- FuncompCGL(y = y_test, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
#' X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
#' PE <- sum((y_test - X_test %*% beta_GIC)^2) / length(y_test)
#' base_GIC$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#' base_GIC$Non.zero_cut <- Non.zero
#' base_GIC <- c(base_GIC,
#'               ERROR_fun(beta_fit = beta_GIC, beta_true = Data$beta,
#'                         basis_fit = B, basis_true = Data$basis.info[,-1],
#'                         sseq = Data$basis.info[, 1],
#'                         m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#'               k = k_opt)
#' base_GIC$coef <- list(beta_C = beta_C, beta_c = tail(beta_GIC, m1))
#'
#' }
#'
#' @export
#'
#'


GIC.FuncompCGL <- function(y, X, Zc = NULL, ref = NULL,
                            lam = NULL, nlam = 100,
                            W = rep(1,times = p - length(ref)),
                            k = 4:10, outer_maxiter = 1e+6,...) {
  y <- drop(y)
  n <- length(y)
  object <- as.list(seq(length(k)))
  this.call <- match.call()
  names(object) <- k
  if(!is.null(lam) || length(k) == 1) {

    # Case I
    if(dim(X)[1] == n ) p = dim(X)[2] / k else p <- dim(X)[2] - 2

    for(i in 1:length(k)){
      ###cat("1", k[i], "\r\n")
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam,
                                W = W, ref = ref,
                                k = k[i],
                                nlam = nlam, outer_maxiter = outer_maxiter,
                                ...)
    }
  } else {

    # Case II: find commom lambda first
    p <- ncol(X) - 2
    for(i in 1:length(k)){

      ###cat("B", k[i], "\n")
      # Caculate integral Z and W matrix (if W is a functoin)
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1,
                                W = W, ref = ref,
                                k = k[i], outer_maxiter = 0, ...)
    }

    lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k

    # Solution path for each df k
    for(i in 1:length(k)) {
      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                k = k[i], ...)
      object[[i]]$call <- as.list(object[[i]]$call)
      object[[i]]$call$k <- k[i]
    }
  }

  GIC.nlam <- sapply(object, function(x) length(drop(x$lam)))
  GIC.nlam_id <- which.min(GIC.nlam)
  GIC.lam <- object[[GIC.nlam_id]]$lam
  GIC.nlam <- GIC.nlam[GIC.nlam_id]
  MSE <-matrix(NA, nrow = length(k), ncol = GIC.nlam)
  for(i in seq(length(k))) {
    predmat <- predict(object = object[[i]], newx = cbind(object[[i]]$Z, Zc))
    MSE[i, ] <- apply(predmat[, 1:GIC.nlam] - y, 2, function(x) mean(x^2))
  }

  result <- list(Funcomp.CGL.fit = object, lam = GIC.lam, MSE = MSE)
  result$class <- "GIC.FuncompCGL"
  result$call <- this.call
  return(result)
}



#' @title
#' GIC criterion selection for FuncompCGL
#'
#' @description
#' Calculate GIC for GIC.FuncompCGL, select beta coefficients vector and df via
#' BIC, AIC, GIC criterion
#'
#' @param p number of compositional predictors
#' @param df_list vector of df's used in \code{GIC_obj}
#' @param GIC_obj an object of \code{\link{GIC.FuncompCGL}}
#' @param GIC_arg argument list used to fit \code{\link{GIC.FuncompCGL}},
#'                need to at least provide information for basis generation,
#'                including \code{basis_fun}, \code{degree} and \code{sseq}.
#' @param cut_type cut method for none zero group - default value is \code{Strict} which
#'                 is no cut-off. Other two methods are \code{Matrix} and \code{Curve}.
#'                 See details.
#' @param lower_tri lower percentage boundary used in cut none zero groups.
#' @param GIC_type variations of GIC's, AIC and BIC. See details.
#' @param refit logical, whehter refit likelihood when small magnitude groups are cut to zeros'
#'              Default value is \code{FALSE}.
#' @param y,Zc when \code{refit=TRUE}, need to privde \code{y} and \code{Zc} to re-calcualte
#'             likelihood.
#' @param method_type mehod used to fitted \code{\link{GIC.FuncompCGL}} \cr
#'                    \code{cgl} - constrained group lasso \cr
#'                    \code{base} - log contract group lasso \cr
#'                    \code{naive} - do not consider the fact of compositional.
#'
#' @return
#' \item{beta}{selected beta vector.}
#' \item{k_opt}{selected \code{df}.}
#' \item{GIC_curve}{a \code{length(df_list)} by \code{nlam} matrix of GIC values}
#' \item{N_zero}{a matrix of numbers of none zero groups selected.}
# list(beta = beta_GIC, k_opt = k_opt, lam_loc = lam_loc,
#      GIC_min = GIC_min, GIC_curve = GIC_curve,
#      N_zero = N_zero, MSE = MSE, scaler = scaler,
#      method_type = method_type, GIC_type = GIC_type, cut_type = cut_type)
#'
#'
#' @details
#' $$GIC(\eqn{\lambda}) = log(MSE) + Selection * \code{alpha}$$, for normal error.
#' If include linear constraints, like in log contract and constraints group lasso model,
#' selection = None zero group - 1; otherwise in naive method without consideration of
#' linear constraints, selection = None zero group. \cr
#'
#' alpha_GIC1 = log(max(p*df, n)) * df / n \cr
#' alpha_GIC2 = log(max(p*df, n)) * df * log(log(n)) / n \cr
#' alpha_GIC3 = log(p*df) *df /n \cr
#' alpha_GIC4 = log(p*df) * df * log(log(n)) / n \cr
#' alpha_BIC = log(n) * df / n \cr
#' alpha_AIC = 2 * df / n.
#'
#' \code{cut_off = 'Curve'}, calculate L2 function norm for curves for each compositional covariates
#' , across \code{df} and \code{lambda}. Consider those betas' with L2 curve norm smaller than
#' sum(L2_j)*lower_tri as none-selected. \cr
#'
#' \code{cut_off = 'Matrix'} calculate L2 norm for coefficient matrix by rows, across
#'  \code{df} and \code{lambda}, condiser these betas' with vector L2 norm that is smaller than
#'  srqt(sum(L2_j^2)) * lower_tri.
#'
#'  \code{cut_off = 'Strict'}, strict none zero group, whenever there is none zero entries considered
#'  as selected.
#'
#'
#' @export
#'



get.GIC <- function(p, df_list, GIC_obj, GIC_arg,
                     cut_type = c("Strict", "Matrix", "Curve"), lower_tri = 0.01,
                     GIC_type = c("GIC1", "GIC2", "GIC3", "GIC4", "AIC", "BIC"),
                     method_type = c("cgl", "naive", "base"),
                     y = NULL, Zc = NULL, refit = FALSE) {
  this.call <- match.call()
  cut_type <- match.arg(cut_type)
  GIC_type <- match.arg(GIC_type)
  method_type <- match.arg(method_type)
  c_const <- 1

  df_len <- length(df_list)
  lam_len <- length(GIC_obj$lam)
  n <- nrow(GIC_obj$Funcomp.CGL.fit[[1]]$Z)

  if(missing(GIC_arg)) GIC_arg <- as.list(GIC_obj$Funcomp.CGL.fit[[1]]$call)
  if(is.null(GIC_arg$basis_fun)) GIC_arg$basis_fun <- "bs"
  if(is.null(GIC_arg$degree)) GIC_arg$degree <- 3
  if(is.null(GIC_arg$sseq)) GIC_arg$sseq = seq(0, 1, length.out = 100)

  #p <- p - as.integer(method_type == "base")
  ## now: p including all compositional predictors



  if(GIC_type %in% c("GIC2", "GIC4")) c_const = log(log(n))
  N_zero <- matrix(NA, nrow = df_len, ncol = lam_len)
  cut_norm <- array(NA, dim = c(df_len, p, lam_len))
  MSE <- GIC_obj$MSE  # matrix(NA, nrow = df_len, ncol = lam_len)
  #print(MSE)

  GIC_curve <- matrix(NA, nrow = df_len, ncol = lam_len)
  GIC_min <- vector()
  temp <- list()
  scaler <- vector()
  # calculate GIC_curve
  for(i in seq(df_len)) {
    dk = df_list[i]
    beta_path <- GIC_obj$Funcomp.CGL.fit[[as.character(dk)]]$beta[, 1:lam_len]
    # all mapped to p predictors
    if(method_type == "base") beta_path <- coef(GIC_obj$Funcomp.CGL.fit[[as.character(dk)]])[, 1:lam_len]

    if (cut_type == "Curve") {
      # GIC_arg: basis_fun, sseq, degree
      # call function Genbasis
      #B <- Genbasis(sseq = GIC_arg$sseq, df = dk, degree = GIC_arg$degree, type = GIC_arg$basis_fun)

      temp <- apply(beta_path, 2,
                    function(x, basis, p, df, lower_tri) {
                      beta_curve = basis %*% matrix(x[1:(p*df)], nrow = df)
                      # equally spaced time sequence
                      norm_L22 <- colSums(beta_curve^2)
                      norm_L22 <- norm_L22 - colSums(beta_curve[c(1, length(GIC_arg$sseq)), ]^2) / 2
                      norm_L22 <- norm_L22 * (GIC_arg$sseq[2] - GIC_arg$sseq[1])
                      # cut by L2 norm percentage
                      a = which(sqrt(norm_L22) > (sqrt(sum(norm_L22)) * lower_tri) ) # which(norm_L22 > (sum(norm_L22) * lower_tri))
                      return(a)
                    },
                    basis = Genbasis(sseq = GIC_arg$sseq, df = dk, degree = GIC_arg$degree, type = GIC_arg$basis_fun),
                    p = p, df = dk, lower_tri = lower_tri)

    } else if (cut_type == "Matrix") {

      temp <-  apply(beta_path, 2,
                     function(x, lower_tri, df) {
                       beta_C <- matrix(x[1:(p*df)], nrow = df)
                       norm_L22 <- apply(beta_C, 2, function(x) sum(x^2) )
                       a = which(sqrt(norm_L22) > (sqrt(sum(norm_L22)) * lower_tri))
                       return(a)
                     },
                     lower_tri = lower_tri, df = dk)

    } else {

      temp <- apply(beta_path, 2,
                    function(x, p, dk) {
                      B = matrix(x[1:(p*dk)], nrow = p, ncol = dk, byrow = TRUE)
                      a = apply(B, 1, function(x) max(abs(x)))
                      return(which(a>0))
                    },
                    p = p, dk = dk)

    }

    N_zero[i, ] <- sapply(temp, length)

    # re-calculate MSE with cut_off betas's
    if(refit) {
      # need user provided y and Zc

      group_index <- matrix(1:(p*dk), nrow = dk)
      for(j in which(N_zero[i, ] > 0)) beta_path[-c(group_index[, temp[[j]]], nrow(beta_path)) ,j] <- 0
      #print(beta_path)
      Xmat <-  cbind2(cbind(GIC_obj$Funcomp.CGL.fit[[as.character(dk)]]$Z, Zc), 1)
      if(method_type == "base"){
        beta_path_orig <- beta_path[-group_index[, GIC_obj$Funcomp.CGL.fit[[1]]$ref], ]
        MSE[i, ] <- apply(beta_path_orig, 2, function(x, ymat, Xmat) mean((ymat - Xmat %*% x)^2),
                          ymat = y, Xmat = Xmat)
      } else {
        MSE[i, ] <- apply(beta_path, 2, function(x, ymat, Xmat) mean((ymat - Xmat %*% x)^2),
                          ymat = y, Xmat = Xmat)
      }
    }
    #print(MSE[i, ])

    alpha = switch(GIC_type,
                   "GIC1" = log(max(p*dk,n)) * dk,
                   "GIC2" = log(max(p*dk,n)) * dk,
                   "GIC3" = log(p*dk) * dk,
                   "GIC4" = log(p*dk) * dk,
                   "AIC" = dk * 2,
                   "BIC" = dk * log(n)) / n
    scaler[i] <- alpha * c_const

    if(method_type == "naive") {
      GIC_curve[i, ] <- log(MSE[i, ]) + N_zero[i, ] * scaler[i]
    } else {
      GIC_curve[i, ] <- log(MSE[i, ]) + (ifelse(N_zero[i, ] == 0, 0, N_zero[i, ] - 1)) * scaler[i]
    }
    #print(GIC_curve[i, ])
    #### ??????? get rid of zero solution
    GIC_min[i] <- min(GIC_curve[i, which(N_zero[i, ] > 0)])
  }
  rm(temp)


  ### ??????? multiple K and lam attains minimum ??????
  #temp <- min(GIC_min)
  #print(GIC_min)
  k_loc <- which.min(GIC_min)
  lam_loc <- which(GIC_curve[k_loc, ] == GIC_min[k_loc])
  #cat("k_loc", k_loc , "\r\n")
  #cat("lam_loc", lam_loc , "\r\n")
  #print(coef(GIC_obj$Funcomp.CGL.fit[[k_loc]])[, lam_loc])
  beta_GIC <- coef(GIC_obj$Funcomp.CGL.fit[[k_loc]])[, lam_loc]
  k_opt <- df_list[k_loc]

  result <- list(beta = beta_GIC, k_opt = k_opt, lam_loc = lam_loc,
                 GIC_min = GIC_min, GIC_curve = GIC_curve,
                 N_zero = N_zero, MSE = MSE, scaler = scaler,
                 method_type = method_type, GIC_type = GIC_type, cut_type = cut_type)
  result$call <- this.call
  return(result)
}










# @title
# GIC cirterion selection for FuncompCGL
#
# @description
# Calculate GIC for compCL, return value of \code{lam}.
#
#@usage
#GIC.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL,
#                           ref = NULL,
#                           W = rep(1,times = p - length(ref)),
#                           k = 4:10, nlam = 100, outer_maxiter = 1e+6,
#                           cut_off = c("Curve","Matrix", "NULL"),
#                           lower_tri = 0.01, ...)
#
#
#
# @inheritParams cv.FuncompCGL
# @param cut_off,lower_tri cut-off None zero group by method cut_off and lower boundary.
#              -default value for cut_off is '\code{Curve}', lower_tri = 0.01
# @param \dots other arguments that could be passed to FuncompCL.
#
# @return an object of class \code{\link{GIC.FuncompCGL}} is returned.
# \item{Funcomp.CGL.fit}{a list, length of \code{k},
#                        of fitted \code{\link{FuncompCGL}} object for the full data.
#                        objects with S3 calss \code{\link{FuncompCGL}}}
# \item{lam}{the values of \code{lam} used in the fits}
# \item{loglike}{loglikelihood, a \code{k} by nlam matrix}
# \item{Nzero}{a \code{k} by nlam matrix for Nzero group cut-off by \code{cut_off} and \code{lower_tri}}
#
# @details
# $$GIC(\eqn{\lambda}) = log(MSE) + Selection * \code{alpha}$$, for normal error.
# If include linear constraints, like in log contract and constraints group lasso model,
# selection = None zero group - 1; otherwise with consideration of linear constraints,
# selection = None zero group.
#
# BIC:alpha = log(n)*df, GIC:alpha=log(log(n)) /n * log(max(p*df,n)) * df
#
# \code{cut_off = 'Curve'}, calculate L2 function norm for curves for each compositional covariates
# , across \code{df} and \code{lambda}. Consider those betas' with L2 curve norm smaller than
# sum(L2_j)*lower_tri as none-selected.
#
# \code{cut_off = 'Matrix'} calculate L2 norm for coefficient matrix by rows, across
#  \code{df} and \code{lambda}, condiser these betas' with vector L2 norm that is smaller than
#  srqt(sum(L2_j^2)) * lower_tri.
#
#  \code{cut_off = 'NULL'}, strict none zero group, whenever there is none zero entries considered
#  as selected.
#
# @examples
#
# df_beta = 5
# p = 30
# beta_C_true = matrix(0, nrow = p, ncol = df_beta)
# beta_C_true[3, ] <- c(-1, 0, 0, 0, -0.5)
# beta_C_true[1, ] <- c(1, 0, 1 , 0, -0.5)
# beta_C_true[2, ] <- c(0, 0,  -1,  0,  1)
#
# nfolds = 10
# k_list <- c(4,5)
# n_train = 100
# n_test = 500
#
# Data <- Model2(n = 100, p = p, m = 0, intercept = TRUE,
#                SNR = 2, sigma = 2,
#                rho_X = 0, rho_W = 0.5,
#                Corr_X = "CorrCS", Corr_W = "CorrAR",
#                df_W = 5, df_beta = df_beta,
#                ns = 20, obs_spar = 1, theta.add = FALSE, #c(0,0,0),
#                beta_C = as.vector(t(beta_C_true)))
# y <- drop(Data$data$y)
# n <- length(y)
# X <- Data$data$Comp
# Zc <- Data$data$Zc
# intercept <- Data$data$intercept
# m <- ifelse(is.null(Zc), 0, dim(Zc)[2]) #+ as.integer(intercept)
# m1 <- m + as.integer(intercept)
# sseq <- Data$basis.info[,1]
# beta_C.true <- matrix(Data$beta[1:(p*(df_beta))],
#                       nrow = p, ncol = df_beta, byrow = TRUE)
# beta_curve.true <- Data$basis.info[,-1] %*% t(beta_C.true)
# Non_zero.true <- (1:p)[apply(beta_C.true, 1, function(x) max(abs(x)) > 0)]
# foldid <- sample(rep(seq(nfolds), length = n))
#
# arg_list <- as.list(Data$call)[-1]
# arg_list$n <- n_test
# Test <- do.call(Model2, arg_list)
# y_test <- drop(Test$data$y)
#
# # y = y
# # X = X
# # Zc = Zc
# # intercept = intercept
# # W = rep(1, p)
# # k = k_list
# # nfolds = 10
# # trim = 0
# # tol = 0
# # inner_eps = 1e-6
# # inner_maxiter = 1E3
# # dfmax = 20
# # lambda.factor = 1e-20
# # mu_ratio = 1
# # outer_eps = 1e-6
# # keep = TRUE
# # Trange = c(0,1)
#
#
#
# GIC_m1 <- GIC.FuncompCGL( y = y, X = X, Zc = Zc, ref = NULL,
#                           k = k_list, cut_off = "Curve", lower_tri = 0.01)
#
#
#
# GIC_curve <- matrix(NA, nrow = length(k_list), ncol = length(GIC_m1$lam))
# ##### 1
# for(i in 1:length(k_list)) {
# dk = k_list[i]
# alpha <- log(log(n)) /n * log(max(p*dk,n)) * dk
# GIC_curve[i, ] <- GIC_m1$loglike[i, ] + (GIC_m1$Nzero[i, ] - 1) * alpha
#
# }
#
# plot.args = list(x = log(GIC_m1$lam),y = GIC_curve[1, ],
# ylim = range(GIC_curve),
# xlab="log(lambda)",
# ylab="GIC",
# type="n")
# do.call("plot",plot.args)
#
# for(i in 1:length(k_list)) {
#
# points(x = log(GIC_m1$lam), y = GIC_curve[i, ], col = rainbow(i))
#
# }
#
# loc  = which(GIC_curve == min(GIC_curve), arr.ind = TRUE)
#
# GIC_min <- vector()
# for(i in 1:length(k_list)) {
# # get rid of zero solution
# GIC_min[i] <- min(GIC_curve[i, (max(which(GIC_m1$Nzero[i, ] == 0))+1) : length(GIC_m1$lam)])
# }
# #cat("GIC_min:", GIC_min, "\r\n")
# k_loc <- which.min(GIC_min)
# lam_loc <- which(GIC_curve == GIC_min[k_loc])
# beta_GIC <-  GIC_m1$Funcomp.CGL.fit[[k_loc]]$beta[, lam_loc]
# k_opt <- k_list[k_loc]
#
#
#
#
# beta_C <- matrix(beta_GIC[1:(p*k_opt)], byrow = TRUE, nrow = p)
# cat("colSums:", colSums(beta_C))
# #Non.zero <- which(abs(beta_C[,1]) > 0)
# Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
# Non.zero <- (1:p)[Non.zero]
# cat("None zero groups:", Non.zero)
# #vet(beta, p = p, k = k_opt)
#
# par(mfrow=c(1,4))
# do.call("plot",plot.args)
# for(i in 1:length(k_list)) {
# points(x = log(GIC_m1$lam), y = GIC_curve[i, ], col = rainbow(i))
# text(log(GIC_m1$lam), GIC_curve[i, ], labels=paste(k_list[i]), cex= 0.7, pos=3, col = rainbow(i))
# }
# matplot(sseq, beta_curve.true,
#         ylab = "coeffcients curve", xlab = "TIME", #main = "TRUE",
#         ylim = range(Data$beta[1:(p*df_beta)]),
#         type = "l")
# abline(a = 0, b = 0, col = "grey", lwd = 2)
# title("TRUE", line = 0.5)
# text(0, beta_curve.true[1, Non_zero.true], labels = paste(Non_zero.true))
#
# B <- splines::bs(Data$basis.info[,1], df = k_opt, intercept = TRUE)
# beta_curve <- B %*% t(beta_C)
# matplot(sseq, beta_curve,
#         ylab = "coef", xlab = "TIME", #main = "ESTI",
#         ylim = range(Data$beta[1:(p*df_beta)])#,
#         #type = "l"
# )
# abline(a = 0, b = 0, col = "grey", lwd = 2)
# title("Estimate", line = 0.5)
# text(0, beta_curve[1, Non.zero], labels = paste(Non.zero))
# text(tail(sseq, 1), beta_curve[dim(beta_curve)[1], Non.zero], labels = paste(Non.zero))
# plot(apply(abs(beta_C),1,sum))
# title(paste0("k=", k_opt), line = 0.5)
# title(paste0("Method cgl"), outer=TRUE, line = -2)
# par(mfrow=c(1,1))
# ##set a cutoff when you compute nonzeros
# Non.zero <- apply(beta_C, 1, function(x)
#                  ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
# Non.zero <- (1:p)[Non.zero]
# Non.zero
#
# \dontrun{
# cgl_GIC <- list()
# # MSE <- crossprod(y -  cbind2(cbind(GIC_m1$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
# # Zc), 1) %*% beta_GIC) / length(y)
# X_train <- cbind2(cbind(GIC_m1$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z, Zc), 1)
# MSE <- sum((y -  X_train %*% beta_GIC)^2) / length(y)
# # R_sqr <- 1 - crossprod(y -  cbind2(cbind(GIC_m1$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
# # Zc), 1) %*% beta_GIC) / crossprod(y -  mean(y))
# R_sqr <- sum((y -  X_train %*% beta_GIC)^2)
# R_sqr <- 1 - R_sqr / crossprod(y -  mean(y))
#
# obj <- FuncompCGL(y = y_test, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
# # PE <- sum((Test$data$y - cbind2(cbind(obj$Z, Test$data$Zc), 1) %*% beta_GIC)^2)
# # / length(drop(Test$data$y))
# X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
# PE <- sum((y_test - X_test %*% beta_GIC)^2) / length(y_test)
# cgl_GIC$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#
# cgl_GIC$Non.zero_cut <- Non.zero
# cgl_GIC <- c(cgl_GIC,
#              ERROR_fun(beta_fit = beta_GIC, beta_true = Data$beta,
#                        basis_fit = B, basis_true = Data$basis.info[,-1],
#                        sseq = Data$basis.info[, 1],
#                        m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#              k = k_opt)
# cgl_GIC$coef <- list(beta_C = beta_C, beta_c = tail(beta_GIC, m1))
#
# }
#
# @export

# GIC.FuncompCGL2 <- function(y, X, Zc = NULL, lam = NULL, ref = NULL,
#                            W = rep(1,times = p - length(ref)),
#                            k = 4:10,
#                            nlam = 100, outer_maxiter = 1e+6,
#                            cut_off = c("Curve","Matrix", "NULL"), lower_tri = 0.01,
#                            ...) {
#   y <- drop(y)
#   n <- length(y)
#   object <- as.list(seq(length(k)))
#   names(object) <- k
#   cut_off <- match.arg(cut_off)
#   if(!is.null(lam) || length(k) == 1) {
#
#     # Case I
#     if(dim(X)[1] == n ) p = dim(X)[2] / k else p <- dim(X)[2] - 2
#
#     for(i in 1:length(k)){
#       ###cat("1", k[i], "\r\n")
#       object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam,
#                                 W = W, ref = ref,
#                                 k = k[i],
#                                 nlam = nlam, outer_maxiter = outer_maxiter,
#                                 ...)
#     }
#   } else {
#     # Case II: find commom lambda first
#     p <- ncol(X) - 2
#     ###cat(p)
#     ###cat(W)
#     ###cat("length(W)", length(W), "\r\n")
#     ###cat("is.missing(ref)", is.null(ref), "\r\n")
#     for(i in 1:length(k)){
#       ###cat("B", k[i], "\n")
#       # Caculate integral Z and W matrix (if W is a functoin)
#       object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1,
#                                 W = W, ref = ref,
#                                 k = k[i], outer_maxiter = 0,
#                                 ...)
#     }
#
#     lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k
#
#     # Solution path for each df k
#     for(i in 1:length(k)) {
#       object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
#                                 W = object[[i]]$W, ref = ref,
#                                 lam = lam0, nlam = nlam,
#                                 outer_maxiter = outer_maxiter,
#                                 k = k[i], ...)
#     }
#   }
#
#   GIC.nlam <- sapply(object, function(x) length(drop(x$lam)))
#   GIC.nlam_id <- which.min(GIC.nlam)
#   GIC.lam <- object[[GIC.nlam_id]]$lam
#   GIC.nlam <- GIC.nlam[GIC.nlam_id]
#   MSE <-matrix(NA, nrow = length(k), ncol = GIC.nlam)
#   for(i in seq(length(k))) {
#     predmat <- predict(object = object[[i]], newx = cbind(object[[i]]$Z, Zc))
#     MSE[i, ] <- apply(predmat - y, 2, function(x) mean(x^2))
#   }
#
#
#
#   ############uncompleted chunk
#   N_zero <- matrix(NA, nrow = length(k), ncol = GIC.nlam)
#   if(cut_off == "Curve") {
#     norm_L2 <- array(NA, dim = c(length(k), p, GIC.nlam))
#     argu <- as.list(object[[1]]$call)
#     argu$basis_fun <- ifelse(is.null(argu$basis_fun), "bs", argu$basis_fun)
#     argu$T.name <- ifelse(is.null(argu$T.name), "TIME", argu$T.name)
#     argu$sseq <- range(unique(sort(X[, argu$T.name])))
#     argu$sseq <- range(c(argu$Trange, argu$sseq))
#     if(!is.null(argu$interval) && startsWith(argu$interval, "S")) argu$sseq <- c(0,1)
#     argu$degree <- ifelse(is.null(argu$degree), 3, argu$degree)
#     argu$sseqs <- seq(argu$sseq[1], argu$sseq[2], length.out = 100)
#     for(i in seq(length(k))){
#       df = k[i]
#       B <- Genbasis(sseq = argu$sseqs,
#                df = df, degree = argu$degree, type = argu$basis_fun)
#       norm_L2[i, , ] <- apply(coef(object[[as.character(df)]]), 2, function(x, B, p, sseq, df)
#                               Curve_error(coef_esti = x,  p = p, k_fit = df, basis_fit = B, sseq = sseq)[(p+1):(2*p)],
#                               B = B, p = p , sseq = argu$sseqs, df = df)
#       N_zero[i, ] = apply(norm_L2[i, , ], 2, function(x, lower_tri)
#                             length(which(x > sqrt(sum(x^2)) * lower_tri)),
#                             lower_tri = lower_tri )
#     }
#
#   } else if (cut_off == "Matrix") {
#     norm_L2 <- array(NA, dim = c(length(k), p, GIC.nlam))
#     for(i in seq(length(k))){
#       df = k[i]
#
#       norm_L2[i, , ] <- apply(coef(object[[as.character(df)]]), 2, function(x, p, df) {
#                               B = matrix(x[1:(p*df)], nrow = p, ncol = df, byrow = TRUE)
#                               return(apply(B, 1, function(x) sum(x^2)))
#                               }, p = p, df = df)
#       N_zero[i, ] <-  apply(norm_L2[i, , ], 2, function(x, lower_tri)
#                           length(which( sqrt(x) > sqrt(sum(x)) * lower_tri)),
#                           lower_tri = lower_tri )
#     }
#
#   } else {
#
#     for(i in seq(length(k))) {
#       df = k[i]
#       N_zero[i, ] <- apply(coef(object[[as.character(df)]]), 2, function(x, p, df) {
#                            B = matrix(x[1:(p*df)], nrow = p, ncol = df, byrow = TRUE)
#                            a = apply(B, 1, function(x) max(abs(x)))
#                            return( length(a>0))
#                            }, p = p, df = df)
#     }
#   }
#   ############uncompleted chunk
#
#   result <- list(Funcomp.CGL.fit = object, lam = GIC.lam, loglike = MSE, Nzero = N_zero)
#   result$class <- "GIC.FuncompCGL"
#   return(result)
# }




