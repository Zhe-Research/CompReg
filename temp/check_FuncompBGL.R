library(compReg)
library(splines)
library(MASS)
library(dplyr)
library(plyr)
source("R/auxiliary.R")
source('R/tools.R')
source("integral.R")
m1 <- FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept,
              W = rep(1, p), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
              k = 4,
              #nfolds = 10, trim = 0,
              tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
              dfmax = 30, lambda.factor = 1e-3,
              mu_ratio = 1, outer_eps = 1e-6,
              #keep = TRUE,
              Trange = c(0,1)
              #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
)


m2 <- FuncompBGL(y = y, X = X, Zc = Zc, intercept = intercept,
                 ref = NULL,
                 k = 4,
                 tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
                 dfmax = 30, lambda.factor = 1e-3,
                 mu_ratio = 1, outer_eps = 1e-6,
                 #keep = TRUE,
                 Trange = c(0,1)
                 #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
)

max(abs(m1$beta - m2$beta))
max(abs(m1$lam - m1$lam))
m1$dim - m2$dim

m3 <- FuncompBGL(y = y, X = X, Zc = Zc, intercept = intercept,
                 ref = 1,
                 k = 4,
                 tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
                 dfmax = 30, lambda.factor = 1e-3,
                 mu_ratio = 1, outer_eps = 1e-6,
                 #keep = TRUE,
                 Trange = c(0,1)
                 #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
)

m4 <- FuncompBGL(y = y, X = X, Zc = Zc, intercept = intercept,
                 ref = 4,
                 k = 4,
                 tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
                 dfmax = 30, lambda.factor = 1e-3,
                 mu_ratio = 1, outer_eps = 1e-6,
                 #keep = TRUE,
                 Trange = c(0,1)
                 #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
)

m5 <- FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept,
           W = rep(1, p), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
           k = 4,
           #nfolds = 10, trim = 0,
           tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
           dfmax = 30, lambda.factor = 1e-3,
           mu_ratio = 1, outer_eps = 1e-6,
           #keep = TRUE,
           Trange = c(0,1)
           #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
)

m6 <- FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept, ref = 4,
                 W = rep(1, p - 1), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
                 k = 4,
                 #nfolds = 10, trim = 0,
                 tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
                 dfmax = 30, lambda.factor = 1e-3,
                 mu_ratio = 1, outer_eps = 1e-6,
                 #keep = TRUE,
                 Trange = c(0,1)
                 #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
)
P_bl <- c(4, 1)
a = 3
k_list = 4
object <- list()
q = 1
  for(l in 1:length(k_list)){
    object[[l]] <- m1$Z
    group_index <- matrix(1:((p)*k_list[l]), nrow = k_list[l])
    Z_ref <- object[[l]][, group_index[, P_bl[q] ]]
    object[[l]] <- object[[l]][, -group_index[, P_bl[q]]]
    group_index <- matrix(1:((p-1)*k_list[l]), nrow = k_list[l])
    for(j in 1:(p-1)) {
      object[[l]][, group_index[, j]] <- object[[l]][, group_index[, j]] - Z_ref
    }

    object[[l]] <- FuncompCGL(y = y, X = object[[l]], Zc = Zc, intercept = intercept,
                              W = rep(1, p - 1), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
                              k = k_list[l],
                              #nfolds = 10, trim = 0,
                              tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
                              dfmax = 30, lambda.factor = 1e-3,
                              mu_ratio = 1, outer_eps = 1e-6,
                              #keep = TRUE,
                              Trange = c(0,1)
                              #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
    )



      # FuncompCGL(y = y, X = object[[l]], Zc = Zc, lam = NULL, intercept = intercept,
      #                         W = rep(1, p-1), #W = function(x){ diag( 1 / apply(x, 2, sd) ) }
      #                         nlam = 1, outer_maxiter = 0, mu_ratio = 0,
      #                         k = k_list[l],
      #                         #dfmax = dfmax -1, pfmax = pfmax - 1,
      #                         tol = tol)
  }

  lam0 <- max(sapply(object, "[[", "lam"))



  for(l in 1:length(k_list)){
    object[[l]] <- FuncompCGL(y = y, X = object[[l]]$Z, Zc = Zc, intercept = intercept,
                              W = object[[l]]$W, #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
                              k = k_list[l],
                              #nfolds = 10, trim = 0,
                              tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
                              dfmax = 30, lambda.factor = 1e-3,
                              mu_ratio = 1, outer_eps = 1e-6,
                              #keep = TRUE,
                              Trange = c(0,1)
                              #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
    )


      # FuncompCGL(y = data$data$y, X = object[[l]]$Z, Zc = data$data$Zc,intercept = intercept,
      #                         W = object[[l]]$W,
      #                         nlam = nlam,   lambda.factor = lambda.factor,
      #                         inner_eps = outer_eps , inner_maxiter = outer_maxiter, mu_ratio = 0,
      #                         lam = lam0,
      #                         k = k_list[l],
      #                         dfmax = dfmax - 1, pfmax = pfmax - 1,
      #                         tol  = tol)

  }

  # cv.nlam <- sapply(object, function(x) length(drop(x$lam)))
  # cv.nlam_id <- which.min(cv.nlam)
  # cv.lam <- object[[cv.nlam_id]]$lam
  # cv.nlam <- cv.nlam[cv.nlam_id]
  #
  #
  # cvm <- matrix(NA, nrow = length(k_list), ncol = max(cv.nlam))
  # rownames(cvm) <- paste0("df=", k_list)
  # colnames(cvm) <- seq(cv.nlam)
  # cvsd <- cvm
  # cvm.trim <- cvm
  # cvsd.trim <- cvm
  #
  #
  # for(l in 1:length(k_list)) {
  #
  #   outlist <- as.list(seq(nfolds))
  #
  #   for(j in 1:nfolds) {
  #     #cat("folds", j, "\r\n")
  #     which <- foldid == j
  #     y_train <- data$data$y[!which]
  #     Z <- object[[l]]$Z
  #     Z_train <- object[[l]]$Z[!which, , drop = FALSE]
  #     Zc_train <- data$data$Zc[!which, , drop = FALSE]
  #     outlist[[j]] <- FuncompCGL(y = y_train, X = Z_train, Zc = Zc_train, k = k_list[l], lam = cv.lam, W = object[[l]]$W,
  #                                intercept = intercept,
  #                                inner_eps = outer_eps , inner_maxiter = outer_maxiter, mu_ratio = 0,
  #                                dfmax = dfmax - 1, pfmax = pfmax - 1, tol  = tol)
  #   }
  #
  #   #X <- cbind(Z, Zc)
  #   cvstuff <- cv.test(outlist, drop(data$data$y), X = cbind(object[[l]]$Z, Zc = data$data$Zc), foldid, lam = cv.lam, trim = 0)
  #   cvm[l, ] <- cvstuff$cvm
  #   cvsd[l, ] <- cvstuff$cvsd
  #
  # }
  #
  # lammin <- ggetmin(lam = drop(cv.lam), cvm = cvm, cvsd =cvsd, k_list = k_list)
  #
  # switch(rule1,
  #        "lam.1se"={
  #          s <- lammin$lam.1se[1]
  #          k_se <- lammin$lam.1se[2]
  #        },
  #        "lam.min" = {
  #          s <- lammin$lam.min[1]
  #          k_se <- lammin$lam.min[2]
  #        })
  #
  #
  # beta3 <- coef(object[[k_se - a]], s = s)
  #
  # MSE <- crossprod(data$data$y -  cbind2(cbind(object[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / n
  # R_sqr_train <- 1 - crossprod(data$data$y -  cbind2(cbind(object[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / crossprod(data$data$y -  mean(data$data$y))
  # beta3 <- c(beta3[ifelse(P_bl[q]==1,0,1):((P_bl[q]-1)*k_se)],
  #            -colSums(matrix(beta3[1:((p-1) *(k_se))], byrow = TRUE, nrow = p-1))
  #            ,beta3[((P_bl[q]-1)*k_se+1):length(beta3)])
  #
  # beta_C <- vet(beta3, p = p, k = k_se)$C
  # if(max(abs(colSums(beta_C)) > 1e-8)) cat("Column sum of coef matrix is nonn zero \r\n")
  # beta_c <- vet(beta3, p = p, k = k_se)$b
  # basis_fit <- Genbasis(sseq = sseq_complete, df = k_se, degree = 3, type = "bs")
  #
  #
  # m2.8 <- FuncompCGL(y = data2$data$y, X = data2$data$Comp,  Zc = data2$data$Zc, intercept = intercept, k = k_se,
  #                    nlam = 1, outer_maxiter = 0)
  #
  # PE <- sum((data2$data$y - cbind2(cbind(m2.8$Z, data2$data$Zc), 1) %*% beta3)^2) / length(drop(data2$data$y))
  # #R_sqr_test <- 1 - sum((data2$data$y - cbind2(cbind(m2.8$Z, data2$data$Zc), 1) %*% beta3)^2) / crossprod(data2$data$y -  mean(data2$data$y))
  #
  #
  #
  #
  # Baseline.1se[[q]][[1]] <- list(pred_error = c(MSE = MSE, PE = PE, Rsqr_train = R_sqr_train))
  #
  # Baseline.1se[[q]][[1]] <- c(Baseline.1se[[q]][[1]],
  #                             ERROR_fun(beta_fit = beta3, beta_true = data$beta,
  #                                       basis_fit = basis_fit, basis_true = basis_true,
  #                                       sseq = sseq_complete,
  #                                       m = m, p = p, Nzero_group = Nzero_group),
  #                             k = k_se)
  # Baseline.1se[[q]][[1]]$coef <- list(beta_C = beta_C, beta_c = beta_c)
  #
  # Non.zero <- apply(beta_C, 1, function(x) ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
  # Non.zero <- (1:p)[Non.zero]
  # group_index <- matrix(1:(p*k_se), nrow = k_se)
  # beta3[group_index[, -Non.zero]] <- 0
  # MSE <- crossprod(data$data$y -  cbind2(cbind(cvm2.6$Funcomp.CGL.fit[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / n
  # R_sqr_train <- 1 - crossprod(data$data$y -  cbind2(cbind(cvm2.6$Funcomp.CGL.fit[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / crossprod(data$data$y -  mean(data$data$y))
  # PE <- sum((data2$data$y - cbind2(cbind(m2.8$Z, data2$data$Zc), 1) %*% beta3)^2) / length(drop(data2$data$y))
  # Baseline.1se[[q]][[2]] <- list(pred_error = c(MSE = MSE, PE = PE, Rsqr_train = R_sqr_train))
  #
  # Baseline.1se[[q]][[2]] <- c(Baseline.1se[[q]][[2]],
  #                             ERROR_fun(beta_fit = beta3, beta_true = data$beta,
  #                                       basis_fit = basis_fit, basis_true = basis_true,
  #                                       sseq = sseq_complete,
  #                                       m = m, p = p, Nzero_group = Nzero_group),
  #                             k = k_se)
  #
  # switch(rule2,
  #        "lam.1se"={
  #          s <- lammin$lam.1se[1]
  #          k_se <- lammin$lam.1se[2]
  #        },
  #        "lam.min" = {
  #          s <- lammin$lam.min[1]
  #          k_se <- lammin$lam.min[2]
  #        })
  #
  #
  # beta3 <- coef(object[[k_se - a]], s = s)
  #
  # MSE <- crossprod(data$data$y -  cbind2(cbind(object[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / n
  # R_sqr_train <- 1 - crossprod(data$data$y -  cbind2(cbind(object[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / crossprod(data$data$y -  mean(data$data$y))
  # beta3 <- c(beta3[ifelse(P_bl[q]==1,0,1):((P_bl[q]-1)*k_se)],
  #            -colSums(matrix(beta3[1:((p-1) *(k_se))], byrow = TRUE, nrow = p-1))
  #            ,beta3[((P_bl[q]-1)*k_se+1):length(beta3)])
  #
  # beta_C <- vet(beta3, p = p, k = k_se)$C
  # if(max(abs(colSums(beta_C)) > 1e-8)) cat("Column sum of coef matrix is nonn zero \r\n")
  # beta_c <- vet(beta3, p = p, k = k_se)$b
  # basis_fit <- Genbasis(sseq = sseq_complete, df = k_se, degree = 3, type = "bs")
  #
  #
  # m2.8 <- FuncompCGL(y = data2$data$y, X = data2$data$Comp,  Zc = data2$data$Zc, intercept = intercept, k = k_se,
  #                    outer_maxiter = 0, nlam = 1)
  #
  # PE <- sum((data2$data$y - cbind2(cbind(m2.8$Z, data2$data$Zc), 1) %*% beta3)^2) / length(drop(data2$data$y))
  #
  # Baseline.min[[q]][[1]] <- list(pred_error = c(MSE = MSE, PE = PE, Rsqr_train = R_sqr_train))
  #
  # Baseline.min[[q]][[1]] <- c(Baseline.min[[q]][[1]],
  #                             ERROR_fun(beta_fit = beta3, beta_true = data$beta,
  #                                       basis_fit = basis_fit, basis_true = basis_true,
  #                                       sseq = sseq_complete,
  #                                       m = m, p = p, Nzero_group = Nzero_group),
  #                             k = k_se)
  # Baseline.min[[q]][[1]]$coef <- list(beta_C = beta_C, beta_c = beta_c)
  #
  # Non.zero <- apply(beta_C, 1, function(x) ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
  # Non.zero <- (1:p)[Non.zero]
  # group_index <- matrix(1:(p*k_se), nrow = k_se)
  # beta3[group_index[, -Non.zero]] <- 0
  # MSE <- crossprod(data$data$y -  cbind2(cbind(cvm2.6$Funcomp.CGL.fit[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / n
  # R_sqr_train <- 1 - crossprod(data$data$y -  cbind2(cbind(cvm2.6$Funcomp.CGL.fit[[k_se - a]]$Z, data$data$Zc), 1) %*% beta3) / crossprod(data$data$y -  mean(data$data$y))
  # PE <- sum((data2$data$y - cbind2(cbind(m2.8$Z, data2$data$Zc), 1) %*% beta3)^2) / length(drop(data2$data$y))
  # Baseline.min[[q]][[2]] <- list(pred_error = c(MSE = MSE, PE = PE, Rsqr_train = R_sqr_train))
  #
  # Baseline.min[[q]][[2]] <- c(Baseline.min[[q]][[2]],
  #                             ERROR_fun(beta_fit = beta3, beta_true = data$beta,
  #                                       basis_fit = basis_fit, basis_true = basis_true,
  #                                       sseq = sseq_complete,
  #                                       m = m, p = p, Nzero_group = Nzero_group),
  #                             k = k_se)

}



object[[1]]$dim
m4$dim

object[[1]]$beta - m4$beta
m4$Z - object[[1]]$Z



k = 4
try_base <- new.env()
try_base$ref <- 4
try_base$p <-  p - 1
try_base$p1 <-  try_base$p * k
try_base$X <- X
try_base$X[, X.names] <- apply(X[, X.names], 2, function(x, ref) x - ref, ref = X[, X.names[try_base$ref], drop = TRUE])
try_base$D <- split(try_base$X[, c(T.name, X.names[-try_base$ref])], try_base$X[, ID.name])
try_base$Z <- matrix(NA, nrow = n, ncol = try_base$p1)
for(i in 1:n) try_base$Z[i, ] <- ITG(try_base$D[[i]], basis, sseq, T.name, interval, insert, method)$integral

try_full <- new.env()
try_full$X <- X
try_full$D <- split(try_full$X[, c(T.name, X.names)], try_full$X[, ID.name])
try_full$Z <- matrix(NA, nrow = n, ncol =  p * k)
for(i in 1:n) try_full$Z[i, ] <- ITG(try_full$D[[i]], basis, sseq, T.name, interval, insert, method)$integral

try_base$object <- m1$Z
try_base$group_index <-  matrix(1:((p)*k), nrow = k)
try_base$Z_ref <- try_base$object[, try_base$group_index[, try_base$ref]]
try_base$object <- try_base$object[, -try_base$group_index[, try_base$ref]]
try_base$group_index <- matrix(1:((p-1)*k), nrow = k)
for(j in 1:(p-1)) {
  #cat("j", j, "\r\n")
  try_base$object[, try_base$group_index[, j]] <- try_base$object[, try_base$group_index[, j]] - try_base$Z_ref
}

# Z <- matrix(NA, nrow = n, ncol = p1)
# for(i in 1:n) Z[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral

#### LOOK in to D
Sub_j <- sample(1:n, 1)
cat("Subject:", Sub_j, "\r\n")
within(try_base$D[[Sub_j]], rm(TIME)) -
  (try_full$D[[Sub_j]][, X.names[-try_base$ref]] - try_full$D[[Sub_j]][, X.names[try_base$ref]#, drop = FALSE                                                                                          ])
                                                                       try_base$D[[Sub_j]]['TIME'] - try_full$D[[Sub_j]]['TIME']


                                                                       #### LOOK into Z
                                                                       max(abs(try_full$Z - m1$Z))
                                                                       max(abs(try_base$object - object[[1]]$Z))
                                                                       max(abs(try_full$Z[Sub_j, -(13:16)] - try_full$Z[Sub_j, (13:16)] - try_base$Z[Sub_j, ]))




                                                                       #### LOOK into object

