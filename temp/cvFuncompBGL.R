cv.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL,
                          W = rep(1,times = p - length(ref)),
                          ref = NULL,
                          k = 4:10,
                          foldid, nfolds = 10, nlam = 100,
                          trim = 0,
                          outer_maxiter = 1e+6, keep = FALSE,
                          ...) {
  y <- drop(y)
  n <- length(y)

  if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else nfolds <- max(foldid)
  if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")


  object <- as.list(seq(length(k))) # list of data for different k
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
    ###cat(p)
    ###cat(W)
    ###cat("length(W)", length(W), "\r\n")
    ###cat("is.missing(ref)", is.null(ref), "\r\n")
    for(i in 1:length(k)){
      ###cat("B", k[i], "\n")
      # Caculate integral Z and W matrix (if W is a functoin)
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1,
                                W = W, ref = ref,
                                k = k[i], outer_maxiter = 0,
                                ...)
    }

    lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k

    # Solution path for each df k
    for(i in 1:length(k)) {
      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                k = k[i], ...)
    }
  }



  cv.nlam <- sapply(object, function(x) length(drop(x$lam))) # different stopping lambda sequence by dfmax/pfmax
  cv.nlam_id <- which.min(cv.nlam)
  cv.lam <- object[[cv.nlam_id]]$lam # shared lambda sequence for different k
  cv.nlam <- cv.nlam[cv.nlam_id]


  cvm <- matrix(NA, nrow = length(k), ncol = max(cv.nlam))
  rownames(cvm) <- paste0("df=", k)
  colnames(cvm) <- seq(cv.nlam)
  cvsd <- cvm

  # deleting later
  if(trim > 0) {
    cvm.trim <- cvm
    cvsd.trim <- cvm
  }

  if(keep) {
    preval <- array(NA, dim = c(n, cv.nlam, length(k)))
  }

  # Cross-validatoin via different folders
  for(l in 1:length(k)) {

    outlist <- as.list(seq(nfolds))

    for(i in 1:nfolds) {
      ###cat("folds", i, "\r\n")
      which <- foldid == i
      y_train <- y[!which]
      Z <- object[[l]]$Z
      Z_train <- Z[!which, , drop = FALSE]
      Zc_train <- Zc[!which, , drop = FALSE]
      outlist[[i]] <- FuncompCGL(y = y_train, X = Z_train, Zc = Zc_train,
                                 k = k[l],
                                 lam = cv.lam, W = object[[l]]$W, ref = ref,
                                 outer_maxiter = outer_maxiter,
                                 ...)
    }


    cvstuff <- cv.test(outlist, y, X = cbind(Z, Zc), foldid, lam = cv.lam, trim = trim, keep = keep) #define diffrent cv.test for GLM
    cvm[l, ] <- cvstuff$cvm
    cvsd[l, ] <- cvstuff$cvsd
    if(keep) preval[, , l] <- cvstuff$fit.preval

    # delet later
    if(trim > 0) {
      cvm.trim[l, ] <- cvstuff$cvmtrim
      cvsd.trim[l, ] <- cvstuff$cvsdtrim
    }

  }


  # Select lambda
  Ftrim = list(cvm = cvm, cvsd = cvsd)
  Ftrim$cvup = Ftrim$cvm + Ftrim$cvsd
  Ftrim$cvlo = Ftrim$cvm - Ftrim$cvsd
  lammin <- ggetmin(lam = cv.lam, cvm = Ftrim$cvm, cvsd = Ftrim$cvsd, k_list = k)
  Ftrim <- c(Ftrim, lammin)

  # Delect later
  if(trim > 0) {
    Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim)
    Ttrim$cvup = Ttrim$cvm + Ttrim$cvsd
    Ttrim$cvlo = Ttrim$cvm - Ttrim$cvsd
    lammin <- ggetmin(lam = cv.lam, cvm = Ttrim$cvm, cvsd = Ttrim$cvsd, k_list = k)
    Ttrim <- c(Ttrim, lammin)
  } else {
    Ttrim <- NULL
  }

  result <- list(Funcomp.CGL.fit = object,
                 lam = cv.lam,
                 Ftrim = Ftrim,
                 Ttrim = Ttrim)
  if(keep) result <- c(result, list(fit.preval = preval, foldid = foldid))
  class(result) <- "cv.FuncompCGL"
  return(result)
}


# cvm2 <- cv.FuncompCGL(y = y, X = X , Zc = Zc, intercept = intercept,
#                       mu_ratio = 0, ref = sample(1:30, 1),
#                       k = 4:6, basis_fun = "bs",
#                       insert = "FALSE", method = "t"
#                       )


# cvm3 <- cv.FuncompCGL(y = y, X = X , Zc = Zc, intercept = intercept, mu_ratio = 0,
#                       k = 4:5, basis_fun = "bs", ref = NULL,
#                       insert = "FALSE", method = "t",
#                       #dfmax = p, tol = 1e-6, #W = rep(1, p-1)
#                       inner_eps = 1, outer_eps = 1
# )
