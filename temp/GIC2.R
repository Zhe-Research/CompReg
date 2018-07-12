GIC.FuncompCGL2 <- function(y, X, Zc = NULL, ref = NULL,
                            lam = NULL, nlam = 100,
                            W = rep(1,times = p - length(ref)),
                            k = 4:10, outer_maxiter = 1e+6,...) {
  y <- drop(y)
  n <- length(y)
  object <- as.list(seq(length(k)))
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
  return(result)
}


get.GIC2 <- function(p, df_list, GIC_obj, GIC_arg,
                     cut_type = c("Strict", "Matrix", "Curve"), lower_tri = 0.01,
                     GIC_type = c("GIC1", "GIC2", "GIC3", "GIC4", "AIC", "BIC"),
                     method_type = c("cgl", "naive", "base"),
                     y = NULL, Zc = NULL, refit = FALSE) {

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
                    basis = Genbasis2(sseq = GIC_arg$sseq, df = dk, degree = GIC_arg$degree, type = GIC_arg$basis_fun),
                    p = p, df = dk, lower_tri = lower_tri)

    } else if (cut_type == "Matrix") {

      temp <-  apply(beta_path, 2,
                     function(x, lower_tri) {
                       beta_C <- matrix(x[1:(p*df)], nrow = df)
                       norm_L22 <- apply(beta_C, 2, function(x) sum(x^2) )
                       a = which(sqrt(norm_L22) > (sqrt(sum(norm_L22)) * lower_tri))
                       return(a)
                     },
                     lower_tri = lower_tri)

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
      for(j in which(N_zero[i, ] > 0)) beta_path[-c(group_index[, A[[j]]], nrow(beta_path)) ,j] <- 0
      if(method_type == "base") beta_path_orig <- beta_path[-group_index[, GIC_obj$Funcomp.CGL.fit[[1]]$ref], ]
      Xmat <-  cbind2(cbind(GIC_obj$Funcomp.CGL.fit[[as.character(dk)]]$Z, Zc), 1)
      MSE[i, ] <- apply(beta_path_orig, 2, function(x, ymat, Xmat) mean((ymat - Xmat %*% x)^2),
                        ymat = y, Xmat = Xmat)
    }


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

    #### ??????? get rid of zero solution
    GIC_min[i] <- min(GIC_curve[i, which(N_zero[i, ] > 0)])
  }
  rm(temp)


  ### ??????? multiple K and lam attains minimum ??????
  #temp <- min(GIC_min)
  k_loc <- which.min(GIC_min)
  lam_loc <- which(GIC_curve[k_loc, ] == GIC_min[k_loc])
  beta_GIC <- GIC_obj$Funcomp.CGL.fit[[k_loc]]$beta[, lam_loc]
  k_opt <- k_list[k_loc]

  result <- list(beta = beta_GIC, k_opt = k_opt, lam_loc = lam_loc,
                 GIC_min = GIC_min, GIC_curve = GIC_curve,
                 N_zero = N_zero, MSE = MSE, scaler = scaler,
                 method_type = method_type, GIC_type = GIC_type, cut_type = cut_type)
}


Genbasis2 <- function(sseq, df, degree, type = c("bs", "fourier", "OBasis")) {

  type <- match.arg(type)
  interval <- range(sseq)

  nknots <- df - (degree + 1)
  if(nknots > 0) {
    knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
  } else {
    knots <- NULL
  }

  basis <- switch(type,
                  "bs" = splines::bs(x = sseq, df =  df, degree = degree,
                                     Boundary.knots = interval, intercept = TRUE),

                  "fourier" = fda::eval.basis(sseq,
                                              basisobj = fda::create.fourier.basis(rangeval = interval, nbasis = df)
                  ),

                  "OBasis" = orthogonalsplinebasis::evaluate(
                    orthogonalsplinebasis::OBasis(orthogonalsplinebasis::expand.knots(c(interval[1], knots, interval[2])),
                                                  order = degree + 1),
                    sseq)
  )
  return(basis)
}
