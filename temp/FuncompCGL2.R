FuncompCGL2 <- function(y, X, Zc = NULL, intercept = TRUE, k,
         T.name = "TIME", ID.name = "Subject_ID",
         degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
         insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
         interval = c("Original", "Standard"), Trange,
         W = rep(1, times = p), #pf = rep(1, times = p),
         dfmax = p, pfmax = min(dfmax * 1.5, p),
         lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
         tol = 0, mu_ratio = 1.01,
         outer_maxiter = 1e+08, outer_eps = 1e-8,
         inner_maxiter = 1e+4, inner_eps = 1e-8
         #, ...
) {

  n <- length(y)
  this.call <- match.call()


  basis_fun <- match.arg(basis_fun)
  method <- match.arg(method)
  insert <- match.arg(insert)
  interval <- match.arg(interval)

  if(dim(X)[1] == n) {
    # already take integral
    Z <- X
    p1 <- dim(X)[2]
    p <- p1 / k
  } else {
    # take integral
    X <- as.data.frame(X)
    X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
    p <- length(X.names)
    if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
    p1 <- p * k
    #if(missing(Time)) Time <- X[, T.name]
    Time <- X[, T.name]
    #if(missing(Subject_ID))
    Subject_ID <- unique(X[, ID.name])
    X[, ID.name] <- factor(X[, ID.name], levels = Subject_ID )

    if(missing(Trange)) Trange <- range(Time)
    if(interval == "Standard") {
      interval <- c(0,1)
      X[, T.name] <- (Time - min(Trange[1])) / diff(Trange)
      #X[, T.name] <- (Time - min(Time)) / diff(range(Time))
    } else {
      interval <- Trange
    }

    sseq <- round(sort(unique(c(Trange, as.vector(X[, T.name])))), 10)
    if(insert != FALSE) sseq <- round(seq(from = interval[1], to = interval[2], by = min(diff(sseq))/20), 10)



    nknots <- k - (degree + 1)
    if(nknots > 0) {
      knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
    } else {
      knots <- NULL
    }

    basis <- switch(basis_fun,

                    "bs" = bs(x = sseq, df = k, degree = degree,
                              Boundary.knots = interval, intercept = TRUE),

                    "fourier" = eval.basis(sseq,
                                           basisobj = create.fourier.basis(rangeval = interval, nbasis = k )),

                    "OBasis" = evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
                                               order = degree + 1),
                                        sseq)

    )


     cat("1")


    X[, X.names] <- log(X[, X.names]/ rowSums(X[, X.names]))
    D <- split(X[, c(T.name, X.names)], X[, ID.name])

    Z <- matrix(NA, nrow = n, ncol = p1)
    for(i in 1:n) {
      #cat(i)
      d <- D[[i]]
      Z[i, ] <- ITG(d, basis, sseq, T.name, interval, insert, method)$integral

      #Z <- sapply(D, function(x, basis, sseq, T.name, interval, insert, method)
      #           ITG(x, basis, sseq, T.name, interval, insert, method)
      #          ,sseq = sseq, basis = basis, T.name = T.name, interval = interval, insert = insert, method = method)
      #Z <- t(Z)

    }



  }


  Z <- as.matrix(Z)
  group.index <- matrix(1:p1, nrow = k)
  if( is.vector(W) ) {
    if(length(W) != p) stop("W shoulde be a vector of length p")
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimensions of Weights matrix')
  } else if( is.function(W) ){
    W_fun <- W
    W <- diag(p1)
    for(i in 1:p) W[group.index[, i], group.index[, i]] <- W_fun(Z[, group.index[, i]])
  }

  cat(dim(W))
  cat(length(W))
  output <- cglasso(y = y, Z = Z, Zc = Zc, k = k, W = W, intercept = intercept,
                    mu_ratio = mu_ratio,
                    lam = lam, nlam = nlam, lambda.factor = lambda.factor,
                    dfmax = dfmax, pfmax = pfmax,
                    tol = tol,
                    outer_maxiter = outer_maxiter, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, inner_eps = inner_eps
                    #, ...
  )
  output$Z <- Z
  output$W <- W
  output$call <- this.call
  class(output) <- "FuncompCGL"
  output$dim <- dim(output$beta)
  return(output)
}



cglasso2 <- function(y, Z, Zc = NULL, k,
         #W_fun = NULL,
         W = rep(1, times = p),
         intercept = TRUE,
         A =  kronecker(matrix(1, ncol = p), diag(k)), b = rep(0, times = k),
         u = 1, mu_ratio = 1.01,
         lam = NULL, nlam = 100,lambda.factor = ifelse(n < p1, 0.05, 0.001),
         dfmax = p, pfmax = min(dfmax * 1.5, p),
         #pf = rep(1, times = p),
         tol = 0,
         outer_maxiter = 3e+08, outer_eps = 1e-8,
         inner_maxiter = 1e+6, inner_eps = 1e-8
         #,estar = 1 + 1e-8
) {


  y <- drop(y)
  Z <- as.matrix(Z)
  n <- length(y)
  p1 <- dim(Z)[2]
  p <- p1 / k
  group.index <- matrix(1:p1, nrow = k)



  inter <- as.integer(intercept)
  Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }
  beta_c <- as.vector(Zc_proj %*% y)
  beta.ini <- c(rep(0, times = p1), beta_c)

  if( is.vector(W) ){
    if(length(W) != p) stop("W should be a vector of length p")
    W_inver <- diag(p1)
    pf <- W
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimension of Weights matrix')
    pf <- rep(1, times = p)
    W_inver <- W
  }

  Z <- Z %*% W_inver
  cat("A \r\n")
  cat(dim(A), "\r\n")
  A <- A %*% W_inver


  if(is.null(lam)) {
    lam0 <- lam.ini(Z = Z, y = y - Zc %*% beta_c, ix = group.index[1, ], iy = group.index[k ,], pf = pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  } else if(length(lam) == 1) {
    lam <- exp(seq(from=log(lam), to=log(lambda.factor * lam),length=nlam))
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

  output <- ALM_GMD(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini, lambda = lam,
                    pf = pf, dfmax = dfmax, pfmax = pfmax, A = A, b = b,
                    group_index = group.index, u_ini = u, mu_ratio = mu_ratio,
                    inner_eps = inner_eps, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter, tol = tol
                    #, estar = estar
  )

  output$beta[1:p1, ] <-  W_inver %*% output$beta[1:p1, ]
  return(output)

}
