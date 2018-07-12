FuncompCGL <- function(y, X, Zc = NULL, intercept = TRUE, ref = NULL,
                        k, degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
                        insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
                        interval = c("Original", "Standard"), Trange,
                        T.name = "TIME", ID.name = "Subject_ID",
                        W = rep(1,times = p - length(ref)),
                        dfmax = p - length(ref), pfmax = min(dfmax * 1.5, p - length(ref)),
                        lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
                        tol = 0, mu_ratio = 1.01,
                        outer_maxiter = 1e+08, outer_eps = 1e-8,
                        inner_maxiter = 1e+4, inner_eps = 1e-8) {

  n <- length(y)
  this.call <- match.call()
  basis_fun <- match.arg(basis_fun)
  method <- match.arg(method)
  insert <- match.arg(insert)
  interval <- match.arg(interval)

  if(dim(X)[1] == n) {
    # Case I: already take integral
    ## 1.
    Z <-  as.matrix(X)
    p1 <- ncol(Z)
    p <- p1 / k
    ## 2.Subtracting baseline integral
    if(is.numeric(ref)) {
      if( !(ref %in% 1:p) ) stop("Reference variable is out of range")
      #mu_ratio = 0
      group.index <- matrix(1:p1, nrow = k)
      Z_ref <- Z[, group.index[, ref]]
      Z <- Z[, -group.index[, ref]]
      for(j in 1:(p-1)) Z[, group.index[, j]] <- Z[, group.index[, j]] - Z_ref
      p1 <- ncol(Z)
    }
  } else {
    # Case II: take integral
    X <- as.data.frame(X)
    X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
    if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
    X[, X.names] <- log(X[, X.names] / rowSums(X[, X.names]))
    p <- length(X.names)
    if(is.numeric(ref) && !(ref %in% 1:p)) stop("Reference variable is out of range")
    # Time range provided in case that sample do not cover the entire range
    if(missing(Trange)) Trange <- range(X[, T.name])
    switch(interval,
           "Standard" = {
             #mapping time sequence on to [0,1]
             X[, T.name] = (X[, T.name] - min(Trange[1])) / diff(Trange)
             interval = c(0, 1)
           },
           "Original" = {
             interval = Trange
           })



    # In case that X is not represented in numerical order of Subject_ID
    X[, ID.name] <- factor(X[, ID.name], levels = unique(X[, ID.name]))

    # generate discrete basis matrix
    sseq <- round(sort(unique(c(Trange, as.vector(X[, T.name])))), 10)
    if(insert != "FALSE") sseq <- round(seq(from = interval[1], to = interval[2],
                                            #by = length(sseq) * 2,
                                            by = min(diff(sseq))/10), #### digit = 10
                                        10)

    # generate knots equally
    nknots <- k - (degree + 1)
    if(nknots > 0) {
      knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
    } else  knots <- NULL


    basis <- switch(basis_fun,
                    "bs" = bs(x = sseq, df = k, degree = degree,
                              Boundary.knots = interval, intercept = TRUE),
                    "fourier" = eval.basis(sseq,
                                           basisobj = create.fourier.basis(rangeval = interval, nbasis = k )),
                    "OBasis" = evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
                                               order = degree + 1),
                                        sseq))




    X[, X.names] <- apply(X[, X.names], 2, function(x, ref) x - ref,
                          ref = if(is.null(ref)) 0 else X[, X.names[ref], drop = TRUE])
    D <- split(X[, c(T.name, X.names[-ifelse(is.null(ref), p + 1, ref)])], X[, ID.name])
    p1 <- (p - length(ref)) * k
    Z <- matrix(NA, nrow = n, ncol = p1)
    for(i in 1:n) Z[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral

  }



  group.index <- matrix(1:p1, nrow = k)
  if( is.vector(W) ) {
    ###cat(" ,length(W) = ", length(W))
    if(length(W) != (p1 / k) ) stop("W shoulde be a vector of length p=", p1/ k)
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimensions of Weights matrix')
  } else if( is.function(W) ){
    W_fun <- W
    W <- diag(p1)
    for(i in 1:(p1/k)) W[group.index[, i], group.index[, i]] <- W_fun(Z[, group.index[, i]])
  }

  ###cat(" ,is.null(ref) = ", is.null(ref))
  ###cat(" ,dfmax = ", dfmax, "\r\n")
  output <- cglasso(y = y, Z = Z, Zc = Zc, k = k, W = W, intercept = intercept,
                    mu_ratio = mu_ratio,
                    lam = lam, nlam = nlam, lambda.factor = lambda.factor,
                    dfmax = dfmax, pfmax = pfmax,
                    tol = tol,
                    outer_maxiter = outer_maxiter, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, inner_eps = inner_eps)
  output$Z <- Z
  output$W <- W
  output$call <- this.call
  output$ref <- ref
  class(output) <- "FuncompCGL"
  output$dim <- dim(output$beta)
  return(output)
}










FuncompBGL <- function(y, X, Zc = NULL, intercept = TRUE, ref = NULL,
                       k, degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
                       insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
                       interval = c("Original", "Standard"), Trange,
                       T.name = "TIME", ID.name = "Subject_ID",
                       W = rep(1,times = p),
                       dfmax = p, pfmax = min(dfmax * 1.5, p),
                       lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
                       tol = 0, mu_ratio = 1.01,
                       outer_maxiter = 1e+08, outer_eps = 1e-8,
                       inner_maxiter = 1e+4, inner_eps = 1e-8
) {

  n <- length(y)

  this.call <- match.call()
  basis_fun <- match.arg(basis_fun)
  method <- match.arg(method)
  insert <- match.arg(insert)
  interval <- match.arg(interval)

  if(dim(X)[1] == n) {
    # Case I: already take integral
    Z <- X
    p1 <- dim(X)[2]
    p <- p1 / k
  } else {
    # Case II: take integral
    X <- as.data.frame(X)
    X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
    if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
    X[, X.names] <- log(X[, X.names] / rowSums(X[, X.names]))
    p <- length(X.names)
    #p1 <- p * k   ###?????

    Time <- X[, T.name]
    #Time range provided in case that sub-sample do not cover the entire range
    if(missing(Trange)) Trange <- range(Time)
    if(interval == "Standard") {
      interval <- c(0,1)
      #mapping time sequence on to [0,1]
      X[, T.name] <- (Time - min(Trange[1])) / diff(Trange) #X[, T.name] <- (Time - min(Time)) / diff(range(Time))
    } else {
      interval <- Trange
    }

    #In case that X is not represented in numerical order of Subject_ID
    Subject_ID <- unique(X[, ID.name])
    X[, ID.name] <- factor(X[, ID.name], levels = Subject_ID )



    #generate discrete basis matrix
    sseq <- round(sort(unique(c(Trange, as.vector(X[, T.name])))), 10)
    if(insert != "FALSE") sseq <- round(seq(from = interval[1], to = interval[2],
                                            #by = length(sseq) * 2,
                                            by = min(diff(sseq))/20),
                                            10)

    nknots <- k - (degree + 1)
    if(nknots > 0) {
      knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
    } else  knots <- NULL


    basis <- switch(basis_fun,
                    "bs" = bs(x = sseq, df = k, degree = degree,
                              Boundary.knots = interval, intercept = TRUE),
                    "fourier" = eval.basis(sseq,
                                           basisobj = create.fourier.basis(rangeval = interval, nbasis = k )),
                    "OBasis" = evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
                                               order = degree + 1),
                                        sseq)
                    )




    if(is.numeric(ref)) {
      ###cat("ref = ", ref)
      ###cat(" ,p = ", p)
      if(! ref %in% 1:p) stop("Reference variable is out of range")
      p <-  p - 1
      p1 <-  p * k
      D_ref <- split(X[, c(T.name, X.names[ref])], X[, ID.name])
      Z_ref <- matrix(NA, nrow = n, ncol = k)
      for(i in 1:n)  Z_ref[i, ] <- ITG(D_ref[[i]], basis, sseq, T.name, interval, insert, method)$integral

      X[, X.names] <- apply(X[, X.names], 2, function(x, ref) x - ref, ref = X[, X.names[ref], drop = TRUE])
      D <- split(X[, c(T.name, X.names[-ref])], X[, ID.name])

    } else {
      ###cat("ref = ", ref)
      ###cat(" ,p = ", p)
      D <- split(X[, c(T.name, X.names)], X[, ID.name])
      p1 <-  p * k
    }

    #D <- split(X[, c(T.name, X.names[-ifelse(is.null(ref), p + 1, ref)])], X[, ID.name])
    Z <- matrix(NA, nrow = n, ncol = p1)
    for(i in 1:n) Z[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral


    # D <- split(X[, c(T.name, X.names)], X[, ID.name])
    # Z <- matrix(NA, nrow = n, ncol = length(X.names)*k)
    # for(i in 1:n) Z[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral
    # for(j in (1:p)[-ref]) Z[, ((j-1)*k + 1) : (j*k)] <- Z[, ((j-1)*k + 1) : (j*k)] - Z[, ((ref*k - k +1):(ref*k))]
    # Z <- Z[, -((ref*k - k +1):(ref*k))]

  }

  ###cat(", p = ", p)
  Z <- as.matrix(Z)
  group.index <- matrix(1:p1, nrow = k)
  if( is.vector(W) ) {
    ###cat(" ,length(W) = ", length(W))
    if(length(W) != p) stop("W shoulde be a vector of length p=", p)
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimensions of Weights matrix')
  } else if( is.function(W) ){
    W_fun <- W
    W <- diag(p1)
    for(i in 1:p) W[group.index[, i], group.index[, i]] <- W_fun(Z[, group.index[, i]])
  }

  ###cat(" ,is.null(ref) = ", is.null(ref))
  ###cat(" ,dfmax = ", dfmax, "\r\n")
  output <- cglasso(y = y, Z = Z, Zc = Zc, k = k, W = W, intercept = intercept,
                    mu_ratio = mu_ratio,
                    lam = lam, nlam = nlam, lambda.factor = lambda.factor,
                    dfmax = dfmax, pfmax = pfmax,
                    tol = tol,
                    outer_maxiter = outer_maxiter, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, inner_eps = inner_eps)
  output$Z <- Z
  output$W <- W
  output$call <- this.call
  output$ref <- ref
  class(output) <- "FuncompCGL"
  output$dim <- dim(output$beta)
  return(output)
}




# m2 <- FuncompBGL(y = y, X = X , Zc = Zc, intercept = intercept, ref = NULL,
#                  mu_ratio = 0,
#                  k = 4, basis_fun = "bs",
#                  insert = "FALSE", method = "t",
#                  #dfmax = p,
#                  tol = 1e-6)
#
# m3 <- FuncompBGL(y = y, X = X , Zc = Zc, intercept = intercept, ref = 3,
#                  mu_ratio = 0,
#                  k = 4, basis_fun = "bs",
#                  insert = "FALSE", method = "t",
#                  #dfmax = p,
#                  tol = 1e-6)
