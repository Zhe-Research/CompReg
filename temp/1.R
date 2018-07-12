f1 <- function(W=rep(1,p),X, dfmax = ifelse(is.null(ref), p, p-1), ref = NULL ){
  #p = ncol(X) - 2
  this.call <- match.call()

  X <- as.data.frame(X)
  X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
  if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
  X[, X.names] <- log(X[, X.names] / rowSums(X[, X.names]))
  p <- length(X.names)
  if(is.numeric(ref)) {
    cat("ref = ", ref)
    cat(" ,p = ", p)
    if(! ref %in% 1:p) stop("Reference variable is out of range")
    p <-  p - 1
    p1 <-  p * k
    X[, X.names] <- apply(X[, X.names], 2, function(x, ref) x - ref, ref = X[, X.names[ref], drop = TRUE])
    D <- split(X[, c(T.name, X.names[-ref])], X[, ID.name])
  } else {
    cat("ref = ", ref)
    cat(" ,p = ", p)
    D <- split(X[, c(T.name, X.names)], X[, ID.name])
    p1 <-  p * k
  }
  cat("p,",p)
  cat(dfmax)
  cat("is.null(ref)", is.null(ref), "\r\n")
  cat(dfmax, "\r\n")
  return(W)
}


f1(X = X, ref=1)
