ERROR_fun <- function(beta_fit, beta_true,
                      basis_fit, basis_true, sseq,
                      m, p, Nzero_group, tol = 0) {

  error.list <- list()

  pos_group <- 1:Nzero_group
  neg_group <- (1:p)[-pos_group]
  df_fit <- (length(beta_fit) - 1 - m) / p
  df_true <- (length(beta_true) -1 -m) / p

  Non.zero_select <- Nzero(X = beta_fit, p = p, k = df_fit)

  error.list$class_error <- Class_error(pos_select = Non.zero_select,
                                        pos_group = pos_group, neg_group = neg_group)


  error.list$coef_error <- Coef_error(coef_true = beta_true, coef_esti = beta_fit,
                                      p = p, k_fit = df_fit, k_true = df_true)


  error.list$curve_error <- Curve_error(coef_true = beta_true, coef_esti = beta_fit, p = p,
                                        k_fit = df_fit, k_true = df_true,
                                        basis_true = basis_true, basis_fit = basis_fit, sseq = sseq)

  error.list$Non.zero <- Non.zero_select
  return(error.list)
}

##classification error
Class_error <- function(pos_select, pos_group, neg_group) {
  pos <- length(pos_group)
  neg <- length(neg_group)
  p <- pos + neg
  FP_select <- pos_select[which( pos_select %in% neg_group )]
  FN_select <- pos_group[which( !(pos_group %in% pos_select) )]
  TP_select <- pos_select[which( pos_select %in% pos_group )]
  TN_select <- neg_group[which( !(neg_group %in% pos_select) )]

  FP <- length(FP_select)
  FN <- length(FN_select)
  TP <- length(TP_select)
  TN <- length(TN_select)

  FPR <- FP / neg # false positive rate
  FNR <- FN / pos # false negative rate
  Sensitivity <- TP / pos # true positive rate
  Specificity <- TN / neg # true negative rate
  PPV <- TP / (TP + FP)   # precision/positive predictive value
  NPV <- TN/ (TN + FN)    # negative predictive value
  FDR <- 1 - PPV          # false discovery rate
  ACC <- (TP + TN) / p
  F1_score <- (2 * TP) / (2 * TP + FP + FN)
  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) # Matthews correlation coefficient
  Informedness <- Sensitivity + Specificity - 1 # Youden's J statistic
  Markedness <- PPV + NPV - 1

  error <- c(FP = FP, FN = FN, TP = TP, TN = TN,
             FPR = FPR, FNR = FNR, Sensi = Sensitivity, Speci = Specificity,
             PPV = PPV, NPV = NPV, FDR = FDR, ACC = ACC, F1_score = F1_score,
             MCC = MCC, Informedness = Informedness, Markedness = Markedness)

  return(error)
}

##Coefficients estimation Norm error
Coef_error <- function(coef_true, coef_esti, p, k_fit, k_true) {
  coef_esti.comp <- vet(X = coef_esti, p = p, k = k_fit)$C
  coef_esti.cont <- vet(X = coef_esti, p = p, k = k_fit)$b
  coef_true.comp <- vet(X = coef_true, p = p, k = k_true)$C
  coef_true.cont <- vet(X = coef_true, p = p, k = k_true)$b
  #cat("1")

  error <- vector()
  error <- c(L2.cont = sqrt(sum( (coef_esti.cont - coef_true.cont)^2 )) )
  if(k_fit == k_true) {
    L2_diff.beta <- sqrt(sum( (coef_true - coef_esti)^2 ))
    L2_diff.comp <- apply(coef_esti.comp - coef_true.comp, 1,
                          function(x) sqrt(sum(x^2)) )
    L2_diff_inf.comp <- max(L2_diff.comp)
    L2_diff_L1.comp <- sum(L2_diff.comp)
    L2_diff.comp <- sqrt(sum( (as.vector(coef_esti.comp - coef_true.comp))^2 ))
    error <- c(error, L2.beta = L2_diff.beta, L2.comp = L2_diff.comp,
               L2_inf.comp = L2_diff_inf.comp, L2_L1.comp = L2_diff_L1.comp)
  }
  #cat('2')
  return(error)
}


Curve_error <- function(coef_true, coef_esti, p,
                        k_fit, k_true, basis_true, basis_fit, sseq) {
  coef_esti.comp <- vet(X = coef_esti, p = p, k = k_fit)$C
  coef_true.comp <- vet(X = coef_true, p = p, k = k_true)$C

  curve_true <- coef_true.comp %*% t(basis_true)
  curve_esti <- coef_esti.comp %*% t(basis_fit)


  curve_diff <- abs(curve_esti - curve_true) ## p by length(sseq)

  ns <- length(sseq)
  time_diff <- sseq[2] - sseq[1]
  extra_sum <- rowSums(curve_diff[, c(1, ns)]^2) * time_diff / 2## crossprod(curve_diff[, c(1, ns)]) * time_diff/2
  add_sum <- apply(curve_diff, 1, function(x) sum(x^2)) * time_diff
  ITG <- add_sum - extra_sum
  L2 <- sqrt(ITG)
  L2_L1 <- sum(L2)
  L2_inf <- max(L2)


  extra_sum <- rowSums(curve_diff[, c(1, ns)]) * time_diff / 2## crossprod(curve_diff[, c(1, ns)]) * time_diff/2
  add_sum <- rowSums(curve_diff) * time_diff
  ITG <- add_sum - extra_sum
  L1 <- ITG
  L1_L1 <- sum(L1)
  L1_inf <- max(L1)

  L1_each.max <- max(curve_diff)
  error <- c(L1_each.max = L1_each.max,
             L1_inf = L1_inf, L1_L1 = L1_L1,
             L2_inf = L2_inf, L2_L1 = L2_L1 )
  return(error)
}
