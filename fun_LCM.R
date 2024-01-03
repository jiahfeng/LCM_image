#####################################################
# fun_LCM.R
# Jiahui Feng
# update Jan. 2, 2024
# 
# support functions of the EM algorithm for a latent-calss model
# 
# Functions:
#
# bh_smooth            # kernel smooth of the baseline hazard function
# phi                  # Gaussian function as the kernel, changeable
# Base_Surv            # estimation function for the baseline hazard
# mixture_cox_em       # EM alogrithm for the latent-class model
# soft_threshold       # soft threshold function for L1 regularization
# multi_logistic_lasso # coordinate descent for multinomial logistic regression with a Lasso penalty
# mixture_cox_gener_em # generalized EM alogrithm for the latent-class model
# 
#####################################################


bernstein <- function(Y, V.est, Tr.est, d.est, r, Z, lambda){
  
  # Y: imaging data
  # V.est: vertical points of triangles
  # Tr.est: triangles
  # d.est: degree of piecewise polynomials
  # r: smoothness parameter
  # Z: expanded grid points
  # lambda: tuning parameter
  
  n <- nrow(Y)
  Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
  B <- Bfull.est$B
  ind.inside <- Bfull.est$Ind.inside
  Q2 <- Bfull.est$Q2
  K <- Bfull.est$K
  Y <- matrix(Y[,ind.inside],nrow=n)
  lambda <- as.matrix(lambda)
  t.area <- Bfull.est$tria.all 
  
  this.call <- match.call()
  n <- nrow(Y)
  npix <- ncol(Y)
  J <- ncol(Q2)
  
  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,t(Y))
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  
  gcv_all <- sapply(lambda,FUN=function(Lam){  
    Dlam <- Lam*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    theta <- crossprod(t(lhs.inv),rhs)
    gamma <- crossprod(t(Q2),theta)
    Yhat <- crossprod(t(W),theta)
    res <- t(Y)-Yhat
    sse <- apply(res^2,2,sum)
    if(!flag){
      df <- sum(1/(1+Cval*Lam))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df <- sum(diag(Hmtx))
    }
    gcv <- npix*sse/(npix-df)^2
  })
  gcv_all <- matrix(gcv_all,nrow=n)
  lam.ind <- apply(gcv_all,1,which.min)
  lambdac <- lambda[lam.ind]
  
  theta <- c()
  gamma <- c()
  Yhat <- c()
  df <- c()
  for (i in 1:n){
    lamc.tmp <- lambdac[i]
    Dlam <- lamc.tmp*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    rhs.tmp <- as.matrix(rhs[,i],ncol=1)
    theta.tmp <- crossprod(t(lhs.inv),rhs.tmp)
    theta <- cbind(theta,theta.tmp)
    gamma.tmp <- crossprod(t(Q2),theta.tmp) 
    gamma <- cbind(gamma,gamma.tmp)
    Yhat.tmp <- crossprod(t(W),theta.tmp)
    Yhat <- cbind(Yhat,Yhat.tmp)
    if(!flag){
      df.tmp <- sum(1/(1+Cval*lamc.tmp))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df.tmp <- sum(diag(Hmtx))
    }
    df <- c(df,df.tmp)
  }
  
  return(list(B, gamma, lambdac, t.area, ind.inside))
}


fpca_img <- function(type = "bernstein", Y, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr){
  
  # type: basis function type
  # Y: imaging data
  # lambda: tuning parameter
  # npc: number of PCs
  # V.est: vertical points of triangles
  # Tr.est: triangles
  # d.est: degree of piecewise polynomials
  # r: smoothness parameter
  # Z: expanded grid points
  # ncr: dimension of the input images
  
  if(type == 'bernstein'){
    est <- bernstein(Y= Y, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                     r = r, Z = Z, lambda = lambda)
    basis <- est[[1]];scores <- t(est[[2]]);t.area <- est[[4]];ind.inside <- est[[5]]
  }
  
  org.idx <- c(1:ncr*ncr)
  desMat <- matrix(0, nrow = ncr*ncr, ncol = ncol(scores))
  for(zz in 1:ncol(desMat)){
    desMat[ind.inside,zz] = basis[,zz]
  }
  
  
  B <- aperm(array(desMat, c(ncr, ncr,ncol(scores))), c(3, 1, 2)) 
  mya <- array(dim = c(1, ncr, ncr))
  for(i in 1){
    mya[i,,] <- matrix(Y[i,], ncr, ncr)
  }
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W <- MFPCA:::calcBasisIntegrals(B, 2,g@argvals) 
  
  S <- t(scores)
  sqrM <- function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  
  U <-1/nrow(Y)*W%*%S%*%t(S)%*%W; G <- W 
  halfG_inv  <- sqrM(G)
  tM <- t(halfG_inv)%*%U%*%halfG_inv
  eigen_res <- eigen(tM) 
  
  fd_list <- lapply(1:npc, function(ipc){
    coef_pc<- halfG_inv%*%as.matrix(Real(eigen_res$vectors[,ipc])) 
  })
  
  sup_basis <- NULL
  for(k in 1:npc){
    kth <- matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis <- cbind(sup_basis, t(kth))
  }
  
  
  return(list(sup_basis, eigen_res$values[1:npc], basis, scores, t.area))
}



bh_smooth <- function(time, bw, n, h0){
  
  # time: observed time points
  # bw: bandwidth
  # n: sample size
  # h0: estimates of baseline hazard function
  
  h_vec <- rep(0,n)
  for(i in 1:n){
    h_vec[i] <- 1 / bw * sum(phi((time - time[i])/bw)*h0)
  }
  return(h_vec)
}

phi <- function(x){
  exp(-0.5*x^2)/sqrt(2*pi)
}


#### Baseline Survival 
Base_Surv <- function(covs, data, beta0, z0, num_components){
  
  # covs: Matrix of observed covariates for the survival function
  # data: Data frame containing survival times and event
  # beta0: Matrix of regression coefficients
  # z0: Latent variable matrix, indicating the weight of each observation in each class
  # num_components: Number of mixture classes
  
  n <- nrow(data)
  beta <- beta0
  z <- z0
  obs <- which(data$event == 1)
  time <- data$time
  
  hazard_value <- matrix(NA, n, num_components)
  
  event <- data$event
  for (k in 1:num_components) {
    z_k <- z[,k]
    beta_k <- beta[,k]
    h <- c(rep(0, n)) #baseline hazard function
    for(i in 1:n){
      zexpx <- z_k * exp(covs %*% beta_k)
      riskset <- which(time >= time[i])
      h[i] <- z_k[i] * event[i] / (sum(zexpx[riskset]))
    }
    h[which(is.na(h))] <- 0
    hazard_value[,k] <- h
  }
  
  H0 <- sapply(c(1:num_components), function(x) cumsum(hazard_value[,x]))
  
  smoothed_hazard_value = sapply(c(1:num_components), function(x) bh_smooth(time, bw = IQR(time)*n^(-1/5), n, hazard_value[,x]))
  
  list(base_hazard = smoothed_hazard_value, cum_hazard = H0)
}


# Mixture Cox Proportional Hazards model using EM algorithm
mixture_cox_em <- function(covs, data, num_components, bandwidth, lambda, gamma, max_iter = 100, tol = 1e-6) {
  
  # covs: Matrix of covariates for the survival data
  # data: Data frame with survival time and event indicator
  # num_components: Number of classes in the model
  # bandwidth: Bandwidth parameter for kernel smoothing
  # lambda: Regularization parameter for the Cox model
  # gamma: Shrinkage parameter for mixture weights
  # max_iter: Maximum number of iterations for the EM algorithm
  # tol: Tolerance level for convergence in the EM algorithm
  
  n <- nrow(covs)
  p <- ncol(covs)
  time <- data$time
  event <- data$event
  
  # Initialize parameters
  z0 <- init_z
  pi0 <- init_pi
  
  
  # Initial estimate for regression coefficients
  beta0 <- matrix(0, p, num_components)
  for (k in 1:num_components) {
    fit <- glmnet(x = covs, y = Surv(time, event), family = "cox", weights = z0[,k], lambda = lambda)
    beta0[,k] <- as.vector(coef(fit))
  }
  
  for (iter in 1:max_iter) {
    ## E-step: Estimate the expected value of the latent variables ##
    surv0 <- Base_Surv(covs = covs, data = data, beta0 = beta0, z0 = z0, num_components = num_components)
    base_hazard <- surv0$base_hazard
    cum_hazard <- surv0$cum_hazard
    
    
    # Compute likelihood for each class
    lik <- matrix(0, n, num_components)
    for (k in 1:num_components) {
      expx_k <- exp(covs %*% beta0[,k])
      lik[,k] <- (base_hazard[,k] * expx_k)^event * exp(-expx_k * cum_hazard[,k])
    }
    
    z_update <- matrix(0, n, num_components)
    temp <- lik %*% diag(pi0, num_components)
    z_update <- do.call(rbind, lapply(1:n, function(i) temp[i,] / sum(temp[i,])))
    z_update[which(is.na(z_update))] <- 0
    
    
    ## M-step: Maximize the expected complete-data log-likelihood ##
    pi_update <- sapply(colMeans(z_update), function(x) max(0, (x - gamma) / (1 - num_components * gamma)))
    ind <- (pi_update>0)
    num_components <- sum(ind)
    pi_update <- pi_update[ind]
    pi_update <- pi_update / sum(pi_update)
    
    lik <- as.matrix(lik[,ind])
    z_update <- matrix(0, n, num_components)
    temp <- lik %*% diag(pi_update, num_components)
    z_update <- do.call(rbind, lapply(1:n, function(i) temp[i,] / sum(temp[i,])))
    z_update[which(is.na(z_update))] <- 0
    z_update[which(z_update < 1e-5)] <- 0
    
    beta_update <- matrix(0, p, num_components)
    for (k in 1:num_components) {
      fit <- glmnet(x = covs, y = Surv(time, event), family = "cox", weights = z_update[,k], lambda = lambda)
      beta_update[,k] <- as.vector(coef(fit))
    }
    
    # Check for convergence
    if (max(abs(pi_update - pi0[ind])) < tol) {
      break
    }
    
    # Update parameters
    pi0 <- pi_update
    beta0 <- beta_update
    z0 <- z_update
  }
  
  m <- length(pi_update)
  temp <- lik %*% diag(pi_update, m)
  BIC <- -2 * sum(log(rowSums(temp))) + (m - 1 + sum(beta_update!=0))*log(n)
  
  
  return(list(BIC = BIC, pi = pi_update, beta = beta_update, z = z_update))
}



# Soft threshold function for L1 regularization
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}


# Coordinate Descent for Multinomial Logistic Regression with Lasso Penalty
multi_logistic_lasso <- function(covs, z0, alpha0, lambda = 0, num_components, max_iter = 100, tol = 1e-3) {
  
  # covs: Matrix of covariates
  # z0: Latent variable matrix, indicating the weight of each observation in each class
  # alpha0: Initial values of the coefficients 
  # data: Data frame with survival time and event indicator
  # num_components: Number of classes in the model
  # lambda: Regularization parameter for the multinomial logistic regression model
  # max_iter: Maximum number of iterations
  # tol: Tolerance level for convergence
  
  
  n <- nrow(covs)
  p <- ncol(covs)
  if (lambda == 0) lambda = c(rep(0, p))
  alpha <- alpha0
  z <- z0
  iter <- 0
  converged <- FALSE
  prev_probG <- matrix(1, nrow = num_components, ncol = n)
  
  
  while (iter < max_iter){
    iter <- iter + 1
    for (g in 1: num_components){
      alpha_g <- alpha[g,]
      eWa <- exp(covs %*% t(alpha))
      probG <- apply(eWa, 1, function(x) x / sum(x))
      w_g <- probG[g,] * (1 - probG[g,])
      w_g[which(1 - probG[g,] < 1e-5)] <- 1e-5
      w_g[which(probG[g,] < 1e-5)] <- 1e-5
      y_g <- covs %*% alpha_g + (z[,g] - probG[g,]) / w_g
      
      iter_inner <- 0
      while (iter_inner < max_iter) {
        iter_inner <- iter_inner + 1
        prev_alpha_g <- alpha_g
        for (j in 1:p) {
          covs_j <- covs[, j]
          covs_without_j <- covs[, -j]
          gradient <- mean(w_g * covs_j * (y_g - covs_without_j %*% alpha_g[-j]))
          alpha_g[j] <- soft_threshold(gradient, lambda[j]) / (mean(w_g * covs_j^2))
        }
        if (max(abs(alpha_g - prev_alpha_g)) < tol) {
          alpha[g,] <- alpha_g
          break
        }
      }
    }
    if (max(abs(probG - prev_probG)) < tol){
      converged <- TRUE
      break
    }
    prev_probG <- probG
  }
  
  
  return(alpha)
}



mixture_cox_gener_em <- function(covs_logistic, covs_cox, data, num_components, lambda, bandwidth, max_iter = 100, tol = 1e-3){
  
  # covs_logistic: Matrix of covariates for logistic regression part of the model
  # covs_cox: Matrix of covariates for Cox proportional hazards part of the model
  # data: Data frame with survival time and event indicator
  # num_components: Number of classes
  # lambda: Regularization parameter for the Cox model.
  # bandwidth: Bandwidth parameter for kernel smoothing in the Cox model
  # max_iter: Maximum number of iterations for the EM algorithm
  # tol: Tolerance level for convergence in the EM algorithm
  
  n <- nrow(covs_logistic)
  time <- data$time
  event <- data$event
  
  # Initialize parameters
  z0 <- init_z
  
  
  # Initial fitting of Cox models for each component
  beta0 <- matrix(0, pX, num_components)
  for (k in 1:num_components) {
    # cv_fit <- cv.glmnet(x = covs_cox, y = Surv(time, event), family = "cox", weights = z0[,k])
    fit <- glmnet(x = covs_cox, y = Surv(time, event), family = "cox", weights = z0[,k], lambda = lambda)
    beta0[,k] <- as.vector(coef(fit))
  }
  
  # Initial logistic regression fitting
  alpha0 <- multi_logistic_lasso(covs = covs_logistic, z0 = z0, alpha0 = matrix(0,nrow=num_components,ncol=pW+1), 
                                 num_components = num_components)
  
  for (iter in 1:max_iter) {
    ## E-step: Estimate the expected value of the latent variables ##
    surv0 <- Base_Surv(covs = covs_cox, data = data, beta0 = beta0, z0 = z0, num_components = num_components)
    base_hazard <- surv0$base_hazard
    cum_hazard <- surv0$cum_hazard
    
    lik <- matrix(0, n, num_components)
    for (k in 1:num_components) {
      expx_k <- exp(covs_cox %*% beta0[,k])
      lik[,k] <- (base_hazard[,k] * expx_k)^event * exp(-expx_k * cum_hazard[,k])
    }
    
    z_update <- matrix(0, n, num_components)
    temp <- lik * exp(covs_logistic %*% t(alpha0))
    z_update <- temp / rowSums(temp)
    z_update[which(is.na(z_update))] <- 0
    z_update[which(z_update < 1e-5)] <- 0
    
    
    ## M-step: Maximize the expected complete-data log-likelihood ##
    alpha_update <- multi_logistic_lasso(covs = covs_logistic, z0 = z_update, alpha0 = alpha0, num_components = num_components)
    
    beta_update <- matrix(0, pX, num_components)
    for (k in 1:num_components) {
      fit <- glmnet(x = covs_cox, y = Surv(time, event), family = "cox", weights = z_update[,k], lambda = lambda)
      beta_update[,k] <- as.vector(coef(fit))
    }
    
    # Check for convergence
    if (max(abs(alpha_update - alpha0)) < tol &
        max(abs(beta_update - beta0)) < tol) {
      break
    }
    
    # Update parameters
    beta0 <- beta_update
    z0 <- z_update
    alpha0 <- alpha_update
  }
  
  eWa <- exp(covs_logistic %*% t(alpha_update))
  probG <- t(apply(eWa, 1, function(x) x/sum(x)))
  s <- sum(alpha_update!=0) + sum(beta_update!=0)
  BIC <- -2 * sum(log(rowSums(probG * lik))) + s * log(n)
  
  return(list(BIC = BIC, alpha = alpha_update, beta = beta_update))
}








