#####################################################
# sim_gener_EM.R
# Jiahui Feng
# update Jan. 2, 2024
# 
# simulation study of the latent-class model with imaging 
# data and survival outcomes using the generalized EM algorithm
#
#####################################################



library(survival)
library(glmnet)
library(MASS)
library(funData)
library(BPST)
library(Triangulation)
library(MFPCA)

source("fun_LCM.R")

n <- 1200
xm <- seq(0, 1, length = 40)
yn <- seq(0, 1, length = 40)
xy_grid <- expand.grid(xm, yn)
xx <- xy_grid$Var1
yy <- xy_grid$Var2
Z <- cbind(xx, yy)


load("Brain.V.rda") # vertices for triangles for simulation setting
load("Brain.Tr.rda") # triangles
d.est <- 2 # degree of piecewise polynomials
r <- 1 # smoothness parameter

Bfull.est <- basis(V,Tr,d.est,r,Z)
B.save <- t(as.matrix(Bfull.est$B)) # DIM: npts * ntriangle*((d+1)(d+2)/2)
ind.inside <- Bfull.est$Ind.inside
Q2 <- Bfull.est$Q2
pX <- nrow(B.save)

# -- inner product --- #
org.idx <- c(1:40*40)
desMat <- matrix(0, nrow = 40*40, ncol = pX)

for(zz in 1:ncol(desMat)){
  desMat[ind.inside,zz] <- B.save[zz,]
}

B <- aperm(array(desMat, c(40, 40, pX)), c(3, 1, 2)) #reorder, 3rd element is the first dimension etc.
ind.outside <- matrix( rep( rep(NA, nrow(Z)), n), nrow = n )
mya <- array(dim = c(1, 40, 40))
for(i in 1){
  mya[i,,] <- matrix(ind.outside[i,], 40, 40)
}

g <- funData(list(c(1:40), c(1:40)), mya) 

M <- MFPCA:::calcBasisIntegrals(B, 2, g@argvals)


pW <- 2
G <- 2
alpha0 <- matrix(0,nrow=G,ncol=pW+1)
alpha0[1,] <- c(-1.1,0.5,-0.5)
# alpha0[2,] <- c(-0.25,0.5,0.5)
if (G==2) lambda0g <- c(0.5,0.25)
if (G==3) lambda0g <- c(0.5,0.25,1)


beta0 <- matrix(0, nrow = G, ncol = pX)
beta0[1, 1:2] <- c(-3, -2)
beta0[2, 1:2] <- c(1, 1)

lambda_seq <- seq(0.03, 0.06, length.out = 5)
Ghat <- c(2:4)
para_combi <- expand.grid(Ghat, lambda_seq)




t <- 100 # number of replicates

res <- list()
for (j in 1:t) {
  set.seed(j)
  
  # generate the imaging data
  coef_image <- matrix(rnorm(n*pX), nrow = n, ncol = pX)
  image_data <- coef_image %*% M
  image_data <- apply(image_data, 2, scale)
  
  ## FPC
  Y <- matrix(NA, n, 40*40)
  Y[,ind.inside] <- coef_image %*% B.save
  fpc_fun <- fpca_img(type = "bernstein", Y = Y, lambda = c(0, 1, 10, 10^2, 10^3, 10^6), npc = 30, 
                      V.est = V, Tr.est = Tr, d.est = d.est, r = r, Z = Z, ncr = 40)
  
  o_tfs_fpc = as.matrix(fpc_fun[[1]]) ## FPC
  
  b.basis_fpc = as.matrix(fpc_fun[[3]]); b.scores_fpc = as.matrix(fpc_fun[[4]])
  
  score_fpc = t(o_tfs_fpc) %*% (b.basis_fpc%*%t(b.scores_fpc)) ## score

  X <- image_data
  W <- cbind(1, t(score_fpc[1:2,]))
  eWa <- matrix(nrow=n,ncol=G)
  probG0 <- matrix(nrow=n,ncol=G)
  C <- rep(NA,n)
  T <- rep(0,n)
  D <- rep(NA,n)
  
  for (i in 1:n){
    set.seed(i)
    u <- runif(1)
    c <- runif(1,0,8)
    for (g in 1:G) eWa[i,g] <- exp(W[i,]%*%as.vector(alpha0[g,]))
    probG0[i,] <- eWa[i,] / sum(eWa[i,])
    C[i] <- which(rmultinom(1, 1, probG0[i,])==1)
    if (C[i] == 2) t <- log(1 - lambda0g[C[i]] * log(u) / exp(sum(beta0[C[i],] * X[i,]))) / lambda0g[C[i]] else t <- -log(u) / (exp(sum(beta0[C[i],] * X[i,])) * lambda0g[C[i]])
    T[i] <- min(c,t)
    D[i] <- as.numeric(t<c) 
  }
  
  times <- T
  event <- D
  index <- C
  
  data_all <- data.frame(index,times, event, X, W)
  data_all <- data_all[order(data_all$times, -data_all$event), ]
  
  surv_data <- data.frame(time = data_all$times, event = data_all$event)
  
  alpha_temp <- list()
  beta_temp <- list()
  BIC <- c()
  for (k in 1:nrow(para_combi)){
    set.seed(k)
    para <- para_combi[k,]
    G_hat <- as.numeric(para[1])
    init_z <- matrix(0, n, G_hat)
    init_z <- t(apply(init_z, 1, function(x){
      x[sample(G_hat, 1)] = 1
      return(x)
    }))
    
    result <- mixture_cox_gener_em(covs_logistic = data.matrix(data_all[,-c(1:(3+pX))]), covs_cox = data.matrix(data_all[,c(4:(3+pX))]), 
                             data = surv_data, num_components = G_hat, lambda = as.numeric(para[2]), max_iter = 100, tol = 1e-3)
    BIC <- c(BIC, result$BIC)
    alpha_temp[[k]] <- result$alpha
    beta_temp[[k]] <- result$beta
  }
  
  index <- which.min(BIC)
  alpha_est <- alpha_temp[[index]]
  beta_est <- beta_temp[[index]]
  num_G <- nrow(alpha_est)
  
  all_result <- c(num_G, t(alpha_est), beta_est)
  
  res[[j]] <- all_result
}

