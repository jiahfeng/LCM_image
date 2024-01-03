#####################################################
# sim_EM.R
# Jiahui Feng
# update Jan. 2, 2024
# 
# simulation study of the latent-class model with imaging 
# data and survival outcomes using the EM algorithm
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

B <- aperm(array(desMat, c(40, 40, pX)), c(3, 1, 2)) #reorder, 3rd element is the first dimension etc
ind.outside <- matrix( rep( rep(NA, nrow(Z)), n), nrow = n )
mya <- array(dim = c(1, 40, 40))
for(i in 1){
  mya[i,,] <- matrix(ind.outside[i,], 40, 40)
}

g <- funData(list(c(1:40), c(1:40)), mya) 

M <- MFPCA:::calcBasisIntegrals(B, 2, g@argvals)


# initial values of the latent variable
G <- 10
init_z <- matrix(0, n, G)
init_z <- t(apply(init_z, 1, function(x){
  x[sample(G, 1)] = 1
  return(x)
}))
init_pi <- colMeans(init_z)


# tuning parameters
gamma_seq <- seq(sqrt(log(n)/n)*0.4, sqrt(log(n)/n)*0.6,length.out = 5)
lambda_seq <- seq(0.02, 0.06, length.out = 5)
para_combi <- expand.grid(gamma_seq, lambda_seq)

beta0 <- matrix(0, nrow = G, ncol = pX)
beta0[1, 1:2] <- c(-3, -2)
beta0[2, 1:2] <- c(1, 1)
lambda0g <- c(0.5, 0.25)



t <- 100 # number of replicates

res <- list()
for (j in 1:t) {
  set.seed(j)
  # generate the imaging data
  coef_image <- matrix(rnorm(n*pX), nrow = n, ncol = pX)
  image_data <- coef_image %*% M
  image_data <- apply(image_data, 2, scale)
  
  
  C <- c(rep(1, n/3), rep(2, 2*n/3)) # index of class
  T <- rep(0,n)
  D <- rep(NA,n)
  X <- image_data
  
  for (i in 1:n){
    u <- runif(1)
    if (C[i] == 2) t <- log(1 - lambda0g[C[i]] * log(u) / exp(sum(beta0[C[i],] * X[i,]))) / lambda0g[C[i]] else t <- -log(u) / (exp(sum(beta0[C[i],] * X[i,])) * lambda0g[C[i]])
    c <- runif(1,0,9)
    T[i] <- min(c,t)
    D[i] <- as.numeric(t<c) 
  }
  
  times <- T
  event <- D
  index <- C
  
  data_all <- data.frame(index, times, event, X)
  data_all <- data_all[order(data_all$times, -data_all$event), ]
  
  surv_data <- data.frame(time = data_all$times, event = data_all$event)
  
  pi_temp <- list()
  BIC_temp <- c()
  beta_temp <- list()
  for (k in 1:nrow(para_combi)) {
    para <- para_combi[k,]
    result <- mixture_cox_em(covs = data.matrix(data_all[,-c(1:3)]), data = surv_data, num_components = G, 
                             lambda = as.numeric(para[2]), gamma = as.numeric(para[1]), max_iter = 100, tol = 1e-3)
    BIC_temp <- c(BIC_temp, result$BIC)
    pi_temp[[k]] <- result$pi
    beta_temp[[k]] <- result$beta
  }
  
  index <- which.min(BIC_temp)
  
  order_pi <- order(pi_temp[[index]])
  pi_est <- c(rep(0, G))
  pi_est[1:length(pi_temp[[index]])] <- pi_temp[[index]][order_pi]
  beta_est <- beta_temp[[index]][,order_pi]
  
  all_result <- c(para_combi[index,], pi_est, beta_est)
  res[[j]] <- all_result
  
  # write.table(t(all_result), file = "mixture_EM.txt", row.names = FALSE, col.names = FALSE, append = TRUE)
}
