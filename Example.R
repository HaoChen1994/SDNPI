rm(list=ls())

library(DensParcorr)
library(mvtnorm)
library(MASS)
library(huge)
library(Matrix)
library(expm)
library(mnormt)
library(foreach)
library(doParallel)
library(igraph)
library(mvtnorm)
source("pathlasso.R")
source("EstMediator.R")

p <- 10
q <- 10
N <- 5
times <- 3
sparsity <- 0.8
XM_related <- floor((1-sparsity)*(p*(p-1)/2))
standardize<-TRUE
max.itr<-1000
tol<-1e-2
trace<-FALSE
rho<-1
rho.increase<-FALSE
nu1<-2
kappa1=kappa2<-10^c(seq(-5,-1,length.out=3))
mu.prod<-c(10)
opt_TPR=opt_TNR=TPR_hima=TNR_hima=TPR_rhima=TNR_rhima=opt_TPR_oracle=opt_TNR_oracle=array()
TPR_all=TNR_all=TPR_all_oracle=TNR_all_oracle <- matrix(nrow=times,ncol=(length(kappa1)*length(mu.prod)))

sigmaT = "autreg"

if(sigmaT=="autreg"){
  sigmaT1 <- 0.4^abs(outer(1:q,1:q,"-"))
  sigmaT2 <- 0.6^abs(outer(1:q,1:q,"-"))
} else if(sigmaT=="band"){
  sigmaT1 <- sigmaT2 <- 1/(abs(outer(1:q,1:q,"-"))+1)
  sigmaT1[abs(row(sigmaT1)-col(sigmaT1))>4] <- 0
  sigmaT2[abs(row(sigmaT2)-col(sigmaT2))>6] <- 0
} 

#group data generate
G1 <- huge.generator(
  n = 10,
  d = p,
  graph = "hub",
  g = 5,
  prob = NULL
)
Theta1 <- G1$theta
omega1.total = Theta1 * sample(c(-1,1), p * p, replace = TRUE) * runif(p * p, 0.3,0.5)
omega1.total[lower.tri(omega1.total, diag = FALSE)] <- 0
omega1.total <- omega1.total + t(omega1.total)
diag(omega1.total) = abs(min(eigen(omega1.total)$values)) + 0.5
sigma1.total = solve(omega1.total)

omega1.total = solve(sigma1.total)
omega1.total[abs(omega1.total) < 10 ^ -4] <- 0

omega2.total <- omega1.total
omega2.total[1:(p/5*2),1:(p/5*2)] <- -1*omega2.total[1:(p/5*2),1:(p/5*2)]

diag(omega2.total) <- diag(omega1.total)
sigma2.total <- solve(omega2.total)

delta <- omega1.total-omega2.total
delta <- as.matrix(delta)
delta_vec <- as.vector(delta[upper.tri(delta, diag = FALSE)])

omega1_N <- list()
omega2_N <- list()

for(i in 1:N){
  omega1_N[[i]] <- omega1.total + Theta1 * matrix(rnorm(p*p, mean = 0, sd = 0.2),nrow=p,ncol=p)
  omega2_N[[i]] <- omega2.total + Theta1 * matrix(rnorm(p*p, mean = 0, sd = 0.2),nrow=p,ncol=p)
} 

X1_w <- list()
X2_w <- list()

cores <- detectCores()
if(cores>50){
  cores <- 50
}
cl <- makeCluster(cores)
registerDoParallel(cl)

X1_w <- foreach(i = 1:N, .errorhandling = 'pass', .packages = c("Matrix", "mnormt","expm")) %dopar% {
  #control group
  omega1 <- as.matrix(omega1_N[[i]])
  omega1[lower.tri(omega1, diag = FALSE)] <- 0
  omega1 <- omega1 + t(omega1)
  diag(omega1) = abs(min(eigen(omega1)$values)) + 1
  sigma1 = solve(omega1)
  # 
  omega1 = solve(sigma1)
  omega1[abs(omega1) < 10 ^ -4] <- 0
  # 
  Z = rmnorm(n=p, mean=0, sigmaT1)
  SQS = sqrtm(sigma1)
  X1_w = SQS%*%Z
  rm(Z)
  gc()
  rm(SQS)
  gc()
  return(X1_w)
}


X2_w <- foreach(i = 1:N, .errorhandling = 'pass', .packages = c("Matrix", "mnormt","expm")) %dopar% {
  #case group
  omega2 <- omega2_N[[i]]
  omega2[lower.tri(omega2, diag = FALSE)] <- 0
  omega2 <- omega2 + t(omega2)
  diag(omega2) = abs(min(eigen(omega2)$values)) + 1
  sigma2 = solve(omega2)
  
  omega2 = solve(sigma2)
  omega2[abs(omega2) < 10 ^ -4] <- 0
  
  Z = rmnorm(n=p, mean=0, sigmaT2)
  SQS = sqrtm(sigma2)
  X2_w = SQS%*%Z
  rm(Z)
  gc()
  rm(SQS)
  gc()
  return(X2_w)
}

parallel::stopCluster(cl)

Mediators<-Estimate_Mediators(X1_w,X2_w,p,cores,0.5)
X1_vec<-Mediators$X1_vec
X2_vec<-Mediators$X2_vec

M1_1 <- X1_vec
M1_2 <- X2_vec
M1 <- rbind(M1_1,M1_2)
alpha <- c(rep(10,XM_related),rep(0,(p*(p-1)/2-XM_related)))
# X
eps.M1<-rnorm(2*N,mean=0,sd=0.1)
X <- M1%*%alpha+eps.M1
# Z
Z <- array(c(rep(0,N),rep(1,N)),dim = c(2*N,1))

# Obtain the results for all hyper-parameters
cores <- detectCores()
if(cores>50){
  cores <- 50
}
cl <- makeCluster(cores)
registerDoParallel(cl)
re<-vector("list",length=length(mu.prod))
for(ss in 1:length(mu.prod))
{
  re[[ss]]<-vector("list",length=length(kappa1))
  mu1<-mu.prod[ss]*kappa1
  re[[ss]] <- foreach(i = 1:length(kappa1), .errorhandling = 'pass') %dopar% {
    re[[ss]]<-pathlasso.2b(X,M1,Z,kappa1=kappa1[i],kappa2=kappa2[i],nu1=nu1,mu1=mu1[i],rho=rho,standardize=standardize,
                           max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=FALSE,beta0=NULL,theta0=NULL)
    return(re[[ss]])
  }
}
parallel::stopCluster(cl)

optimal_result<-pathlasso.aic_opt(X1_w, X2_w, p=p, re, mu.prod=mu.prod, kappa1=kappa1)

cond_diff_matrix<-optimal_result[[2]]
