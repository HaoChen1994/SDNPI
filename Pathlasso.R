#########################################
# soft-thresholding function
soft.thred<-function(mu,omega)
{
  if(length(mu)==1)
  {
    return(max(abs(mu)-omega,0)*sign(mu))
  }else
  {
    return(pmax(abs(mu)-omega,rep(0,length(mu)))*sign(mu))
  }
}

# solution lemma for pathway lasso
pathlasso.sol<-function(lambda,omega,phi1,phi2,mu1,mu2)
{
  if(lambda==0)
  {
    a<-soft.thred(mu1,omega)/phi1
    b<-soft.thred(mu2,omega)/phi2
  }else
  {
    C1<-((phi2*mu1-lambda*mu2)>omega*(phi2-lambda))&((phi1*mu2-lambda*mu1)>omega*(phi1-lambda))
    C2<-((phi2*mu1+lambda*mu2)>omega*(phi2-lambda))&((phi1*mu2+lambda*mu1)<(-omega*(phi1-lambda)))
    C3<-((phi2*mu1+lambda*mu2)<(-omega*(phi2-lambda)))&((phi1*mu2+lambda*mu1)>omega*(phi1-lambda))
    C4<-((phi2*mu1-lambda*mu2)<(-omega*(phi2-lambda)))&((phi1*mu2-lambda*mu1)<(-omega*(phi1-lambda)))
    C5<-((abs(mu1)>omega)&((phi1*abs(mu2)-lambda*abs(mu1)<=omega*(phi1-lambda))))
    C6<-((abs(mu2)>omega)&((phi2*abs(mu1)-lambda*abs(mu2)<=omega*(phi2-lambda))))
    
    x1<-phi2*(mu1-omega)-lambda*(mu2-omega)
    x2<-phi2*(mu1-omega)+lambda*(mu2+omega)
    x3<-phi2*(mu1+omega)+lambda*(mu2-omega)
    x4<-phi2*(mu1+omega)-lambda*(mu2+omega)
    
    y1<-phi1*(mu2-omega)-lambda*(mu1-omega)
    y2<-phi1*(mu2+omega)+lambda*(mu1-omega)
    y3<-phi1*(mu2-omega)+lambda*(mu1+omega)
    y4<-phi1*(mu2+omega)-lambda*(mu1+omega)
    
    de<-phi1*phi2-lambda^2
    
    if(C1)
    {
      a<-x1/de
      b<-y1/de
    }else
      if(C2)
      {
        a<-x2/de
        b<-y2/de
      }else
        if(C3)
        {
          a<-x3/de
          b<-y3/de
        }else
          if(C4)
          {
            a<-x4/de
            b<-y4/de
          }else
            if(C5)
            {
              a<-(abs(mu1)-omega)*sign(mu1)/phi1
              b<-0
            }else
              if(C6)
              {
                a<-0
                b<-(abs(mu2)-omega)*sign(mu2)/phi2
              }else
              {
                a<-0
                b<-0
              }
  }
  
  return(list(a=a,b=b))
}
#########################################

#########################################
# two-block pathway lasso
library("glmnet")

# log-likelihood function
log.Lik<-function(X,M1,Y,beta,theta,delta)
{
  l1<-sum(diag(t(M1-X%*%beta)%*%(M1-X%*%beta)))*(-1/2)
  l2<-(t(Y-X*delta-M1%*%theta)%*%(Y-X*delta-M1%*%theta))[1,1]*(-1/2)
  
  re<-data.frame(M1=l1,Y=l2,sum=l1+l2)
  rownames(re)<-"logLik"
  
  return(re)
}

# objective function
obj.func<-function(X,M1,Y,beta,theta,delta,kappa1,kappa2,nu1,mu1)
{
  l1<-sum(diag(t(M1-X%*%beta)%*%(M1-X%*%beta)))
  l2<-(t(Y-X*delta-M1%*%theta)%*%(Y-X*delta-M1%*%theta))[1,1]
  
  l<-l1+l2
  
  p1<-kappa1*(sum(abs(beta)*abs(theta))+sum(nu1*(beta^2))+sum(nu1*(theta^2)))
  p2<-kappa2*abs(delta)
  p3<-mu1*(sum(abs(beta))+sum(abs(theta)))
  
  pen<-p1+p2+p3
  
  return(l/2+pen)
}

# data standardized
pathlasso.2b.std<-function(X,M1,Z,kappa1,kappa2,nu1=1,mu1=0,rho=1,max.itr=10000,tol=1e-6,
                           rho.increase=FALSE,trace=FALSE,beta0=NULL,theta0=NULL)
{
  n<-length(X)
  p1<-ncol(M1)
  
  if(trace)
  {
    beta.trace<-NULL
    theta.trace<-NULL
    
    td.beta.trace<-NULL
    td.theta.trace<-NULL

    delta.trace<-NULL
    
    tau1.trace<-NULL
    tau2.trace<-NULL
    
    objfunc<-NULL
  }
  
  if(rho.increase)
  {
    rho0<-rho
  }else
  {
    rho0<-0
  }
  
  # set initial values
  if(is.null(beta0))
  {
    td.beta0=beta0<-matrix(rep(0,p1),nrow=1)
  }else
  {
    td.beta0<-beta0
  }
  if(is.null(theta0))
  {
    td.theta0=theta0<-rep(0,p1)  
  }else
  {
    td.theta0<-theta0
  }
  
  delta0<-0
  
  tau1=tau2<-rep(0,p1)
  
  s<-0
  diff<-100
  
  pi.x  <- 1/(1+ exp( -X*delta0-M1%*%td.theta0))
  for ( i in 1:n )
  { 
    if (pi.x[i] <0.001) pi.x[i]=0.001
    else if (pi.x[i]>0.999) pi.x[i]=0.999
  }
  Y <- X*delta0+M1%*%td.theta0 + (Z-pi.x)/(1/4)
  
  Y <- Y-mean(Y)
  
  time<-system.time(
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      # update beta
      beta.new<-(t(X)%*%M1-t(tau1)+rho*td.beta0)/((t(X)%*%X)[1,1]+rho)
      # update theta
      theta.new<-solve(t(M1)%*%M1+rho*diag(rep(1,p1)))%*%(t(M1)%*%(Y-X*delta0)-tau2+rho*td.theta0)
    
      # update delta
      delta.new<-soft.thred((t(X)%*%(Y-M1%*%theta.new))[1,1],kappa2)/((t(X)%*%X)[1,1])
      
      # update tilde.beta and tilde.theta
      td.beta.new<-matrix(rep(NA,p1),nrow=1)
      td.theta.new<-rep(NA,p1)
      for(k in 1:p1)
      {
        sol.tmp<-pathlasso.sol(lambda=kappa1,omega=mu1,phi1=2*kappa1*nu1+rho,phi2=2*kappa1*nu1+rho,
                               mu1=tau1[k]+rho*beta.new[k],mu2=tau2[k]+rho*theta.new[k])
        td.beta.new[k]<-sol.tmp$a
        td.theta.new[k]<-sol.tmp$b
      }
      
      # update tau1, tau2, tau3, tau4
      tau1<-tau1+rho*c(beta.new-td.beta.new)
      tau2<-tau2+rho*c(theta.new-td.theta.new)
      
      if(trace)
      {
        beta.trace<-cbind(beta.trace,c(beta.new))
        theta.trace<-cbind(theta.trace,c(theta.new))
        
        td.beta.trace<-cbind(td.beta.trace,c(td.beta.new))
        td.theta.trace<-cbind(td.theta.trace,c(td.theta.new))
        
        delta.trace<-c(delta.trace,delta.new)
        
        tau1.trace<-cbind(tau1.trace,tau1)
        tau2.trace<-cbind(tau2.trace,tau2)
        
        objfunc<-c(objfunc,obj.func(X,M1,Y,td.beta.new,td.theta.new,delta.new,kappa1,kappa2,nu1,mu1))
      }
      
      diff.beta<-max(abs(beta.new-beta0))
      diff.theta<-max(abs(theta.new-theta0))
      diff.delta<-abs(delta.new-delta0)
      
      diff<-max(c(diff.beta,diff.theta,diff.delta))
      
      # print(data.frame(diff=diff,beta=diff.beta,theta=diff.theta,zeta=diff.zeta,pi=diff.pi,Lambda=diff.Lambda,delta=diff.delta,obj=objfunc[s]))
      
      beta0<-beta.new
      theta0<-theta.new
      delta0<-delta.new
      td.beta0<-td.beta.new
      td.theta0<-td.theta.new
      
      rho<-rho+rho0
      
      pi.x  <- 1/(1+ exp( -X*delta0-M1%*%theta0))
      for ( i in 1:n )
      { 
        if (pi.x[i] <0.001) pi.x[i]=0.001
        else if (pi.x[i]>0.999) pi.x[i]=0.999
      }
      Y <- X*delta0+M1%*%theta0 + (Z-pi.x)/(1/4)
      
      Y <- Y-mean(Y)
    })
  
  if(s>max.itr)
  {
    warning("Method does not converge!")
  }
  
  constraint1=constraint2<-matrix(NA,1,2)
  colnames(constraint1)=colnames(constraint2)<-c("beta=td.beta","theta=td.theta")
  constraint1[1,1]<-(max(abs(beta.new-td.beta.new))<tol)
  constraint1[1,2]<-(max(abs(theta.new-td.theta.new))<tol)
  constraint2[1,1]<-max(abs(beta.new-td.beta.new))
  constraint2[1,2]<-max(abs(theta.new-td.theta.new))
  constraint<-cbind(data.frame(t(constraint1)),data.frame(t(constraint2)))
  colnames(constraint)<-c("Satisfied","value")
  
  beta.est<-td.beta.new
  colnames(beta.est)<-paste0("M1.",1:p1)
  rownames(beta.est)<-"X"
  theta.est<-matrix(td.theta.new,ncol=1)
  rownames(theta.est)<-paste0("M1.",1:p1)
  colnames(theta.est)<-"Y"
  
  net.mat<-matrix(NA,1+p1,p1+1)
  rownames(net.mat)<-c("X",paste0("M1.",1:p1))
  colnames(net.mat)<-c(paste0("M1.",1:p1),"Y")
  net.mat[1,]<-c(beta.est,delta.new)
  net.mat[2:(p1+1),(p1+1):(p1+1)]<-cbind(theta.est)
  
  IE.M1.est<-c(beta.est)*theta.est
  names(IE.M1.est)<-paste0("M1.",1:p1)
  
  if(trace)
  {
    re<-list(beta=beta.est,theta=theta.est,
             delta=delta.new,para.mat=net.mat,IE.M1=IE.M1.est,
             logLik=log.Lik(X,M1,Y,td.beta.new,td.theta.new,delta.new),
             converge=(s<=max.itr),constraint=constraint,time=time,
             beta.trace=beta.trace,theta.trace=theta.trace,
             td.beta.trace=td.beta.trace,td.theta.trace=td.theta.trace,
             delta.trace=delta.trace,objfunc=objfunc,Y=Y)
  }else
  {
    re<-list(beta=beta.est,theta=theta.est,delta=delta.new,para.mat=net.mat,
             logLik=log.Lik(X,M1,Y,td.beta.new,td.theta.new,delta.new),
             converge=(s<=max.itr),constraint=constraint,time=time,Y=Y)
  }
  
  return(re)
}

pathlasso.2b<-function(X,M1,Z,kappa1,kappa2,nu1=1,mu1=0,rho=1,
                       standardize=TRUE,max.itr=10000,tol=1e-6,rho.increase=FALSE,trace=FALSE,
                       beta0=NULL,theta0=NULL)
{
  n<-length(X)
  p1<-ncol(M1)
  
  if(standardize)
  {
    # standardize data
    X.sd<-sd(X)
    M1.sd<-apply(M1,2,sd)
    M1.sd[M1.sd==0] <- 10^-3
    
    X.std<-scale(X,center=TRUE,scale=TRUE)
    M1.std<-scale(M1,center=TRUE,scale=TRUE)
    M1.std[is.nan(M1.std)] <- 10^-3
    
    re.std<-pathlasso.2b.std(X.std,M1.std,Z,kappa1=kappa1,kappa2=kappa2,
                             nu1=nu1,mu1=mu1,rho=rho,max.itr=max.itr,tol=tol,
                             rho.increase=rho.increase,trace=trace,
                             beta0=beta0,theta0=theta0)
    Y <- re.std$Y
    Y.sd<-sd(Y)
    beta.est<-re.std$beta*(M1.sd/X.sd)
    theta.est<-re.std$theta*(Y.sd/M1.sd)
    delta.est<-re.std$delta*(Y.sd/X.sd)
    
    
    net.mat<-matrix(NA,1+p1,p1+1)
    rownames(net.mat)<-c("X",paste0("M1.",1:p1))
    colnames(net.mat)<-c(paste0("M1.",1:p1),"Y")
    net.mat[1,]<-c(beta.est,delta.est)
    net.mat[2:(p1+1),(p1+1):(p1+1)]<-cbind(theta.est)
    
    IE.M1.est<-c(c(beta.est)*theta.est)
    names(IE.M1.est)<-paste0("M1.",1:p1)
    
    re<-list(beta=beta.est,theta=theta.est,
             delta=delta.est,para.mat=net.mat,IE.M1=IE.M1.est,
             logLik=log.Lik(X,M1,Y,beta.est,theta.est,delta.est),converge=re.std$converge,constraint=re.std$constraint,
             time=re.std$time,out.scaled=re.std,Y=Y)
  }else
  {
    re<-pathlasso.2b.std(X,M1,Z,kappa1=kappa1,kappa2=kappa2,
                         nu1=nu1,mu1=mu1,rho=rho,max.itr=max.itr,tol=tol,
                         rho.increase=rho.increase,trace=trace,
                         beta0=beta0,theta0=theta0,Y=Y)
  }
  
  return(re)
}

###############################################################################################################################

pathlasso.aic_opt <- function(X1, X2,  p , re, mu.prod, kappa1) {
  #----------------------------------------------------------------------#
  # Input: 
  #       X1, list, 
  #           control group data;
  #       X2, list, 
  #           case group data;
  #       p, numeric, 
  #           the spatial dimension of imaging(fMRI) data ;
  #       re, list
  #           the estimators under different tuning parameters.
  # Output:
  #       opt_aic, numeric
  #           the indices (row and column) of the optimal combination of tuning parameters that minimizes the AIC.
  #       cond_diff_matrix, matrix,
  #           conditional differential matrix
  #----------------------------------------------------------------------#
  
  n1 <- length(X1)
  n2 <- length(X2)
  p<-p
  
  # Initialize vectors
  bt <- vector("list", length = length(mu.prod))
  sum_bt <- vector("list", length = length(mu.prod))
  sum_aic <- vector("list", length = length(mu.prod))
  
  ## Compute AIC
  for (i in 1:length(mu.prod)) {
    bt[[i]] <- vector("list", length = length(kappa1))
    sum_bt[[i]] <- vector("list", length = length(kappa1))
    sum_aic[[i]] <- vector("list", length = length(kappa1))
    
    for (j in 1:length(kappa1)) {
      for (k in 1:(p * (p - 1) / 2)) {
        bt[[i]][[j]][k] <- re[[i]][[j]]$beta[k] * re[[i]][[j]]$theta[k]
        sum_bt[[i]][[j]] <- sum(bt[[i]][[j]] != 0)
        sum_aic[[i]][[j]] <- log(-2 * re[[i]][[j]]$logLik$sum) + sum_bt[[i]][[j]] + 
          (sum_bt[[i]][[j]]^2 + sum_bt[[i]][[j]]) / ((n1+n2) - sum_bt[[i]][[j]] - 1)
      }
    }
  }
  
  ##AIC select optimal,get the indices of the minimum AIC value
  final_sum <- matrix(nrow = length(mu.prod), ncol = length(kappa1))
  for (i in 1:length(mu.prod)) {
    for (j in 1:length(kappa1)) {
      final_sum[i, j] <- sum_aic[[i]][[j]]
    }
  }
  opt_aic <- which(final_sum == min(final_sum), arr.ind = TRUE)
  
  
  # Obtain the conditional differential matrix
  cond_diff_vec<-array()
  cond_diff_vec_rawnum <- bt[[opt_aic[1, 1]]][[opt_aic[1, 2]]]
  cond_diff_vec <- ifelse(cond_diff_vec_rawnum != 0, 1, 0)
  
  cond_diff_matrix <- matrix(0, nrow = p, ncol = p)
  cond_diff_matrix[upper.tri(cond_diff_matrix, diag = FALSE)] <- cond_diff_vec
  cond_diff_matrix[lower.tri(cond_diff_matrix, diag = FALSE)] <- t(cond_diff_matrix)[lower.tri(cond_diff_matrix, diag = FALSE)]
  
  final_re <- list(opt_aic,cond_diff_matrix)
  return(final_re)
}
