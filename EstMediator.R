Estimate_Mediators <- function(X1, X2, p,cores = NULL, dens.level = 0.5) {
  #----------------------------------------------------------------------#
  # Input: 
  #       X1, list, 
  #           control group data;
  #       X2, list, 
  #           case group data;
  #       p, numeric, 
  #           the spatial dimension of imaging(fMRI) data ;
  #       cores, numeric, 
  #           the cores of Parallel;
  #       dens.level, numeric,
  #           the parameters in the DensParcorr function.
  # Output:
  #       X1_vec, n1 * (p * (p - 1) / 2)matrix, X2_vec, n2 * (p * (p - 1) / 2) matrix
  #           the estimated mediator matrix.
  #----------------------------------------------------------------------#

  
  n1 <- length(X1)
  n2 <- length(X2)
  p <- p
  
  ## Initialize 
  X1_omega <- list()
  X2_omega <- list()
  W1 <- list()
  W2 <- list()
  X1_vec <- matrix(nrow = n1, ncol = p * (p - 1) / 2)
  W1_vec <- array()
  X2_vec <- matrix(nrow = n2, ncol = p * (p - 1) / 2)
  W2_vec <- array()
  
  
  ##Parallel: Speed up the execution
  if(is.null(cores)){
    cores <- detectCores()
  }else{
    cores <- cores
  }

  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  ##Compute the precision matrix: X1_omega
  X1_omega <- foreach(i = 1:n1, .errorhandling = 'pass', .packages = "DensParcorr") %dopar% {
    X1.clime <- DensParcorr(X1[[i]], dens.level = dens.level, select = TRUE)
    X1_omega <- X1.clime$selected.precision
    rm(X1.clime)
    gc()
    return(X1_omega)
  }
  print("X1_omega Done!")
  
  ##Compute the precision matrix: X2_omega
  X2_omega <- foreach(i = 1:n2, .errorhandling = 'pass', .packages = "DensParcorr") %dopar% {
    X2.clime <- DensParcorr(X2[[i]], dens.level = dens.level, select = TRUE)
    X2_omega <- X2.clime$selected.precision
    rm(X2.clime)
    gc()
    return(X2_omega)
  }
  print("X2_omega Done!")
  
  parallel::stopCluster(cl)
  
  
  ##Fisher transformation of the partial correlation
  mut_in_f <- function(x) {
    return(1/2 * log((1 + x) / (1 - x)))
  }
  
  W1 <- lapply(lapply(X1_omega, cov2cor), mut_in_f)
  W2 <- lapply(lapply(X2_omega, cov2cor), mut_in_f)
  
  ##Vectorization and obtain the estimated mediator matrix
  for (i in 1:n1) {
    W1_vec <- as.vector(W1[[i]][upper.tri(W1[[i]], diag = FALSE)])
    X1_vec[i, ] <- c(W1_vec)
    
  }
  
  for (i in 1:n2) {
    W2_vec <- as.vector(W2[[i]][upper.tri(W2[[i]], diag = FALSE)])
    X2_vec[i, ] <- c(W2_vec)
  }
  
  ##Clean up the memory
  rm(X1_omega, X2_omega, W1, W2, W1_vec, W2_vec)
  gc()
  
  # 返回结果
  return(list(X1_vec = X1_vec, X2_vec = X2_vec))
}



