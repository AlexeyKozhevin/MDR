library(doParallel)


mdr.efe <- function(dataset, r, K=5, type, prob = NA, parallel = F)
{
  fun <- switch(type,
         indep = mdr.indep(dataset,r,K,prob,parallel),
         fix = mdr.fix(dataset,r,K,prob,parallel)
    )  
  return (fun)
}
  
mdr.indep <- function(dataset,r,K,prob,parallel)
{  
  submodels <- t(combn(x = 1:n, m = r))  
  errors <- rep(0,nrow(submodels))
  N <- nrow(dataset)
  folds <- KFold(1:N,K)
  cases <- matrix(0, ncol = K, nrow = 3^r)
  controls <- matrix(0, ncol = K, nrow = 3^r)
  
  errors <- rep(0,nrow(submodels))
  
  if (parallel == F)
  {
    for (comb in 1:nrow(submodels))
    {
      data <- dataset[,c(1,submodels[comb,]+1)]
      for (k in 1:K)
      {
        res <- counter(data,folds,k,r)
        controls[,k] <- res[,1]
        cases[,k] <- res[,2]
      }
      
      Error <- 0
      
      for (k in 1:K)
      {
        Error <- Error + KFold.EFE.indep(controls,cases,k,prob)
      }
      
      errors[comb] <- Error/K
    }
  }
  else
  {
    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    registerDoParallel(cl)
    
    ls <- foreach (comb = 1:nrow(submodels),.packages="stats",.export=c("counter","KFold.EFE.indep","number_by_code")) %dopar%
    {
      data <- dataset[,c(1,submodels[comb,]+1)]
      for (k in 1:K)
      {
        res <- counter(data,folds,k,r)
        controls[,k] <- res[,1]
        cases[,k] <- res[,2]
      }
      
      Error <- 0
      
      for (k in 1:K)
      {
        Error <- Error + KFold.EFE.indep(controls,cases,k,prob)
      }
      Error/K
    }
    stopCluster(cl)
    for (i in 1:nrow(submodels))
    {
      errors[i] <- ls[[i]]
    }
  }
  fin <- submodels[which(errors == min(errors)),]
  
  if (is.numeric(dim(fin))) return (rep(NA,r))
  else return (fin)
}

mdr.fix <- function(dataset,r,K,prob,parallel)
{  
  submodels <- t(combn(x = 1:n, m = r))  
  errors <- rep(0,nrow(submodels))
  N <- nrow(dataset)
  folds_controls <- KFold(1:(N/2),K)
  folds_cases <- KFold(1:(N/2)+N/2,K)
  
  cases <- matrix(0, ncol = K, nrow = 3^r)
  controls <- matrix(0, ncol = K, nrow = 3^r)
  
  errors <- rep(0,nrow(submodels))
  
  if (parallel == F)
  {
    for (comb in 1:nrow(submodels))
    {
      data <- dataset[,c(1,submodels[comb,]+1)]
      for (k in 1:K)
      {
        res <- counter.fix(data,folds_controls,folds_cases,k,r)
        controls[,k] <- res[,1]
        cases[,k] <- res[,2]
      }
      
      Error <- 0
      
      for (k in 1:K)
      {
        Error <- Error + KFold.EFE.fix(controls,cases,k,prob)
      }
      
      errors[comb] <- Error/K
    }
  }
  else
  {
    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    registerDoParallel(cl)
    
    ls <- foreach (comb = 1:nrow(submodels),.packages="stats",.export=c("counter.fix",
                                                                        "KFold.EFE.fix","number_by_code")) %dopar%
    {
      data <- dataset[,c(1,submodels[comb,]+1)]
      for (k in 1:K)
      {
        res <- counter.fix(data,folds_controls,folds_cases,k,r)
        controls[,k] <- res[,1]
        cases[,k] <- res[,2]
      }
      
      Error <- 0
      
      for (k in 1:K)
      {
        Error <- Error + KFold.EFE.fix(controls,cases,k,prob)
      }
      
      Error/K
    }
    stopCluster(cl)
    for (i in 1:nrow(submodels))
    {
      errors[i] <- ls[[i]]
    }
    
  }
  
  fin <- submodels[which(errors == min(errors)),]
  
  if (is.numeric(dim(fin))) return (rep(NA,r))
  else return (fin)
}

KFold.EFE.indep <- function(controls,cases,k,prob)
{
  train_cases <- apply(cases[,-k],1,sum)
  train_controls <- apply(controls[,-k],1,sum)
  
  total_train_cases <- sum(train_cases)
  total_train_controls <- sum(train_controls)
  
  test_cases <- cases[,k]
  test_controls <- controls[,k]
  
  total_test_cases <- sum(test_cases)
  total_test_controls <- sum(test_controls)
  
  if (is.na(prob))
  {
    prob_y <- total_train_cases/(total_train_cases+total_train_controls)  
  } else
  {
    prob_y <- prob
  }
  
  #if psi is not for overline
  #psi1 <- 1/(total_test_cases/(total_test_cases+total_test_controls))
  
  #if psi is for overline
  psi1 <- 1/prob_y
  psi0 <- 1/(1-prob_y)
  
  f_PA <- 2*( train_cases/(train_cases+train_controls) > prob_y ) - 1
  f_PA[is.na(f_PA)] = 1
  
  Err <- (test_cases*psi1*(f_PA != 1) +
                    test_controls*psi0*(f_PA != -1))
  
  #print(Err)
  Err <- sum(Err,na.rm = T)
  
  return (Err/(total_test_cases+total_test_controls))
  
}

KFold.EFE.fix <- function(controls,cases,k,prob)
{
  train_cases <- apply(cases[,-k],1,sum)
  train_controls <- apply(controls[,-k],1,sum)
  
  total_train_cases <- sum(train_cases)
  total_train_controls <- sum(train_controls)
  
  test_cases <- cases[,k]
  test_controls <- controls[,k]
  
  total_test_cases <- sum(test_cases)
  total_test_controls <- sum(test_controls)
  
  f_PA <- 2*(train_cases/total_train_cases*prob/(train_controls/total_train_controls*(1-prob) + 
                                                         train_cases/total_train_cases*prob) > prob) - 1
  f_PA[is.na(f_PA)] = 1
  
  Err <- (test_cases*(f_PA != 1)/total_test_cases +
                    test_controls*(f_PA != -1)/total_test_controls)
  
  Err <- sum(Err,na.rm = T)
  
  return (Err)
  
}

counter <- function(dataset, folds, k, r)
{
  #counts and returns the number of cases and controls for k-th fold
  data <- dataset[folds[[k]],]
  ctrl_case <- matrix(0,nrow = 3^r, ncol = 2)
  
  indexes <- number_by_code(data[,-1])
  
  for (i in 1:nrow(data))
  {
    col <- (data[i,1]+1)/2+1
    ctrl_case[indexes[i],col] <- ctrl_case[indexes[i],col] + 1
  }
  return (ctrl_case)
}


counter.fix <- function(dataset, folds_controls, folds_cases, k, r)
{
  #counts and returns the number of cases and controls for k-th fold
  data_cases <- dataset[folds_cases[[k]],]
  data_controls <- dataset[folds_controls[[k]],]
  ctrl_case <- matrix(0,nrow = 3^r, ncol = 2)
  
  indexes_cases <- number_by_code(data_cases[,-1])
  indexes_controls <- number_by_code(data_controls[,-1])
  
  for (i in 1:nrow(data_cases))
  {
    ctrl_case[indexes_controls[i],1] <- ctrl_case[indexes_controls[i],1] + 1
    ctrl_case[indexes_cases[i],2] <- ctrl_case[indexes_cases[i],2] + 1
  }
  return (ctrl_case)
}

KFold <- function(array,K)
{
  N  <- length(array)
  len_fold <- (N%/%K)
  Folds <- list()
  if (K <= N)
  {
    for (k in 1:(K-1))
    {
      Folds[[k]] <- array[(k-1)*len_fold+(1:len_fold)]
    }
    Folds[[K]] <- array[((K-1)*len_fold+1):N]
    return (Folds)
  }
  else 
  {
    print("Not enough observations")
    return (NULL)
  }
}

number_by_code <- function(array)
{
  number <- 0
  power  <-  0
  for (i in 1:ncol(array))
  {
    number <- number + 3^power*array[,i]
    power  <- power+1
  }
  return (number+1)
}