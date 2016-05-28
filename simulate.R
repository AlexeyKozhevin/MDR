library("MASS")

diseaseXOR1 <- function(coef.number, features)
{
  p <- 0.05
  y <- 2*(rbinom(1,1,(sum(features[coef.number])%%2)*p))-1
  
  return (y)
}


diseaseXOR2 <- function(coef.number, features)
{
  p <- 0.1
  y <- 2*(rbinom(1,1,(sum(features[coef.number])%%2)*p))-1
  
  return (y)
}


diseaseXOR3 <- function(coef.number, features)
{
  p <- 0.2
  y <- 2*(rbinom(1,1,(sum(features[coef.number])%%2)*p))-1
  
  return (y)
}


list_of_diseases <- list(diseaseXOR1,diseaseXOR2,diseaseXOR3)
list_of_MAFs <- list(rep(0.5,3),rep(0.5,3),rep(0.5,3))

simulated_SNPs <- function(dim,numberOfObs,cost,coef.number,model,ro = 0,size=rep(1,dim),type)
{
  fun <- switch(type,
    indep = simulated_SNPs.indep(dim,numberOfObs,coef.number,model,ro = 0,size=rep(1,dim)),
    fix = simulated_SNPs.fix(dim,numberOfObs,coef.number,model,ro = 0,size=rep(1,dim)),
    comb = simulated_SNPs.comb(dim,numberOfObs,cost,coef.number,model,ro = 0,size=rep(1,dim))
    )
  return (fun)
}

simulated_SNPs.indep <- function(dim,numberOfObs,coef.number,model,ro = 0,size=rep(1,dim))
{
  
  disease = function(coef.number, data)
  {
    if (is.vector(data)) data <- matrix(data,nrow = 1)
    resp <- rep(0,nrow(data))
    for (i in 1:nrow(data))
    {
      list_of_diseases[[model]](coef.number, data[i,]) -> resp[i]
    }
    return (resp)
  }
  
  MAF <- runif(n,min = 0.05,max = 0.5)
  MAF[coef.number] <- list_of_MAFs[[model]]
  
  leveldown <- rep(0,dim)
  levelup <- rep(0,dim)
  
  for (pos in 1:dim)
  {
    leveldown[pos] <- qnorm((1-MAF[pos])^2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    levelup[pos] <- qnorm(MAF[pos]^2, mean = 0, sd = 1, lower.tail = F, log.p = FALSE)
  }
  
  Sigma <- matrix (NA,size[1],size[1])
  
  for (i in 1:size[1])
  {
    for (j in 1:size[1])
    {
      ifelse(i == j,Sigma[i,j] <- 1, Sigma[i,j] <- ro)  
    }
  }
  dim <- size[1]
  
  for (i in 2:length(size))
  {
    
    nullmatrix <- matrix(0,dim,size[i])
    Sigma <- cbind(Sigma,nullmatrix)
    tmp <- matrix(NA, size[i],size[i])
    for (k in 1:size[i])
    {
      for (j in 1:size[i])
      {
        ifelse(k == j,tmp[k,j] <- 1, tmp[k,j] <- ro)  
      }
    }
    tmp <- cbind(t(nullmatrix),tmp)
    Sigma <- rbind(Sigma, tmp)
    dim <- dim + size[i]
  }
  
  mu <-rep(0,dim)
  x <- mvrnorm(n=numberOfObs, mu, Sigma)
  
  #sequence of SNP
  
  for( i in seq_len(nrow(x)))
  {
    for( j in seq_len(ncol(x)))
    {
      ifelse(x[i,j] > levelup[j], x[i,j] <- 2, ifelse(x[i,j] < leveldown[j], x[i,j] <- 0, x[i,j] <- 1)) 
    }
  }
  
  y <- disease(coef.number,x)
  
  return (cbind(y,x))
}

simulated_SNPs.comb <- function(dim,numberOfObs,cost,coef.number,model,ro = 0,size=rep(1,dim))
{
  
  disease = function(coef.number, data)
  {
    if (is.vector(data)) data <- matrix(data,nrow = 1)
    resp <- rep(0,nrow(data))
    for (i in 1:nrow(data))
    {
      list_of_diseases[[model]](coef.number, data[i,]) -> resp[i]
    }
    return (resp)
  }
  
  
  MAF <- runif(n,min = 0.05,max = 0.5)
  MAF[coef.number] <- list_of_MAFs[[model]]
  
  
  leveldown <- rep(0,dim)
  levelup <- rep(0,dim)
  
  for (pos in 1:dim)
  {
    leveldown[pos] <- qnorm((1-MAF[pos])^2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    levelup[pos] <- qnorm(MAF[pos]^2, mean = 0, sd = 1, lower.tail = F, log.p = FALSE)
  }
  
  Sigma <- matrix (NA,size[1],size[1])
  
  for (i in 1:size[1])
  {
    for (j in 1:size[1])
    {
      ifelse(i == j,Sigma[i,j] <- 1, Sigma[i,j] <- ro)  
    }
  }
  dim <- size[1]
  
  for (i in 2:length(size))
  {
    
    nullmatrix <- matrix(0,dim,size[i])
    Sigma <- cbind(Sigma,nullmatrix)
    tmp <- matrix(NA, size[i],size[i])
    for (k in 1:size[i])
    {
      for (j in 1:size[i])
      {
        ifelse(k == j,tmp[k,j] <- 1, tmp[k,j] <- ro)  
      }
    }
    tmp <- cbind(t(nullmatrix),tmp)
    Sigma <- rbind(Sigma, tmp)
    dim <- dim + size[i]
  }
  
  #print(Sigma)
  
  mu <-rep(0,dim)
  
  number_of_1 <- c(0,0)
  number_of_2 <- c(0,0)
  
  cntrl_cases_1 <- list(matrix(0,ncol = dim+1,nrow = numberOfObs/2),matrix(0,ncol = dim+1,nrow = numberOfObs/2))
  cntrl_cases_2 <- list(matrix(0,ncol = dim+1,nrow = cost/2),matrix(0,ncol = dim+1,nrow = cost/2))
  independent_sample <- matrix(0,ncol = dim+1,nrow = cost)
  
  prob <- 0
  
  
  condition <-  T
  
  N <- 0

  while (condition)
  {
    N <- N+1
    x <- mvrnorm(1, mu, Sigma)
    x <- 1 + 1*(x > levelup) - 1*(x < leveldown)
    y <- disease(coef.number,x)
    
    if (N <= cost) independent_sample[N,] <- c(y,x)
    
    ind <- (y+2)/2+1
    
    if (number_of_1[ind] < numberOfObs/2)
    {
        number_of_1[ind] <- number_of_1[ind]+1
        cntrl_cases_1[[ind]][number_of_1[ind],] <- c(y,x)
        if (all(number_of_1 == numberOfObs/2))
        {
          tildeN1 <- N
        }
    }
    
    if (number_of_2[ind] < cost/2)
    {
      number_of_2[ind] <- number_of_2[ind]+1
      cntrl_cases_2[[ind]][number_of_2[ind],] <- c(y,x)
      if (all(number_of_2 == cost/2))
      {
        tildeN2 <- N
      }
    }
    
    condition <- (N < cost)||(any(number_of_1 < numberOfObs/2))||(any(number_of_2 < cost/2))
    prob <- prob + (y == 1)
  }
  
  prob <- prob/N
  
  return (list(independent_sample,rbind(cntrl_cases_1[[1]],cntrl_cases_1[[2]]),
               rbind(cntrl_cases_2[[1]],cntrl_cases_2[[2]]),prob,tildeN1,tildeN2))
}