
negbin = function(m,N,p)
{
  prob <- choose(m+N-1,m)*p^m*(1-p)^(N)
  return (prob)
}

sstr = function(C,w,p)
{
  N  <- 2
  totalprob <- 1
  while (totalprob > 0.95)
  {
    t_prev <- totalprob
    N1  <- N%/%2
    N_1 <- N - N1
    prob  <- 0
    m  <- N1
    while (m <= (w+1)/w*C-N/w-N_1)
    {
      prob = prob + negbin(m,N_1,p)
      m <- m+1
    }
    m <- N_1
    while (m <= (w+1)/w*C-N/w-N1)
    {
      prob = prob + negbin(m,N1,1-p)
      m <- m+1
    }  
    totalprob  <- prob
    #print(totalprob)
    N  <- N + 1
  }
  return (c(N-2,t_prev,N-1,totalprob))
}