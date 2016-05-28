power <- function(truemodel,estmodel)
{
  array <- c(truemodel,estmodel)
  array <- unique(array)
  power <- c(0,0)
  power <- as.numeric(length(array)==length(truemodel))
  return (power)
}