nPr=function(n,r)
{
  return(factorial(n)/factorial(n-r))
}

nCr=function(n,r)
{
  return(choose(n,r))
}

SimPerm=function(n_vec)
{
  return(factorial(sum(n_vec))/prod(factorial(n_vec)))
}

