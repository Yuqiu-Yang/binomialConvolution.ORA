
estimate_mom <- function(passage_data)
{
  Y = passage_data$Y
  X = passage_data$X
  N = passage_data$N[1]
  n = nrow(passage_data)

  true_positive_prob = mean(Y)/N + cov(X,Y)/var(X)*(1-mean(X)/N)
  true_negative_prob = (1-mean(Y)/N) + cov(X,Y)/var(X)*mean(X)/N
  prob_est = c(true_positive_prob, true_negative_prob)

  mX = mean(X)
  mY = mean(Y)
  vX = var(X)
  cXY = cov(X,Y)
  data_omega = cbind(X,Y,(X-mX)^2,(X-mX)*(Y-mY))
  Omega = cov(data_omega)

  A = rbind(
    c(-1/N*cXY/vX, 1/N, -cXY/vX^2*(1-mX/N),1/vX*(1-mX/N)),
    c(1/N*cXY/vX, -1/N, -cXY/vX^2*mX/N,1/vX*mX/N)
  )

  CV = A %*% Omega %*% t(A) / n

  SE = sqrt(diag(CV))

  result = c(prob_est, SE)
  return(result)
}



estimate_gmm_helper <- function(par,passage_data,Emp.Mom,invOmega)
{

  # data must be structured as a list with array items X and Y
  # and scalar items N and n

  mu.X = par[1]
  var.X = par[2]
  pi.tp = par[3]
  pi.tn = par[4]

  mu.Y = mu.X * pi.tp + (N - mu.X) * (1 - pi.tn)
  var.Y = mu.X * pi.tp * (1 - pi.tp) + (N - mu.X) * pi.tn * (1 - pi.tn) +
    var.X * (pi.tp + pi.tn - 1)^2
  cov.XY = var.X * (pi.tp + pi.tn - 1)

  Y = passage_data$Y
  X = passage_data$X
  N = passage_data$N[1]
  n = nrow(passage_data)

  if (missing(Emp.Mom)) {
    Emp.Mom = c(mean(X),var(X),mean(Y),var(Y),cov(X,Y))
    #Emp.Mom = c(mean(X),var(X),mean(Y),cov(X,Y))
  }

  if (missing(invOmega)) {
    mX = mean(X)
    mY = mean(Y)
    n = length(X)
    data.omega = cbind(X,(X-mX)^2,Y,(Y-mY)^2,(X-mX)*(Y-mY))
    Omega = cov(data.omega)
    invOmega = solve(Omega)
  }

  Th.Mom = c(mu.X, var.X, mu.Y, var.Y, cov.XY)

  D = n * as.numeric( t(Emp.Mom - Th.Mom) %*% invOmega %*% (Emp.Mom - Th.Mom) )

  return(D)

}


estimate_gmm <- function(passage_data)
{
  mom_est = estimate_mom(passage_data)
  gmm_est = optim(c(mean(X),var(X),mom.ests$pi.hat),
                    fn = estimate_gmm_helper,
                    passage_data = passage_data,
                    method = "L-BFGS-B",
                    lower = c(0,0,0,0),
                    upper = c(passage_data$N[1],Inf,1,1),
                    hessian = TRUE)
  result = c(gmm_est$par[c(3,4)], sqrt(diag(solve(gmm_est$hessian)))[c(3,4)])
  return(result)
}






