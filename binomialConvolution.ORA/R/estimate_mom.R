library(dplyr)
library(pracma)

estimate_mom <- function(passage_data)
{
  Y = passage_data$Y
  X = passage_data$X
  P = length(unique(passage_data$passage))

  passage_data$passage = factor(passage_data$passage)

  # For each passage we compute the passage length and sample size
  Np = aggregate(passage_data$N, by=list(passage=passage_data$passage), FUN=function(x)x[1])$x
  np = aggregate(passage_data$passage, by=list(passag=passage_data$passage), FUN=function(x)length(x))$x
  # based on this we can compute the weights
  wp = Np*np/(sum(Np*np))

  # Then we compute for each passage the sample mean variance covariance
  mXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=mean)$x
  vXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=var)$x
  mYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=mean)$x
  # vYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=var)$x
  cov = passage_data %>% group_by(passage) %>% summarise(sig=cov(X,Y))
  cXYp = cov$sig

  mX = sum(mXp * wp)
  vX = sum(vXp * wp)
  cXY = sum(cXYp * wp)
  N = sum(Np * wp)
  mY = sum(mYp * wp)

  true_positive_prob = mY/N + cXY/vX*(1-mX/N)
  true_negative_prob = (1-mY/N) + cXY/vX*mX/N

  prob_est = c(true_positive_prob, true_negative_prob)

  Omega = matrix(0, nrow=4, ncol=4)
  counter = 1
  for(p in levels(passage_data$passage))
  {
    ind = which(passage_data$passage == p)
    data_omega = cbind(X[ind], Y[ind],
                       (X[ind] - mXp[counter])^2,
                       (X[ind]-mXp[counter])*(Y[ind]-mYp[counter]))
    Omega = Omega + cov(data_omega)/(np[counter]) * (wp[counter]^2)
    counter = counter + 1
  }

  A = rbind(
    c(-1/N*cXY/vX, 1/N, -cXY/vX^2*(1-mX/N),1/vX*(1-mX/N)),
    c(1/N*cXY/vX, -1/N, -cXY/vX^2*mX/N,1/vX*mX/N)
  )

  CV = A %*% Omega %*% t(A)
  SE = sqrt(diag(CV))
  return(list("pi.hat"=prob_est,
              "pi.hat.est"=SE))
}



estimate_gmm_helper <- function(par,passage_data,Emp.Mom,invOmega)
{

  # data must be structured as a list with array items X and Y
  # and scalar items N and n
  Y = passage_data$Y
  X = passage_data$X
  passage_data$passage = factor(passage_data$passage)

  Np = aggregate(passage_data$N, by=list(passage=passage_data$passage), FUN=function(x)x[1])$x
  np = aggregate(passage_data$passage, by=list(passag=passage_data$passage), FUN=function(x)length(x))$x
  # based on this we can compute the weights
  wp = Np*np/(sum(Np*np))

  # Then we compute for each passage the sample mean variance covariance
  mXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=mean)$x
  vXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=var)$x
  mYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=mean)$x
  vYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=var)$x
  cov = passage_data %>% group_by(passage) %>% summarise(sig=cov(X,Y))
  cXYp = cov$sig

  mX = sum(mXp * wp)
  vX = sum(vXp * wp)
  cXY = sum(cXYp * wp)
  N = sum(Np * wp)
  mY = sum(mYp * wp)
  vY = sum(vYp * wp)

  n_nuisance_par = length(par)-2
  mu.X = unname(par[1:(n_nuisance_par/2)])
  var.X = unname(par[(n_nuisance_par/2+1):n_nuisance_par])
  pi.tp = par[n_nuisance_par+1]
  pi.tn = par[n_nuisance_par+2]

  mu.Y = mu.X * pi.tp + (N - mu.X) * (1 - pi.tn)
  var.Y = mu.X * pi.tp * (1 - pi.tp) + (N - mu.X) * pi.tn * (1 - pi.tn) +
    var.X * (pi.tp + pi.tn - 1)^2
  cov.XY = var.X * (pi.tp + pi.tn - 1)

  # n = nrow(passage_data)

  if (missing(Emp.Mom)) {
    # Then we compute for each passage the sample mean variance covariance
    Emp.Mom = c(mX, vX, mY, vY, cXY)
    # Emp.Mom = c(mean(X),var(X),mean(Y),var(Y),cov(X,Y))
    #Emp.Mom = c(mean(X),var(X),mean(Y),cov(X,Y))
  }

  if (missing(invOmega)) {
    Omega = matrix(0, nrow=5, ncol=5)
    counter = 1
    for(p in levels(passage_data$passage))
    {
      ind = which(passage_data$passage == p)
      data_omega = cbind(X[ind], (X[ind] - mXp[counter])^2,
                         Y[ind], (Y[ind] - mYp[counter])^2,
                         (X[ind]-mXp[counter])*(Y[ind]-mYp[counter]))
      Omega = Omega + cov(data_omega)/(np[counter]) * (wp[counter]^2)
      counter = counter + 1
    }
    # data.omega = cbind(X,(X-mX)^2,Y,(Y-mY)^2,(X-mX)*(Y-mY))
    # Omega = cov(data.omega)
    invOmega = solve(Omega)
  }

  Th.Mom = c(sum(mu.X*wp), sum(var.X*wp),
             sum(mu.Y*wp), sum(var.Y*wp),
             sum(cov.XY*wp))

  D = as.numeric(t(Emp.Mom - Th.Mom) %*% invOmega %*% (Emp.Mom - Th.Mom))

  return(D)

}


estimate_gmm <- function(passage_data)
{
  mom.est = estimate_mom(passage_data)
  passage_data$passage = factor(passage_data$passage)
  P = length(unique(passage_data$passage))

  gmm_est = optim(c(aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=mean)$x,
                    aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=var)$x,
                    mom.est$pi.hat),
                    fn = estimate_gmm_helper,
                    passage_data = passage_data,
                    method = "L-BFGS-B",
                    lower = c(rep(0, 2*P),0,0),
                    upper = c(aggregate(passage_data$N, by=list(passage=passage_data$passage), FUN=function(x)x[1])$x,
                              rep(Inf, P),1,1),
                    hessian = TRUE)

  return(list("pi.hat"=gmm_est$par[c(2*P+1,2*P+2)],
              "pi.hat.est"=sqrt(diag(pinv(gmm_est$hessian)))[c(2*P+1,2*P+2)]))
}






