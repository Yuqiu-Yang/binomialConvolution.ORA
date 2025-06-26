library(pracma)


estimate_mom <- function(passage_data,
                         passage_moments,
                         significance_level=0.05,
                         return_ci=TRUE)
{
  Y = passage_data$Y
  X = passage_data$X
  P = length(unique(passage_data$passage))

  mX = passage_moments$mX
  vX = passage_moments$vX
  cXY = passage_moments$cXY
  N = passage_moments$N
  mY = passage_moments$mY

  true_positive_prob = mY/N + cXY/vX*(1-mX/N)
  true_negative_prob = (1-mY/N) + cXY/vX*mX/N

  prob_est = c(true_positive_prob, true_negative_prob)
  prob_est = pmax(pmin(prob_est, 1),0)
  ul = ll = SE = NA
  if(return_ci)
  {
    Omega = matrix(0, nrow=4, ncol=4)
    counter = 1
    for(p in levels(passage_data$passage))
    {
      ind = which(passage_data$passage == p)
      data_omega = cbind(X[ind], Y[ind],
                         (X[ind] - passage_moments$mXp[counter])^2,
                         (X[ind]-passage_moments$mXp[counter])*(Y[ind]-passage_moments$mYp[counter]))
      Omega = Omega + cov(data_omega)/(passage_moments$np[counter]) * (passage_moments$wp[counter]^2)
      counter = counter + 1
    }

    A = rbind(
      c(-1/N*cXY/vX, 1/N, -cXY/vX^2*(1-mX/N),1/vX*(1-mX/N)),
      c(1/N*cXY/vX, -1/N, -cXY/vX^2*mX/N,1/vX*mX/N)
    )

    CV = A %*% Omega %*% t(A)
    SE = sqrt(diag(CV))

    ul = prob_est + qnorm(1-significance_level/2)*SE
    ll = prob_est - qnorm(1-significance_level/2)*SE
    ul = pmax(pmin(ul, 1),0)
    ll = pmax(pmin(ll, 1),0)
  }
  return(list("pi.hat"=prob_est,
              "pi.hat.ul"=ul,
              "pi.hat.ll"=ll,
              "SE"=SE))
}



estimate_gmm_helper <- function(par,
                                passage_data,
                                passage_moments,
                                Emp.Mom,
                                invOmega)
{
  Y = passage_data$Y
  X = passage_data$X

  n_nuisance_par = length(par)-2
  mu.X = unname(par[1:(n_nuisance_par/2)])
  var.X = unname(par[(n_nuisance_par/2+1):n_nuisance_par])
  pi.tp = par[n_nuisance_par+1]
  pi.tn = par[n_nuisance_par+2]

  mu.Y = mu.X * pi.tp + (passage_moments$N - mu.X) * (1 - pi.tn)
  var.Y = mu.X * pi.tp * (1 - pi.tp) + (passage_moments$N - mu.X) * pi.tn * (1 - pi.tn) +
    var.X * (pi.tp + pi.tn - 1)^2
  cov.XY = var.X * (pi.tp + pi.tn - 1)

  # n = nrow(passage_data)

  if (missing(Emp.Mom)) {
    # Then we compute for each passage the sample mean variance covariance
    Emp.Mom = c(passage_moments$mX,
                passage_moments$vX,
                passage_moments$mY,
                passage_moments$vY,
                passage_moments$cXY)
    # Emp.Mom = c(mean(X),var(X),mean(Y),var(Y),cov(X,Y))
    #Emp.Mom = c(mean(X),var(X),mean(Y),cov(X,Y))
  }

  wp = passage_moments$wp
  if (missing(invOmega)) {
    Omega = matrix(0, nrow=5, ncol=5)
    counter = 1
    for(p in levels(passage_data$passage))
    {
      ind = which(passage_data$passage == p)
      data_omega = cbind(X[ind], (X[ind] - passage_moments$mXp[counter])^2,
                         Y[ind], (Y[ind] - passage_moments$mYp[counter])^2,
                         (X[ind]-passage_moments$mXp[counter])*(Y[ind]-passage_moments$mYp[counter]))
      Omega = Omega + cov(data_omega)/(passage_moments$np[counter]) * (wp[counter]^2)
      counter = counter + 1
    }
    # data.omega = cbind(X,(X-mX)^2,Y,(Y-mY)^2,(X-mX)*(Y-mY))
    # Omega = cov(data.omega)
    invOmega = pinv(Omega)
  }

  Th.Mom = c(sum(mu.X*wp), sum(var.X*wp),
             sum(mu.Y*wp), sum(var.Y*wp),
             sum(cov.XY*wp))

  D = as.numeric(t(Emp.Mom - Th.Mom) %*% invOmega %*% (Emp.Mom - Th.Mom))

  return(D)

}


estimate_gmm <- function(passage_data,
                         passage_moments,
                         significance_level=0.05,
                         return_ci=TRUE)
{
  mom.est = estimate_mom(passage_data=passage_data,
                         passage_moments=passage_moments,
                         significance_level=significance_level,
                         return_ci=FALSE)
  P = length(unique(passage_data$passage))

  gmm_est = optim(c(passage_moments$mXp,
                    passage_moments$vXp,
                    mom.est$pi.hat),
                    fn = estimate_gmm_helper,
                    passage_data = passage_data,
                    passage_moments = passage_moments,
                    method = "L-BFGS-B",
                    lower = c(rep(0, 2*P),0,0),
                    upper = c(passage_moments$Np,
                              rep(Inf, P),1,1),
                    hessian = TRUE)

  prob_est = gmm_est$par[c(2*P+1,2*P+2)]
  prob_est = pmax(pmin(prob_est, 1),0)
  ul = ll = SE = NA
  if(return_ci)
  {
    SE = sqrt(diag(pinv(gmm_est$hessian)))[c(2*P+1,2*P+2)]

    ul = prob_est + qnorm(1-significance_level/2)*SE
    ll = prob_est - qnorm(1-significance_level/2)*SE
    ul = pmax(pmin(ul, 1),0)
    ll = pmax(pmin(ll, 1),0)
  }
  return(list("pi.hat"=prob_est,
              "pi.hat.ul"=ul,
              "pi.hat.ll"=ll,
              "SE"=SE))
}






