library(dplyr)



get_passage_moments <- function(passage_data,
                                moments_to_compute=c("mXp", "vXp", "mYp", "vYp", "cXYp", "Np", "np"),
                                flat=FALSE)
{
  result = vector(mode='list', length = length(moments_to_compute))
  names(result) = moments_to_compute

  if("mXp" %in% moments_to_compute)
  {
    result$mXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=mean)$x
  }
  if("vXp" %in% moments_to_compute)
  {
    result$vXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=var)$x
  }
  if("mYp" %in% moments_to_compute)
  {
    result$mYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=mean)$x
  }
  if("vYp" %in% moments_to_compute)
  {
    result$vYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=var)$x
  }
  if("Np" %in% moments_to_compute)
  {
    result$Np = aggregate(passage_data$N, by=list(passage=passage_data$passage), FUN=function(x)x[1])$x
  }
  if("np" %in% moments_to_compute)
  {
    result$np = aggregate(passage_data$passage, by=list(passag=passage_data$passage), FUN=function(x)length(x))$x
  }
  if('cXYp' %in% moments_to_compute)
  {
    cov = passage_data %>% group_by(passage) %>% summarise(sig=cov(X,Y))
    result$cXYp = cov$sig
  }
  if(flat)
  {
    return(unlist(result))
  }else{
    return(result)
  }

}


estimate_mom <- function(passage_data)
{
  Y = passage_data$Y
  X = passage_data$X
  P = length(unique(passage_data$passage))
  H = diag(P)/P - matrix(1, nrow=P, ncol=P)/(P^2)

  passage_sample_moments = get_passage_moments(passage_data=passage_data,
                                               moments_to_compute=c("mXp", "vXp", 'Np', "np"))

  mX = mean(passage_sample_moments$mXp)
  N = mean(passage_sample_moments$Np)

  vX = mean(passage_sample_moments$vXp)
  sigma1 = as.numeric(vX + t(passage_sample_moments$mXp) %*% H %*% (passage_sample_moments$mXp))
  sigma2 = as.numeric(t(passage_sample_moments$mXp) %*% H %*% (passage_sample_moments$Np-passage_sample_moments$mXp) - vX)

  cXY = cov(X,Y)
  mY = mean(Y)

  true_negative_prob = 1 - (mY-cXY*mX/sigma1)/(N-mX*(1+sigma2/sigma1))
  true_positive_prob = cXY/sigma1 - (sigma2/sigma1)*(1-true_negative_prob)

  # true_positive_prob = mY/N + cXY/vX*(1-mX/N)
  # true_negative_prob = (1-mY/N) + cXY/vX*mX/N
  prob_est = c(true_positive_prob, true_negative_prob)

  if(P > 1)
  {
    bootstrap_sample = bootstrap_passages(passage_data=passage_data,
                                          true_positive_prob=NA,
                                          true_negative_prob=NA,
                                          sample_prob=1)
    mXs = bind_rows(lapply(bootstrap_sample, FUN=get_passage_moments, moments_to_compute=c("mXp"), flat=T))
    temp = lapply(bootstrap_sample, FUN=get_passage_moments, moments_to_compute=c("vXp"), flat=T)
    vXs = sapply(temp, FUN=mean)
    mYs = sapply(bootstrap_sample, FUN=function(df) {return(mean(df$Y))})
    cXYs = sapply(bootstrap_sample, FUN=function(df) {return(cov(df$X, df$Y))})
    data_omega = cbind(mXs, vXs, mYs, cXYs)
    Omega = cov(data_omega)

    partial_tn_mX=-(-cXY*(1/(P*sigma1)-2*mX*(t(passage_sample_moments$mXp)%*%H)/(sigma1^2))/(N-mX*(1-sigma2/sigma1)) +
      (mY-cXY*mX/sigma1)/(-(N-mX*(1-sigma2/sigma1))^2)*(-(1/P*(1+sigma2/sigma1)+mX*((t(passage_sample_moments$Np)%*%H-2*t(passage_sample_moments$mXp)%*%H)/sigma1-
                                                                                      2*sigma2/(sigma1)^2*t(passage_sample_moments$mXp)%*%H))))
    partial_tp_mX=-2*cXY*(t(passage_sample_moments$mXp)%*%H)/(sigma1^2) - (((t(passage_sample_moments$Np)%*%H-2*t(passage_sample_moments$mXp)%*%H)/sigma1-
                                                                               2*sigma2/(sigma1)^2*t(passage_sample_moments$mXp)%*%H)*((mY-cXY*mX/sigma1)/(N-mX*(1+sigma2/sigma1)))+
                                                                               (sigma2/sigma1)*(-partial_tn_mX))

    partial_tn_vX=-(cXY*mX/(sigma1^2)/(N-mX*(1+sigma2/sigma1))-((mY-cXY*mX/sigma1)/(N-mX*(1+sigma2/sigma1))^2)*(mX*(1/sigma1+sigma2/(sigma1^2))))
    partial_tp_vX=-cXY/(sigma1^2)-((-1/sigma1-sigma2/(sigma1)^2)*(mY-cXY*mX/sigma1)/(N-mX*(1+sigma2/sigma1)) + (sigma2/sigma1)*(-partial_tn_vX))

    partial_tn_mY = -1/(N-mX*(1+sigma2/sigma1))
    partial_tp_mY = sigma2/sigma1*partial_tn_mY

    partial_tn_cXY = mX/sigma1/(N-mX*(1+sigma2/sigma1))
    partial_tp_cXY = 1/sigma1+sigma2/sigma1*partial_tn_cXY

    A = rbind(c(partial_tp_mX, partial_tp_vX, partial_tp_mY, partial_tp_cXY),
              c(partial_tn_mX, partial_tn_vX, partial_tn_mY, partial_tn_cXY))
  }else{
    data_omega = cbind(X,Y,(X-mX)^2,(X-mX)*(Y-mY))
    Omega = cov(data_omega)/np
    A = rbind(
      c(-1/N*cXY/vX, 1/N, -cXY/vX^2*(1-mX/N),1/vX*(1-mX/N)),
      c(1/N*cXY/vX, -1/N, -cXY/vX^2*mX/N,1/vX*mX/N)
    )
  }
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
  N = unname(get_passage_moments(passage_data, "Np", flat=T))


  n_nuisance_par = length(par)-2
  mu.X = unname(par[1:(n_nuisance_par/2)])
  var.X = unname(par[(n_nuisance_par/2+1):n_nuisance_par])
  pi.tp = par[n_nuisance_par+1]
  pi.tn = par[n_nuisance_par+2]

  mu.Y = mu.X * pi.tp + (N - mu.X) * (1 - pi.tn)
  var.Y = mu.X * pi.tp * (1 - pi.tp) + (N - mu.X) * pi.tn * (1 - pi.tn) +
    var.X * (pi.tp + pi.tn - 1)^2
  cov.XY = var.X * (pi.tp + pi.tn - 1)



  n = nrow(passage_data)

  if (missing(Emp.Mom)) {
    Emp.Mom = unname(get_passage_moments(passage_data,
                                  moments_to_compute=c("mXp", "vXp", "mYp", 'vYp', 'cXYp'), flat=T))
    # Emp.Mom = c(mean(X),var(X),mean(Y),var(Y),cov(X,Y))
    #Emp.Mom = c(mean(X),var(X),mean(Y),cov(X,Y))
  }

  if (missing(invOmega)) {

    if(P > 1)
    {
      bootstrap_sample = bootstrap_passages(passage_data=passage_data,
                                            true_positive_prob=NA,
                                            true_negative_prob=NA,
                                            sample_prob=1)
      data.omega = bind_rows(lapply(bootstrap_sample,
                                    FUN=get_passage_moments,
                                    moments_to_compute=c("mXp", "vXp", "mYp", 'vYp', 'cXYp'),
                                    flat=T))
      Omega = cov(data.omega)
    }else{
      mX = mean(X)
      mY = mean(Y)
      n = length(X)
      data.omega = cbind(X,(X-mX)^2,Y,(Y-mY)^2,(X-mX)*(Y-mY))
      Omega = cov(data.omega)
    }
    invOmega = solve(Omega)
  }

  Th.Mom = c(mu.X, var.X, mu.Y, var.Y, cov.XY)

  D = n * as.numeric(t(Emp.Mom - Th.Mom) %*% invOmega %*% (Emp.Mom - Th.Mom))

  return(D)

}


estimate_gmm <- function(passage_data)
{
  mom.est = estimate_mom(passage_data)

  P = length(unique(passage_data$passage))

  gmm_est = optim(c(get_passage_moments(passage_data, "mXp",T),
                    get_passage_moments(passage_data, "vXp",T),
                    mom.est$pi.hat),
                    fn = estimate_gmm_helper,
                    passage_data = passage_data,
                    method = "L-BFGS-B",
                    lower = c(rep(0, 2*P),0,0),
                    upper = c(get_passage_moments(passage_data, "Np",T),rep(Inf, P),1,1),
                    hessian = TRUE)

  return(list("pi.hat"=gmm_est$par[c(2*P+1,2*P+2)],
              "pi.hat.est"=sqrt(diag(solve(gmm_est$hessian)))[c(2*P+1,2*P+2)]))
}






