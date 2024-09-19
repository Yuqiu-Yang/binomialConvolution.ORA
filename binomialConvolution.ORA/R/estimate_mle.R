library(binomialConvolution)


passage_likelihood <- function(par,
                               passage_data)
{
  probs = 1/(1+exp(-par))
  sum_log = 0
  for(r in 1 : nrow(passage_data))
  {
    n_trials = as.numeric(passage_data[r, c("n_positive", "n_negative")])
    if(n_trials[2] < 1)
    {
      keep = 1
    }else{
      keep = c(1,2)
    }
    if(length(keep) == 1)
    {
      lg = dbinom(x=as.numeric(passage_data[r, "Y"]),
                  size=n_trials[keep],
                  prob=probs[keep],
                  log=TRUE)
    }else{
      lg = log_likelihood(success_probs=probs[keep],
                          n_trials = n_trials[keep],
                          samples = as.numeric(passage_data[r, "Y"]))
    }
    sum_log = sum_log + lg
  }
  return(sum_log)
}


passage_profile_likelihood <- function(par,
                                       passage_data,
                                       fixed_success_probs,
                                       fixed_success_probs_index)
{
  sum_log = 0
  for(r in 1 : nrow(passage_data))
  {
    n_trials = as.numeric(passage_data[r, c("n_positive", "n_negative")])
    if(n_trials[2] < 1)
    {
      keep = 1
    }else{
      keep = c(1,2)
    }
    if(length(keep) > 1)
    {
      lg = log_likelihood_equality_constraint(par,
                                              samples = as.numeric(passage_data[r, "Y"]),
                                              n_trials = n_trials[keep],
                                              fixed_success_probs = fixed_success_probs,
                                              fixed_success_probs_index = fixed_success_probs_index)
    }else if(keep == fixed_success_probs_index){
      # if keep the first one and the first one is fixed
      lg = dbinom(x=as.numeric(passage_data[r, "Y"]),
                  size=n_trials[keep],
                  prob=fixed_success_probs,
                  log=TRUE)
    }else{
      # if keep the first one but the second is fixed
      lg = dbinom(x=as.numeric(passage_data[r, "Y"]),
                  size=n_trials[keep],
                  prob=1/(1+exp(-par)),
                  log=TRUE)
    }
    sum_log = sum_log + lg
  }
  return(sum_log)
}


profile_likelihood_ci <- function(fixed_success_probs,
                                  passage_data,
                                  fixed_success_probs_index,
                                  mle_llh,
                                  significance_level=0.05)
{
  test_llh = optimize(f=passage_profile_likelihood,
                     interval=c(-10,10),
                     passage_data=passage_data,
                     fixed_success_probs=fixed_success_probs,
                     fixed_success_probs_index=fixed_success_probs_index,
                     maximum=T)$objective
  lrt = -2 * (test_llh - mle_llh)
  return(lrt - qchisq(1-significance_level, 1))
}



estimate_mle <- function(passage_data,
                         passage_moments,
                         significance_level=0.05)
{
  mom.est = estimate_mom(passage_data=passage_data,
                         passage_moments=passage_moments,
                         significance_level=significance_level)
  probs=pmax(pmin(mom.est$pi.hat,0.99),0.001)
  probs[2]=1-probs[2]
  prob_est_optim = optim(par=log(probs/(1-probs)),
                        fn=passage_likelihood,
                        passage_data=passage_data,
                        control=list(fnscale=-1),
                        method = "BFGS")
  prob_est = 1/(1+exp(-prob_est_optim$par))
  llh = prob_est_optim$value
  # Profile confidence intervals based on the identical prob assumption
  p1_ul = tryCatch(
    expr = {
      uniroot(f=profile_likelihood_ci,
              lower=prob_est[1], upper=0.9999999,
              passage_data=passage_data,
              fixed_success_probs_index=1,
              mle_llh=llh,
              significance_level=significance_level)$root
    },
    error = function(e){
      return(1)
    }
  )
  p1_ll = tryCatch(
    expr = {
      uniroot(f=profile_likelihood_ci,
              lower=0.0000001, upper=prob_est[1],
              passage_data=passage_data,
              fixed_success_probs_index=1,
              mle_llh=llh,
              significance_level=significance_level)$root
    },
    error = function(e){
      return(0)
    }
  )

  p2_ll = tryCatch(
    expr = {
      1-uniroot(f=profile_likelihood_ci,
                lower=prob_est[2], upper=0.9999999,
                passage_data=passage_data,
                fixed_success_probs_index=2,
                mle_llh=llh,
                significance_level=significance_level)$root
    },
    error = function(e){
      return(0)
    }
  )


  p2_ul = tryCatch(
    expr = {
      1-uniroot(f=profile_likelihood_ci,
                lower=0.0000001, upper=prob_est[2],
                passage_data=passage_data,
                fixed_success_probs_index=2,
                mle_llh=llh,
                significance_level=significance_level)$root
    },
    error = function(e){
      return(1)
    }
  )


  return(list("pi.hat"=c(prob_est[1], 1-prob_est[2]),
              "pi.hat.ul"=c(p1_ul, p2_ul),
              "pi.hat.ll"=c(p1_ll, p2_ll)))
}


