library(binomialConvolution)


same_prob_likelihood <- function(par,
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


same_prob_profile_likelihood <- function(par,
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
                                              samples = as.numeric(passage_data[r, "sample"]),
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
                  prob=par_to_probs(par),
                  log=TRUE)
    }
    sum_log = sum_log + lg
  }
  return(sum_log)
}



estimate_mle <- function(passage_data,
                         passage_moments,
                         significance_level=0.05)
{
  mom.est = estimate_mom(passage_data=passage_data,
                         passage_moments=passage_moments,
                         significance_level=significance_level)
  probs=mom.est$pi.hat
  same_prob_est = optim(par=log(probs/(1-probs)),
                        fn=same_prob_likelihood,
                        passage_data=passage_data,
                        control=list(fnscale=-1),
                        method = "BFGS")
  prob_est = 1/(1+exp(-same_prob_est$par))
  # Profile confidence intervals based on the identical prob assumption
  p_p1_grid = seq(0.9, 0.99999, length.out = 100)
  p_p1_nllh = numeric(100)
  for(k in 1 : 100)
  {
    p_p1_nllh[k] = -optimize(f=same_prob_profile_likelihood,
                             interval=c(-10,10),
                             passage_data=passage_data,
                             fixed_success_probs=p_p1_grid[k],
                             fixed_success_probs_index=1,
                             maximum=T)$objective
  }
  p_p2_grid = seq(0.1, 0.5, length.out = 100)
  p_p2_nllh = numeric(100)
  for(k in 1 : 100)
  {
    p_p2_nllh[k] = -optimize(f=same_prob_profile_likelihood,
                             interval=c(-10,10),
                             passage_data=passage_data,
                             fixed_success_probs=p_p2_grid[k],
                             fixed_success_probs_index=2,
                             maximum=T)$objective
  }

  2*(final_est$Y.Human$p_p1.nllh - min(final_est$Y.Human$p_p1.nllh))


}


