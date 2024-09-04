library(binomialConvolution)


same_prob_likelihood <- function(par, passage_data)
{
  probs = par_to_probs(par)
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
      lg = dbinom(x=as.numeric(df[r, "sample"]),
                  size=n_trials[keep],
                  prob=probs[keep],
                  log=TRUE)
    }else{
      lg = log_likelihood(success_probs=probs[keep],
                          n_trials = n_trials[keep],
                          samples = as.numeric(df[r, "sample"]))
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
                                              samples = as.numeric(df[r, "sample"]),
                                              n_trials = n_trials[keep],
                                              fixed_success_probs = fixed_success_probs,
                                              fixed_success_probs_index = fixed_success_probs_index)
    }else if(keep == fixed_success_probs_index){
      # if keep the first one and the first one is fixed
      lg = dbinom(x=as.numeric(df[r, "sample"]),
                  size=n_trials[keep],
                  prob=fixed_success_probs,
                  log=TRUE)
    }else{
      # if keep the first one but the second is fixed
      lg = dbinom(x=as.numeric(df[r, "sample"]),
                  size=n_trials[keep],
                  prob=par_to_probs(par),
                  log=TRUE)
    }
    sum_log = sum_log + lg
  }
  return(sum_log)
}
