rm(list=ls())

load("data/real_data/human_ORA.rda")
load("data/real_data/ai_ORA.rda")
human_ORA$passage = factor(human_ORA$passage)
ai_ORA$passage = factor(ai_ORA$passage)

## Remove two influential points
human_outlier = which(human_ORA$passage == 'Passage_08')[15]
ai_outlier = which(ai_ORA$passage == 'Passage_09')[31]

human_ORA = human_ORA[-human_outlier, ]
ai_ORA = ai_ORA[-ai_outlier, ]

codefolder = "code"
source(paste0(codefolder, "/utility.R"))
source(paste0(codefolder, "/estimate_mom.R"))
source(paste0(codefolder, "/estimate_regression.R"))
source(paste0(codefolder, "/estimate_mle.R"))

estimation_procedure = function(passage_data,
                                p_name=NA,
                                significance_level=0.05,
                                n_bootstrap=2000,
                                estimation_method="reg",
                                ...)
{
  if(estimation_method=="mom")
  {
    est_fun = estimate_mom
  }else if(estimation_method=="reg"){
    est_fun = estimate_linear_regression
  }else if(estimation_method=='gmm'){
    est_fun = estimate_gmm
  }else{
    est_fun = estimate_mle
  }
  result = list()
  if(!is.na(p_name))
  {
    passage_data = droplevels(passage_data[which(passage_data$passage == p_name), ])
  }
  if(estimation_method == "reg")
  {
    est = est_fun(passage_data=passage_data,
                  significance_level=significance_level,
                  return_ci=TRUE)
  }else{
    passage_moments = get_passage_moments(passage_data=passage_data)
    est = est_fun(passage_data=passage_data,
                  passage_moments=passage_moments,
                  significance_level=significance_level,
                  return_ci=TRUE,
                  ...)
  }

  mn_boot_1 = bootstrap_passages(passage_data=passage_data,
                                 true_positive_prob=NA,
                                 true_negative_prob=NA,
                                 n_bootstrap=n_bootstrap,
                                 sample_prob=1)

  mn_boot_1_est = matrix(0, nrow=n_bootstrap, ncol=2)

  ##### Regular
  for(n_boot in 1 : n_bootstrap)
  {
    if(estimation_method == "reg")
    {
      est_boot = est_fun(passage_data=mn_boot_1[[n_boot]],
                         significance_level=significance_level,
                         return_ci=FALSE)
    }else{
      passage_moments_boot = get_passage_moments(passage_data=mn_boot_1[[n_boot]])
      est_boot = est_fun(passage_data=mn_boot_1[[n_boot]],
                         passage_moments=passage_moments_boot,
                         significance_level=significance_level,
                         return_ci=FALSE, ...)
    }
    mn_boot_1_est[n_boot, ] = est_boot$pi.hat
  }

  result = list("est"=est,
                "mn_boot_1_est"=mn_boot_1_est)
  return(result)
}



estimation = function(passage_data,
                      estimation_method,
                      significance_level=0.05,
                      n_bootstrap=50)
{
  result = list("same_est"=list(),
                "diff_est"=list())

  result$same_est = estimation_procedure(passage_data=passage_data,
                                         p_name=NA,
                                         significance_level = significance_level,
                                         n_bootstrap = n_bootstrap,
                                         estimation_method=estimation_method)

  for(p_name in levels(passage_data$passage))
  {
    if(estimation_method=="mle")
    {
      result$diff_est[[p_name]] = estimation_procedure(passage_data=passage_data,
                                                       p_name=p_name,
                                                       significance_level = significance_level,
                                                       n_bootstrap=n_bootstrap,
                                                       estimation_method=estimation_method,
                                                       initial_prob=result$same_est$est$pi.hat)
    }else{
      result$diff_est[[p_name]] = estimation_procedure(passage_data=passage_data,
                                                       p_name=p_name,
                                                       significance_level = significance_level,
                                                       n_bootstrap=n_bootstrap,
                                                       estimation_method=estimation_method)
    }

  }
  return(result)
}

significance_level = 0.05
n_bootstrap = 2000
set.seed(42)
for(method in c("reg", "gmm", "mle"))
{
  result = list('human'=list(),
                "ai"=list())
  result$human = estimation(passage_data = human_ORA,
                            estimation_method = method,
                            significance_level = significance_level,
                            n_bootstrap=n_bootstrap)
  result$ai = estimation(passage_data = ai_ORA,
                         estimation_method = method,
                         significance_level = significance_level,
                         n_bootstrap=n_bootstrap)
  saveRDS(result, file=paste0("./real_result_ci/",method, "_result.rds"))
}

# For MLE, we further compute AIC and BIC

result = readRDS("./real_result_ci/mle_result.rds")

diff_prob_likelihood = function(par, passage_data)
{
  counter = 1
  sum_log = 0
  for(p_name in levels(passage_data$passage))
  {
    lg = passage_likelihood(par[(2*counter - 1) : (2*counter)],
                            droplevels(passage_data[which(passage_data$passage == p_name), ]))
    print(lg)
    counter = counter + 1
    sum_log = sum_log + lg
  }
  return(sum_log)
}

for(type in c("human", "ai"))
{
  if(type == 'human')
  {
    passage_data = human_ORA
  }else{
    passage_data = ai_ORA
  }
  same_est_prob = result[[type]]$same_est$est$pi.hat
  same_est_prob[2] = 1-same_est_prob[2]
  same_est_prob
  diff_est_prob = c()
  for(p_name in levels(passage_data$passage))
  {
    temp = result[[type]]$diff_est[[p_name]]$est$pi.hat
    temp[2] = 1-temp[2]
    diff_est_prob = c(diff_est_prob, temp)
  }
  diff_est_prob
  same_est_par = log(same_est_prob/(1-same_est_prob))
  # rep(same_est_par, 10)
  diff_est_par = log(diff_est_prob/(1-diff_est_prob))
  same_prob_lk = passage_likelihood(par=same_est_par,
                                    passage_data=passage_data)
  diff_prob_lk = diff_prob_likelihood(par=diff_est_par,
                                      passage_data=passage_data)
  print(diff_prob_lk > same_prob_lk)
  # Need to make sure that same < diff for both human and ai
  # if not change the initial values

  aic_same = -2 * same_prob_lk + 2 * 2
  aic_diff = -2 * diff_prob_lk + 2 * 20
  cat(aic_same, aic_diff)
  bic_same = -2 * same_prob_lk + 2 * log(nrow(passage_data))
  bic_diff = -2 * diff_prob_lk + 20 * log(nrow(passage_data))
  cat(bic_same, bic_diff)

  result[[type]]$summary = list("same_prob_lk"=same_prob_lk,
                                "diff_prob_lk"=diff_prob_lk,
                                "aic_same"=aic_same,
                                "aic_diff"=aic_diff,
                                "bic_same"=bic_same,
                                "bic_diff"=bic_diff)
}


saveRDS(result, file=paste0("./real_result_ci/mle_result_w_ic.rds"))

