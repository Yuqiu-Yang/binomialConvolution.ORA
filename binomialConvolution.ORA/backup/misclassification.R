rm(list=ls())
library(binomialConvolution)
setwd("/Users/YuqiuYang/Desktop/repos/binomialConvolution/results/passages/")
load("./passage_list.RData")
same_prob_likelihood = function(par,
                                df)
{
  probs = par_to_probs(par)
  sum_log = 0
  for(r in 1 : nrow(df))
  {
    n_trials = as.numeric(df[r, c("n1", "n2")])
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
                                         df,
                                         fixed_success_probs,
                                         fixed_success_probs_index)
{
  sum_log = 0
  for(r in 1 : nrow(df))
  {
    n_trials = as.numeric(df[r, c("n1", "n2")])
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
diff_prob_likelihood = function(par, df)
{
  passages = sort(unique(df$passage))
  counter = 1
  sum_log = 0
  for(passage in passages)
  {
    lg = same_prob_likelihood(par[(2*counter - 1) : (2*counter)],
                              df[which(df$passage == passage), ])
    counter = counter + 1
    sum_log = sum_log + lg
  }
  return(sum_log)
}
final_est = list()
final_est[["Y.Human"]] = list()
final_est[["Y.AI"]] = list()
for(type in c("Y.Human", "Y.AI"))
{
  n_rows = sum(sapply(passage_list, FUN = nrow))
  df = data.frame(matrix(NA, nrow=n_rows, ncol=4))
  colnames(df) = c("passage", "sample", "n1", "n2")
  counter = 1
  for(passage in names(passage_list))
  {
    for(r in 1 : nrow(passage_list[[passage]]))
    {
      df[counter, "passage"] = passage
      df[counter, "sample"] = passage_list[[passage]][[type]][r]
      df[counter, "n1"] = passage_list[[passage]]$X[r]
      df[counter, "n2"] = passage_list[[passage]]$N[r] - df[counter, "n1"]
      counter = counter + 1
    }
  }
  same_prob_est = optim(par=c(1,-1),
                        fn=same_prob_likelihood,
                        df=df,
                        control=list(fnscale=-1),
                        method = "BFGS")
  diff_prob_est = optim(par=rep(same_prob_est$par,
                                10),
                        fn=diff_prob_likelihood,
                        df=df,
                        control=list(fnscale=-1),
                        method = "BFGS")
  same_probs = par_to_probs(same_prob_est$par)
  same_probs
  diff_probs = c()
  for(i in 1 : 10)
  {
    p = par_to_probs(diff_prob_est$par[(2*i-1) : (2*i)])
    diff_probs = c(diff_probs, p)
  }
  diff_probs
  aic_same = -2 * same_prob_est$value + 2 * 2
  aic_diff = -2 * same_prob_est$value + 2 * 20
  cat(aic_same, aic_diff)
  bic_same = -2 * diff_prob_est$value + 2 * log(nrow(df))
  bic_diff = -2 * diff_prob_est$value + 20 * log(nrow(df))
  cat(bic_same, bic_diff)
  # Profile confidence intervals based on the identical prob assumption
  p_p1_grid = seq(0.9, 0.99999, length.out = 100)
  p_p1_nllh = numeric(100)
  for(k in 1 : 100)
  {
    p_p1_nllh[k] = -optimize(f=same_prob_profile_likelihood,
                             interval=c(-10,10),
                             df=df,
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
                             df=df,
                             fixed_success_probs=p_p2_grid[k],
                             fixed_success_probs_index=2,
                             maximum=T)$objective
  }
  # Joint confidence region based on the identical prob assumption
  k.grid <- 100
  p1.grid <- seq(0.9,0.9999,length.out=k.grid)
  p2.grid <- seq(0.1,0.5,length.out=k.grid)
  nllh <- array(NA,dim=c(k.grid,k.grid))
  for (i in 1:k.grid) {
    for (j in 1:k.grid) {
      par = probs_to_par(c(p1.grid[i], p2.grid[j]))
      nllh[i,j] <- -same_prob_likelihood(par, df)
    }
  }
  # contour(p1.grid, p2.grid, 2*(nllh-min(nllh)),
  #         levels = qchisq(c(0.5,0.8,0.9,0.95,0.99),2),
  #         xlab = "p1", ylab = "p2", main = "Contour Plot",
  #         xlim=c(0.99,1), ylim=c(0.2,0.5))
  final_est[[type]] = list("same" = same_prob_est,
                           "diff" = diff_prob_est,
                           "same_probs" = same_probs,
                           "diff_probs" = diff_probs,
                           "same_ic" = c(aic_same, bic_same),
                           "diff_ic" = c(aic_diff, bic_diff),
                           "p1.grid" = p1.grid,
                           "p2.grid" = p2.grid,
                           "nllh" = nllh,
                           "p_p1.grid"=p_p1_grid,
                           "p_p1.nllh"=p_p1_nllh,
                           "p_p2.grid"=p_p2_grid,
                           "p_p2.nllh"=p_p2_nllh)
}
contour(final_est$Y.Human$p1.grid, final_est$Y.Human$p2.grid,
        2*(final_est$Y.Human$nllh-min(final_est$Y.Human$nllh)),
        levels = qchisq(c(0.5,0.8,0.9,0.95,0.99),2),
        xlab = "True Positive Probability",
        ylab = "False Positive Probability", main = "Contour Plot",
        xlim=c(0.97,1), ylim=c(0.25,0.45),
        col=2, lwd =2 )
contour(final_est$Y.AI$p1.grid, final_est$Y.AI$p2.grid,
        2*(final_est$Y.AI$nllh-min(final_est$Y.AI$nllh)),
        levels = qchisq(c(0.5,0.8,0.9,0.95,0.99),2),
        col=1, lwd =2 , add=T)
legend("top", col=c(1,2), legend = c("ai", "human"),
       lty=1, lwd = 2)
save(final_est,
     file="/Users/YuqiuYang/Desktop/repos/binomialConvolution/results/passages/error_est.RData")


plot(final_est$Y.Human$p_p1.grid,
     2*(final_est$Y.Human$p_p1.nllh - min(final_est$Y.Human$p_p1.nllh)),
     type='l', ylim=c(0,100),
     xlim = c(0.94, 1),
     xlab = "True Positive probability",
     ylab = "Chisq stat",
     lwd = 2, col = 2)
lines(final_est$Y.AI$p_p1.grid,
      2*(final_est$Y.AI$p_p1.nllh - min(final_est$Y.AI$p_p1.nllh)),
      lwd = 2, col = 1)
abline(h=qchisq(0.99, 1),
       lty = 2, lwd = 2, col='blue')
abline(h=qchisq(0.95, 1),
       lty = 2, lwd = 2, col="green")
legend("topleft", col=c(1,2), legend = c("ai", "human"),
       lty=1, lwd = 2)



plot(final_est$Y.Human$p_p2.grid,
     2*(final_est$Y.Human$p_p2.nllh - min(final_est$Y.Human$p_p2.nllh)),
     type='l', ylim=c(0,100),
     xlim = c(0.1, 0.5),
     xlab = "False Positive probability",
     ylab = "Chisq stat",
     lwd = 2, col = 2)
lines(final_est$Y.AI$p_p2.grid,
      2*(final_est$Y.AI$p_p2.nllh - min(final_est$Y.AI$p_p2.nllh)),
      lwd = 2, col = 1)
abline(h=qchisq(0.99, 1),
       lty = 2, lwd = 2, col='blue')
abline(h=qchisq(0.95, 1),
       lty = 2, lwd = 2, col="green")
legend("top", col=c(1,2), legend = c("ai", "human"),
       lty=1, lwd = 2)
