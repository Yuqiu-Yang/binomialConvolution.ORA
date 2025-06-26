rm(list=ls())
setwd("data")
library(stringr)
methods = c("gmm", "reg")
species = c('human', "ai")

for(method in methods)
{
  method_simulation_result = read.csv(paste0("./results/simulation/se/v_ratio_table_",method,".csv"))
  tp_scaling = mean(method_simulation_result$tp_mn_boot_1_v_ratio)
  tn_scaling = mean(method_simulation_result$tn_mn_boot_1_v_ratio)
  real_result = readRDS(paste0("./real_result_ci/",method,"_result.rds"))

  for(sp in species)
  {
    result = data.frame(matrix(0, nrow=11, ncol=7))
    colnames(result) = c("Passage",
                         paste0("tp_", c("se", "ll", "ul")),
                         paste0("tn_", c("se", "ll", "ul")))
    all_passages = c(paste0("Passage_", str_pad(1:10, 2, pad="0")),
                     "Combined")
    result$Passage = all_passages

    n_row=1
    for(passage in all_passages)
    {
      if(passage == "Combined")
      {
        tp_ci = quantile(real_result[[sp]]$same_est$mn_boot_1_est[,1], probs=c(0.025, 0.975), na.rm=TRUE)
        tn_ci = quantile(real_result[[sp]]$same_est$mn_boot_1_est[,2], probs=c(0.025, 0.975), na.rm=TRUE)

        tp_se = sqrt(var(real_result[[sp]]$same_est$mn_boot_1_est[,1], na.rm=TRUE)/tp_scaling)

        tn_se = sqrt(var(real_result[[sp]]$same_est$mn_boot_1_est[,2], na.rm=TRUE)/tn_scaling)
      }else{
        tp_ci = quantile(real_result[[sp]]$diff_est[[passage]]$mn_boot_1_est[,1], probs=c(0.025, 0.975), na.rm=TRUE)
        tn_ci = quantile(real_result[[sp]]$diff_est[[passage]]$mn_boot_1_est[,2], probs=c(0.025, 0.975), na.rm=TRUE)

        tp_se = sqrt(var(real_result[[sp]]$diff_est[[passage]]$mn_boot_1_est[,1], na.rm=TRUE)/tp_scaling)

        tn_se = sqrt(var(real_result[[sp]]$diff_est[[passage]]$mn_boot_1_est[,2], na.rm=TRUE)/tn_scaling)
      }
      result$tp_se[n_row] = tp_se
      result$tp_ul[n_row] = tp_ci[2]
      result$tp_ll[n_row] = tp_ci[1]
      result$tn_se[n_row] = tn_se
      result$tn_ul[n_row] = tn_ci[2]
      result$tn_ll[n_row] = tn_ci[1]
      n_row = n_row+1
    }
    write.csv(result, file = paste0("./results/real_data/",
                                    method,"_",sp,".csv"),
              row.names = F)
  }
}


#MLE
real_result = readRDS(paste0("./results/real_data/mle_result_w_ic.rds"))

for(sp in species)
{
  result = data.frame(matrix(0, nrow=11, ncol=7))
  colnames(result) = c("Passage",
                       paste0("tp_", c("se", "ll", "ul")),
                       paste0("tn_", c("se", "ll", "ul")))
  all_passages = c(paste0("Passage_", str_pad(1:10, 2, pad="0")),
                   "Combined")
  result$Passage = all_passages

  n_row=1
  for(passage in all_passages)
  {
    if(passage == "Combined")
    {
      tp_ci = c(real_result[[sp]]$same_est$est$pi.hat.ll[1],
                real_result[[sp]]$same_est$est$pi.hat.ul[1])
      tn_ci = c(real_result[[sp]]$same_est$est$pi.hat.ll[2],
                real_result[[sp]]$same_est$est$pi.hat.ul[2])
    }else{
      tp_ci = c(real_result[[sp]]$diff_est[[passage]]$est$pi.hat.ll[1],
                real_result[[sp]]$diff_est[[passage]]$est$pi.hat.ul[1])
      tn_ci = c(real_result[[sp]]$diff_est[[passage]]$est$pi.hat.ll[2],
                real_result[[sp]]$diff_est[[passage]]$est$pi.hat.ul[2])
    }
    tp_se = diff(tp_ci)/2/qnorm(0.975)
    tn_se = diff(tn_ci)/2/qnorm(0.975)
    result$tp_se[n_row] = tp_se
    result$tp_ul[n_row] = tp_ci[2]
    result$tp_ll[n_row] = tp_ci[1]
    result$tn_se[n_row] = tn_se
    result$tn_ul[n_row] = tn_ci[2]
    result$tn_ll[n_row] = tn_ci[1]
    n_row=n_row+1
  }
  write.csv(result, file = paste0("./results/real_data/",
                                  "mle_",sp,".csv"),
            row.names = F)
}
ic_result = data.frame(matrix(0, nrow=2, ncol=7))
colnames(ic_result) = c("Species",
                        "same_prob_lk", "diff_prob_lk",
                        "aic_same", "aic_diff",
                        "bic_same", "bic_diff")
ic_result$Species = species
ic_result[1, 2:7] = unlist(real_result$human$summary)
ic_result[2, 2:7] = unlist(real_result$ai$summary)
write.csv(ic_result, file='./results/real_data/ic_result.csv',
          row.names = F)








