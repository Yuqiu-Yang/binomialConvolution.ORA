rm(list=ls())
setwd("data")

summarize_est <- function(est_result, true_val, n_bootstrap=0)
{
  bias = mean(est_result$est) - true_val
  variance = var(est_result$est)
  rmse = sqrt(mean((est_result$est - true_val)^2))
  coverage = mean((est_result$ll <= true_val) & (est_result$ul >= true_val))
  if(n_bootstrap >= 1)
  {
    semi_par_coverage = mean((est_result$semi_par_ll <= true_val) & (est_result$semi_par_ul >= true_val))
    mn_boot_1_coverage = mean((est_result$mn_boot_1_ll <= true_val) & (est_result$mn_boot_1_ul >= true_val))
    mn_boot_2sqrt_coverage = mean((est_result$mn_boot_2sqrt_ll <= true_val) & (est_result$mn_boot_2sqrt_ul >= true_val))
    mn_boot_23_coverage = mean((est_result$mn_boot_23_ll <= true_val) & (est_result$mn_boot_23_ul >= true_val))
    semi_par_v = mean(est_result$semi_par_v)
    mn_boot_1_v = mean(est_result$mn_boot_1_v)
    mn_boot_2sqrt_v = mean(est_result$mn_boot_2sqrt_v)
    mn_boot_23_v = mean(est_result$mn_boot_23_v)
  }else{
    semi_par_coverage = mn_boot_1_coverage = mn_boot_2sqrt_coverage = mn_boot_23_coverage =
    semi_par_v =
    mn_boot_1_v =
    mn_boot_2sqrt_v =
    mn_boot_23_v = NA
  }
  delta_v = mean(est_result$delta_v, na.rm=TRUE, trim=0.01)
  return(c(bias, variance,
           rmse, coverage,
           semi_par_coverage, mn_boot_1_coverage,
           mn_boot_2sqrt_coverage, mn_boot_23_coverage,
           delta_v, semi_par_v, mn_boot_1_v,
           mn_boot_2sqrt_v, mn_boot_23_v))
}


simulation_setting = read.csv("./simulation/se/simulation_setting.csv")

for(method in c("gmm", "reg", "mle"))
{
  tp_est_summary = tn_est_summary = data.frame(matrix(0, nrow=nrow(simulation_setting),
                                                      ncol=13))
  colnames(tp_est_summary) = colnames(tn_est_summary) = c("bias", "variance", "rmse",
                                                          "coverage",
                                                          "semi_par_coverage",
                                                          "mn_boot_1_coverage",
                                                          "mn_boot_2sqrt_coverage",
                                                          "mn_boot_23_coverage",
                                                          "delta_v", "semi_par_v",
                                                          "mn_boot_1_v",
                                                          "mn_boot_2sqrt_v", "mn_boot_23_v")
  est_result = list()
  for(i_setting in 1 : nrow(simulation_setting))
  {
    est_result[[paste0("setting_", i_setting)]] = list()

    completed_simulation = list.files(paste0("./simulation/setting_", i_setting, "/passage_est"), pattern=method)
    completed_simulation = gsub("passage_", "", completed_simulation)
    completed_simulation = gsub(paste0("_", method), "", completed_simulation)
    completed_simulation = gsub("_result.rds", "", completed_simulation)
    if(length(completed_simulation) < 1)
    {
      next
    }
    completed_simulation = as.numeric(completed_simulation)
    n_completed_simulation = max(completed_simulation)

    tp_est_result = tn_est_result = data.frame(matrix(0,
                                                      nrow = n_completed_simulation,
                                                      ncol=16))
    colnames(tp_est_result) = colnames(tn_est_result) = c("est",
                                                          "ll", "ul",
                                                          "semi_par_ll", "semi_par_ul",
                                                          "mn_boot_1_ll", "mn_boot_1_ul",
                                                          "mn_boot_2sqrt_ll", "mn_boot_2sqrt_ul",
                                                          "mn_boot_23_ll", "mn_boot_23_ul",
                                                          "delta_v", "semi_par_v",
                                                          "mn_boot_1_v",
                                                          "mn_boot_2sqrt_v", "mn_boot_23_v")

    for(i_simulation in 1 : n_completed_simulation)
    {
      # print(i_simulation)
      est = readRDS(paste0("./simulation/se/setting_",
                           i_setting, "/passage_est/passage_",
                           i_simulation,"_",method,"_result.rds"))
      for(prob_type in c("tp", "tn"))
      {
        if(prob_type == "tp")
        {
          ind = 1
        }else{
          ind = 2
        }

        prob_est = est[[paste0(method,"_est")]]$pi.hat[ind]
        prob_est_ll = ifelse(is.na(est[[paste0(method,"_est")]]$pi.hat.ll[ind]),
                             est[[paste0(method,"_est")]]$pi.hat[ind],
                             est[[paste0(method,"_est")]]$pi.hat.ll[ind])
        prob_est_ul = ifelse(is.na(est[[paste0(method,"_est")]]$pi.hat.ul[ind]),
                             est[[paste0(method,"_est")]]$pi.hat[ind],
                             est[[paste0(method,"_est")]]$pi.hat.ul[ind])
        delta_v = (est[[paste0(method,"_est")]]$SE[ind])^2

        if(simulation_setting$n_bootstrap[i_setting] >= 1)
        {
          semi_par_v = var(est$semi_par_est[,ind], na.rm = T)
          semi_par_se = sqrt(semi_par_v)
          semi_par_ll = prob_est - qnorm(0.975) * semi_par_se
          semi_par_ul = prob_est + qnorm(0.975) * semi_par_se

          mn_boot_1_v = var(est$mn_boot_1_est[,ind], na.rm = T)
          mn_boot_1_se = sqrt(mn_boot_1_v)
          mn_boot_1_ll = prob_est - qnorm(0.975) * mn_boot_1_se
          mn_boot_1_ul = prob_est + qnorm(0.975) * mn_boot_1_se

          mn_boot_2sqrt_v = var(est$mn_boot_2sqrt_est[,ind], na.rm = T)*2/sqrt(simulation_setting$n_students[i_setting])
          mn_boot_2sqrt_se = sqrt(mn_boot_2sqrt_v)
          mn_boot_2sqrt_ll = prob_est - qnorm(0.975) * mn_boot_2sqrt_se
          mn_boot_2sqrt_ul = prob_est + qnorm(0.975) * mn_boot_2sqrt_se

          mn_boot_23_v = var(est$mn_boot_23_est[,ind], na.rm = T)*2/3
          mn_boot_23_se = sqrt(mn_boot_23_v)
          mn_boot_23_ll = prob_est - qnorm(0.975) * mn_boot_23_se
          mn_boot_23_ul = prob_est + qnorm(0.975) * mn_boot_23_se
        }else{
          semi_par_v = semi_par_ll = semi_par_ul = NA

          mn_boot_1_v = mn_boot_1_ll = mn_boot_1_ul = NA

          mn_boot_2sqrt_v = mn_boot_2sqrt_ll = mn_boot_2sqrt_ul = NA

          mn_boot_23_v = mn_boot_23_ll = mn_boot_23_ul = NA
        }


        est_summary = c(prob_est, prob_est_ll, prob_est_ul,
                        semi_par_ll, semi_par_ul,
                        mn_boot_1_ll, mn_boot_1_ul,
                        mn_boot_2sqrt_ll, mn_boot_2sqrt_ul,
                        mn_boot_23_ll, mn_boot_23_ul,
                        delta_v, semi_par_v,
                        mn_boot_1_v,
                        mn_boot_2sqrt_v, mn_boot_23_v)

        if(prob_type == "tp")
        {
          tp_est_result[i_simulation, ]=est_summary
        }else{
          tn_est_result[i_simulation, ]=est_summary
        }
      }
    }

    for(i in 1 : 11)
    {
      tp_est_result[,i] = pmin(1, pmax(0, tp_est_result[,i]))
      tn_est_result[,i] = pmin(1, pmax(0, tn_est_result[,i]))
    }
    est_result[[paste0("setting_", i_setting)]][['tp_est_result']] = tp_est_result
    est_result[[paste0("setting_", i_setting)]][['tn_est_result']] = tn_est_result

    tp_true = simulation_setting$true_positive_prob[i_setting]
    tn_true = simulation_setting$true_negative_prob[i_setting]
    tp_est_summary[i_setting, ] = summarize_est(tp_est_result, tp_true, simulation_setting$n_bootstrap[i_setting])
    tn_est_summary[i_setting, ] = summarize_est(tn_est_result, tn_true, simulation_setting$n_bootstrap[i_setting])
  }
  write.csv(tp_est_summary, file=paste0("./results/simulation/se/",method,"_tp_est_summary.csv"),
            row.names = F)
  write.csv(tn_est_summary, file=paste0("./results/simulation/se/",method,"_tn_est_summary.csv"),
            row.names = F)
}

###################################


rm(list=ls())
setwd("data/results/simulation/se/")

for(method in c("gmm", "reg", "mle"))
{
  for(prob_type in c("tp", "tn"))
  {
    file_name = paste0("./", method,"_",prob_type,"_est_summary.csv")
    result = read.csv(file_name)
    v_mode = colnames(result)[grep("_v", colnames(result))]
    for(v_m in v_mode)
    {
      result[[paste0(v_m, "_ratio")]] = result[[v_m]]/result$variance
    }
    write.csv(result, file = file_name, row.names = FALSE)
  }
}

######################################

methods = c("gmm", "reg")
se_methods = c("delta_v_ratio", "semi_par_v_ratio", "mn_boot_1_v_ratio" ,
               "mn_boot_2sqrt_v_ratio", "mn_boot_23_v_ratio")
prob_types = c("tp", "tn")


for(method in methods)
{
  se_ratio_table = data.frame(data.frame(matrix(0, ncol=10, nrow= 32)))
  colnames(se_ratio_table) = paste(rep(prob_types,each=5), se_methods, sep = "_")
  for(prob_type in prob_types)
  {
    result = read.csv(paste0("./results/se/", method,
                             "_", prob_type,"_est_summary.csv"))
    for(se_method in se_methods)
    {
      se_ratio_table[[paste0(prob_type,"_",se_method)]] = result[[se_method]][correct_ind]
    }
  }
  write.csv(se_ratio_table,
            file=paste0("./results/se/v_ratio_table_",method,".csv"),
            row.names = FALSE)
}




