rm(list=ls())
setwd("data")
simulation_setting = read.csv("./simulation/rmse_misspec/simulation_setting.csv")

#####################
rmse_result = data.frame(matrix(0, nrow=nrow(simulation_setting),
                                ncol=6))
colnames(rmse_result) = paste(rep(c("gmm", "reg", "mle"), each=2), c("_tp","_tn") ,sep='')
rmse_result$varying_variable = rep(c("tp_icc","tn_icc", "tp_tn_icc"), each=15)
rmse_result$varying_values = c(seq(from=0, to=0.06, length.out=15),
                               seq(from=0, to=0.06, length.out=15),
                               seq(from=0, to=0.06, length.out=15))


for(i_setting in 1:nrow(simulation_setting))
{
  true_positive_prob = simulation_setting$true_positive_prob[i_setting]
  true_negative_prob = simulation_setting$true_negative_prob[i_setting]
  for(method in c("gmm", "reg", "mle"))
  {
    tp_est = tn_est = numeric(1000)
    for(i_simulation in 1 : 1000)
    {
      est = readRDS(paste0("./simulation/rmse_misspec/setting_",
                           i_setting, "/passage_est/passage_",
                           i_simulation,"_",method,"_result.rds"))
      tp_est[i_simulation] = est[[paste0(method,"_est")]]$pi.hat[1]
      tn_est[i_simulation] = est[[paste0(method,"_est")]]$pi.hat[2]
    }
    tp_rmse = sqrt(mean((tp_est-true_positive_prob)^2))
    tn_rmse = sqrt(mean((tn_est-true_negative_prob)^2))
    rmse_result[i_setting, paste0(method, "_tp")] = tp_rmse
    rmse_result[i_setting, paste0(method, "_tn")] = tn_rmse
  }
}

write.csv(rmse_result, file="./results/simulation/rmse_misspec/rmse_misspec_result.csv",
          row.names = FALSE)
#####################
bias_result = data.frame(matrix(0, nrow=nrow(simulation_setting),
                                ncol=6))
colnames(bias_result) = paste(rep(c("gmm", "reg", "mle"), each=2), c("_tp","_tn") ,sep='')
bias_result$varying_variable = rep(c("tp_icc","tn_icc", "tp_tn_icc"), each=15)
bias_result$varying_values = c(seq(from=0, to=0.06, length.out=15),
                               seq(from=0, to=0.06, length.out=15),
                               seq(from=0, to=0.06, length.out=15))


for(i_setting in 1:nrow(simulation_setting))
{
  true_positive_prob = simulation_setting$true_positive_prob[i_setting]
  true_negative_prob = simulation_setting$true_negative_prob[i_setting]
  for(method in c("gmm", "reg", "mle"))
  {
    tp_est = tn_est = numeric(1000)
    for(i_simulation in 1 : 1000)
    {
      est = readRDS(paste0("./simulation/rmse_misspec/setting_",
                           i_setting, "/passage_est/passage_",
                           i_simulation,"_",method,"_result.rds"))
      tp_est[i_simulation] = est[[paste0(method,"_est")]]$pi.hat[1]
      tn_est[i_simulation] = est[[paste0(method,"_est")]]$pi.hat[2]
    }
    tp_rmse = abs(mean(tp_est-true_positive_prob))
    tn_rmse = abs(mean(tn_est-true_negative_prob))
    bias_result[i_setting, paste0(method, "_tp")] = tp_rmse
    bias_result[i_setting, paste0(method, "_tn")] = tn_rmse
  }
}

write.csv(bias_result, file="./results/simulation/rmse_misspec/bias_misspec_result.csv",
          row.names = FALSE)
#####################
var_result = data.frame(matrix(0, nrow=nrow(simulation_setting),
                                ncol=6))
colnames(var_result) = paste(rep(c("gmm", "reg", "mle"), each=2), c("_tp","_tn") ,sep='')
var_result$varying_variable = rep(c("tp_icc","tn_icc", "tp_tn_icc"), each=15)
var_result$varying_values = c(seq(from=0, to=0.06, length.out=15),
                               seq(from=0, to=0.06, length.out=15),
                               seq(from=0, to=0.06, length.out=15))


for(i_setting in 1:nrow(simulation_setting))
{
  true_positive_prob = simulation_setting$true_positive_prob[i_setting]
  true_negative_prob = simulation_setting$true_negative_prob[i_setting]
  for(method in c("gmm", "reg", "mle"))
  {
    tp_est = tn_est = numeric(1000)
    for(i_simulation in 1 : 1000)
    {
      est = readRDS(paste0("./simulation/rmse_misspec/setting_",
                           i_setting, "/passage_est/passage_",
                           i_simulation,"_",method,"_result.rds"))
      tp_est[i_simulation] = est[[paste0(method,"_est")]]$pi.hat[1]
      tn_est[i_simulation] = est[[paste0(method,"_est")]]$pi.hat[2]
    }
    tp_rmse = var(tp_est)
    tn_rmse = var(tn_est)
    var_result[i_setting, paste0(method, "_tp")] = tp_rmse
    var_result[i_setting, paste0(method, "_tn")] = tn_rmse
  }
}

write.csv(var_result, file="./results/simulation/rmse_misspec/var_misspec_result.csv",
          row.names = FALSE)







