###########################################
# REG
###########################################
rm(list=ls())
set.seed(42)
setwd("/work/DPDS/s205711/ORA/rmse_simulation/")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/utility.R")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/estimate_regression.R")
simulation_setting = read.csv("./simulation_setting.csv")
n_simulation = 1000
significance_level = 0.05
for(i_setting in 1 : nrow(simulation_setting))
{
  setting_folder = paste0("./setting_", i_setting)
  for(i_simulation in 1 : n_simulation)
  {
    file_name = paste0(setting_folder, "/passage_data/passage_", i_simulation)
    passage_data = readRDS(paste0(file_name, ".rds"))

    reg_est = estimate_linear_regression(passage_data=passage_data,
                                         significance_level=significance_level,
                                         return_ci=FALSE)

    result = list("reg_est"=reg_est)
    file_name = paste0(setting_folder, "/passage_est/passage_", i_simulation)
    saveRDS(result, file=paste0(file_name, "_reg_result.rds"))
  }
}
