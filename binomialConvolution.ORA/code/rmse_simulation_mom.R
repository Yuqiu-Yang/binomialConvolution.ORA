###########################################
# MOM
###########################################
rm(list=ls())
set.seed(42)
setwd("/work/DPDS/s205711/ORA/rmse_simulation/")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/utility.R")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/estimate_mom.R")
simulation_setting = read.csv("./simulation_setting.csv")
n_simulation = 1000
significance_level = 0.05
for(i_setting in 1 : nrow(simulation_setting))
{
  print(i_setting)
  setting_folder = paste0("./setting_", i_setting)
  for(i_simulation in 1 : n_simulation)
  {
    file_name = paste0(setting_folder, "/passage_data/passage_", i_simulation)
    passage_data = readRDS(paste0(file_name, ".rds"))

    passage_moments = get_passage_moments(passage_data=passage_data)

    mom_est = estimate_mom(passage_data=passage_data,
                           passage_moments=passage_moments,
                           significance_level=significance_level,
                           return_ci=FALSE)

    result = list("mom_est"=mom_est)
    file_name = paste0(setting_folder, "/passage_est/passage_", i_simulation)
    saveRDS(result, file=paste0(file_name, "_mom_result.rds"))
  }
}
