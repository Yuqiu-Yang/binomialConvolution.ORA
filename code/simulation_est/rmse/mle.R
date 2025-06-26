#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
i_setting_s = as.numeric(args[1])
i_setting_e = as.numeric(args[2])
###########################################
# MLE
###########################################
set.seed(42)
datafolder = "data/simulation/rmse"
codefolder = "code"
source(paste0(codefolder, "/utility.R"))
source(paste0(codefolder, "/estimate_mom.R"))
source(paste0(codefolder, "./estimate_mle.R"))
setwd(datafolder)
simulation_setting = read.csv("./simulation_setting.csv")
n_simulation = 1000
significance_level = 0.05
for(i_setting in i_setting_s : i_setting_e)
{
  setting_folder = paste0("./setting_", i_setting)
  for(i_simulation in 1 : n_simulation)
  {
    if(paste0("passage_", i_simulation, "_mle_result.rds") %in%
       list.files(paste0("./setting_",i_setting, "/passage_est/")))
    {
      next
    }
    file_name = paste0(setting_folder, "/passage_data/passage_", i_simulation)
    passage_data = readRDS(paste0(file_name, ".rds"))

    passage_moments = get_passage_moments(passage_data=passage_data)
    mle_est = estimate_mle(passage_data=passage_data,
                           passage_moments=passage_moments,
                           significance_level=significance_level,
                           return_ci=FALSE)

    result = list("mle_est"=mle_est)
    file_name = paste0(setting_folder, "/passage_est/passage_", i_simulation)
    saveRDS(result, file=paste0(file_name, "_mle_result.rds"))
  }
}


