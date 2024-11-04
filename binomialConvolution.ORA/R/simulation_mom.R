###########################################
# MOM
###########################################
rm(list=ls())
set.seed(42)
# setwd("/work/DPDS/s205711/ORA/simulation/")
setwd("/work/DPDS/s205711/ORA/simulation_misspecification/")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/utility.R")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/estimate_mom.R")
simulation_setting = read.csv("./simulation_setting.csv")
n_simulation = 1000
significance_level = 0.05
n_bootstrap = 0
for(i_setting in 1 : nrow(simulation_setting))
{
  setting_folder = paste0("./setting_", i_setting)
  for(i_simulation in 1 : n_simulation)
  {
    file_name = paste0(setting_folder, "/passage_data/passage_", i_simulation)
    passage_data = readRDS(paste0(file_name, ".rds"))

    passage_moments = get_passage_moments(passage_data=passage_data)

    mom_est = estimate_mom(passage_data=passage_data,
                           passage_moments=passage_moments,
                           significance_level=significance_level,
                           return_ci=TRUE)
    # Bootstrap
    semi_par_boot = bootstrap_passages(passage_data=passage_data,
                                       true_positive_prob=mom_est$pi.hat[1],
                                       true_negative_prob=mom_est$pi.hat[2],
                                       n_bootstrap=n_bootstrap,
                                       sample_prob=NA)
    mn_boot_2sqrt = bootstrap_passages(passage_data=passage_data,
                                       true_positive_prob=NA,
                                       true_negative_prob=NA,
                                       n_bootstrap=n_bootstrap,
                                       sample_prob=2/sqrt(nrow(passage_data)))
    mn_boot_23 = bootstrap_passages(passage_data=passage_data,
                                    true_positive_prob=NA,
                                    true_negative_prob=NA,
                                    n_bootstrap=n_bootstrap,
                                    sample_prob=2/3)
    semi_par_est = mn_boot_2sqrt_est = mn_boot_23_est = matrix(0, nrow=n_bootstrap, ncol=2)
    if(n_bootstrap >= 1)
    {
      for(n_boot in 1 : n_bootstrap)
      {
        passage_moments_boot = get_passage_moments(passage_data=semi_par_boot[[n_boot]])

        mom_est_boot = estimate_mom(passage_data=semi_par_boot[[n_boot]],
                                    passage_moments=passage_moments_boot,
                                    significance_level=significance_level,
                                    return_ci=FALSE)
        semi_par_est[n_boot, ] = mom_est_boot$pi.hat
      }

      for(n_boot in 1 : n_bootstrap)
      {
        passage_moments_boot = get_passage_moments(passage_data=mn_boot_2sqrt[[n_boot]])

        mom_est_boot = estimate_mom(passage_data=mn_boot_2sqrt[[n_boot]],
                                    passage_moments=passage_moments_boot,
                                    significance_level=significance_level,
                                    return_ci=FALSE)
        mn_boot_2sqrt_est[n_boot, ] = mom_est_boot$pi.hat
      }

      for(n_boot in 1 : n_bootstrap)
      {
        passage_moments_boot = get_passage_moments(passage_data=mn_boot_23[[n_boot]])

        mom_est_boot = estimate_mom(passage_data=mn_boot_23[[n_boot]],
                                    passage_moments=passage_moments_boot,
                                    significance_level=significance_level,
                                    return_ci=FALSE)
        mn_boot_23_est[n_boot, ] = mom_est_boot$pi.hat
      }
    }

    result = list("mom_est"=mom_est,
                  "semi_par_est"=semi_par_est,
                  "mn_boot_2sqrt_est"=mn_boot_2sqrt_est,
                  "mn_boot_23_est"=mn_boot_23_est)
    file_name = paste0(setting_folder, "/passage_est/passage_", i_simulation)
    saveRDS(result, file=paste0(file_name, "_mom_result.rds"))
  }
}
