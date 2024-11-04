setwd("/work/DPDS/s205711/ORA/simulation_misspecification/")
# setwd("/work/DPDS/s205711/ORA/simulation/")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/utility.R")
# We first generate simulation settings
# n_students = 50
# n_words = c(44, 69)
# positive_prob = c(0.96, 0.98)
# positive_icc = c(0, 0.01, 0.06)[c(1,3)]
# true_positive_prob = c(0.98, 0.999)[1]
# true_negative_prob = c(0.75, 0.85)[1]
# true_positive_icc = c(0, 0.01, 0.06)[c(2,3)]
# true_negative_icc = c(0, 0.01, 0.06)[c(2,3)]
#
# simulation_setting = expand.grid(n_students,
#                                 n_words,
#                                 positive_prob,
#                                 positive_icc,
#                                 true_positive_prob,
#                                 true_negative_prob,
#                                 true_positive_icc,
#                                 true_negative_icc)
# colnames(simulation_setting) = c("n_students",
#                                  "n_words",
#                                  "positive_prob",
#                                  "positive_icc",
#                                  "true_positive_prob",
#                                  "true_negative_prob",
#                                  "true_positive_icc",
#                                  "true_negative_icc")
#
# write.csv(simulation_setting,
#           file="./simulation_setting.csv", row.names = FALSE)

set.seed(42)
###########################################
# We first generate data
###########################################
simulation_setting = read.csv("./simulation_setting.csv")
n_simulation = 1000
for(i_setting in 1 : nrow(simulation_setting))
{
  setting_folder = paste0("./setting_", i_setting,"/passage_data")
  dir.create(setting_folder, recursive = T)
  for(i_simulation in 1 : n_simulation)
  {
    file_name = paste0(setting_folder, "/passage_", i_simulation)
    passage_data = simulate_passages(positive_prob=simulation_setting$positive_prob[i_setting],
                                     true_positive_prob=simulation_setting$true_positive_prob[i_setting],
                                     true_negative_prob=simulation_setting$true_negative_prob[i_setting],
                                     positive_icc=simulation_setting$positive_icc[i_setting],
                                     true_positive_icc=simulation_setting$true_positive_icc[i_setting],
                                     true_negative_icc=simulation_setting$true_negative_icc[i_setting],
                                     passage_name="1",
                                     n_students=simulation_setting$n_students[i_setting],
                                     n_words=simulation_setting$n_words[i_setting])
    saveRDS(passage_data,
            file=paste0(file_name, ".rds"))
  }
}

###########################################
# Then we generate estimates
###########################################

for(i_setting in 1 : nrow(simulation_setting))
{
  setting_folder = paste0("./setting_", i_setting,"/passage_est")
  dir.create(setting_folder, recursive = T)
}







