setwd("/work/DPDS/s205711/ORA/rmse_simulation/")
source("../binomialConvolution.ORA/binomialConvolution.ORA/R/utility.R")
# We first generate simulation settings
n_students = seq(from=30, to=100, by=5)
n_words = 60
positive_prob = 0.95
positive_icc = seq(from=0, to=0.06, length.out=15)
true_positive_prob = seq(from=0.85, to=0.999, length.out=15)
true_negative_prob = seq(from=0.5, to=0.95, length.out=15)

simulation_setting = NULL
# Number of students
temp = data.frame(n_students = n_students,
                  n_words = n_words,
                  positive_prob = positive_prob,
                  positive_icc = 0,
                  true_positive_prob = 0.98,
                  true_negative_prob = 0.7)
simulation_setting = rbind(simulation_setting, temp)
# Positive ICC
temp = data.frame(n_students = 50,
                  n_words = n_words,
                  positive_prob = positive_prob,
                  positive_icc = positive_icc,
                  true_positive_prob = 0.98,
                  true_negative_prob = 0.7)
simulation_setting = rbind(simulation_setting, temp)
# True positive prob
temp = data.frame(n_students = 50,
                  n_words = n_words,
                  positive_prob = positive_prob,
                  positive_icc = 0,
                  true_positive_prob = true_positive_prob,
                  true_negative_prob = 0.7)
simulation_setting = rbind(simulation_setting, temp)
# True negative prob
temp = data.frame(n_students = 50,
                  n_words = n_words,
                  positive_prob = positive_prob,
                  positive_icc = 0,
                  true_positive_prob = 0.98,
                  true_negative_prob = true_negative_prob)
simulation_setting = rbind(simulation_setting, temp)

write.csv(simulation_setting,
          file="./simulation_setting.csv", row.names = FALSE)

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
    passage_data = simulate_passages(positive_prob=0.95,
                                     true_positive_prob=simulation_setting$true_positive_prob[i_setting],
                                     true_negative_prob=simulation_setting$true_negative_prob[i_setting],
                                     positive_icc=simulation_setting$positive_icc[i_setting],
                                     true_positive_icc=0,
                                     true_negative_icc=0,
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







