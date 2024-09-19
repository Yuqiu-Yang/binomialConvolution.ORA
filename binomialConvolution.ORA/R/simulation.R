setwd("/Users/yuqiuianyang/Desktop/work/potgieter/binomialConvolution.ORA/data/simulation/")

# We first generate simulation settings
# n_students = 50
# n_words = c(44, 69)
# positive_prob = c(0.96, 0.98)
# positive_overdispersion = c(0, 0.01, 0.06)
# true_positive_prob = c(0.98, 0.999)
# true_negative_prob = c(0.75, 0.85)
# true_positive_overdispersion = 0
# true_negative_overdispersion = 0
#
# simulation_setting = expand.grid(n_students,
#                                 n_words,
#                                 positive_prob,
#                                 positive_overdispersion,
#                                 true_positive_prob,
#                                 true_negative_prob,
#                                 true_positive_overdispersion,
#                                 true_negative_overdispersion)
# colnames(simulation_setting) = c("n_students",
#                                  "n_words",
#                                  "positive_prob",
#                                  "positive_overdispersion",
#                                  "true_positive_prob",
#                                  "true_negative_prob",
#                                  "true_positive_overdispersion",
#                                  "true_negative_overdispersion")
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
                                     positive_overdispersion=simulation_setting$positive_overdispersion[i_setting],
                                     true_positive_overdispersion=simulation_setting$true_positive_overdispersion[i_setting],
                                     true_negative_overdispersion=simulation_setting$true_negative_overdispersion[i_setting],
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







