
n_students = 50
n_words = c(44, 69)
positive_prob = c(0.96, 0.98)
positive_overdispersion = c(0, 0.01, 0.06)
true_positive_prob = c(0.98, 0.999)
true_negative_prob = c(0.75, 0.85)
true_positive_overdispersion = 0
true_negative_overdispersion = 0

simulation_setting = expand.grid(n_students,
                                n_words,
                                positive_prob,
                                positive_overdispersion,
                                true_positive_prob,
                                true_negative_prob,
                                true_positive_overdispersion,
                                true_negative_overdispersion)
colnames(simulation_setting) = c("n_students",
                                 "n_words",
                                 "positive_prob",
                                 "positive_overdispersion",
                                 "true_positive_prob",
                                 "true_negative_prob",
                                 "true_positive_overdispersion",
                                 "true_negative_overdispersion")



n_simulation = 100

i_setting = 1

passage_data = simulate_passages(positive_prob=simulation_setting$positive_prob[i_setting],
                                 true_positive_prob=simulation_setting$true_positive_prob[i_setting],
                                 true_negative_prob=simulation_setting$true_negative_prob[i_setting],
                                 positive_overdispersion=simulation_setting$positive_overdispersion[i_setting],
                                 true_positive_overdispersion=simulation_setting$true_positive_overdispersion[i_setting],
                                 true_negative_overdispersion=simulation_setting$true_negative_overdispersion[i_setting],
                                 passage_name="1",
                                 n_students=simulation_setting$n_students[i_setting],
                                 n_words=simulation_setting$n_words[i_setting])

passage_moments = get_passage_moments(passage_data=passage_data)

gmm_est = estimate_gmm(passage_data=passage_data,
                       passage_moments=passage_moments,
                       significance_level=0.05)
reg_est = estimate_linear_regression(passage_data=passage_data,
                                     significance_level=0.05)
mle_est = estimate_mle(passage_data=passage_data,
                       significance_level=0.05)

