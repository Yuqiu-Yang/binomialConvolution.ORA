

curate_data <- function(passage_list)
{
  ora_list = list("human_ORA"=NA,
                  "ai_ORA"=NA)
  passage_names = sort(names(passage_list))

  for(measurement_type in c("human_ORA", "ai_ORA"))
  {
    temp_df = data.frame()
    for(passage_name in passage_names)
    {
      passage_data = passage_list[[passage_name]]

      temp = data.frame("passage"=passage_name,
                        "true_ORA"=passage_data$true_ORA,
                        "ORA"=passage_data[[measurement_type]],
                        "n_words"=passage_data$n_words,
                        "n_positives"=passage_data$true_ORA,
                        "n_negatives"=passage_data$n_words-passage_data$true_ORA)

      temp_df = rbind(temp_df, temp)
    }
    ora_list[[measurement_type]] = temp_df
  }
  return(ora_list)
}



