library(extraDistr)


generate_counts_par <- function(n_students,
                                size,
                                prob,
                                overdispersion=0)
{
  par_list = list("fun_name"=NA,
                  "fun_args"=list("n"=n_students,
                                  "size"=size))
  if(overdispersion > 0)
  {
    par_list$fun_name = "rbbinom"
    par_alpha = prob * (1-overdispersion) / overdispersion
    par_list$fun_args$alpha = par_alpha
    par_list$fun_args$beta = par_alpha * (1/prob - 1)
  }else{
    par_list$fun_name = "rbinom"
    par_list$fun_args$prob = prob
  }
  return(par_list)
}


generate_counts <- function(positive_prob,
                            true_positive_prob,
                            true_negative_prob,
                            positive_overdispersion=0,
                            true_positive_overdispersion=0,
                            true_negative_overdispersion=0,
                            passage_name="1",
                            n_students=40,
                            n_words=50)
{
  positive_par_list = generate_counts_par(n_students=n_students,
                                          size=n_words,
                                          prob=positive_prob,
                                          overdispersion=positive_overdispersion)
  n_positive = do.call(what=positive_par_list$fun_name,
                       args=positive_par_list$fun_args)

  true_positive_par_list = generate_counts_par(n_students=n_students,
                                               size=n_positive,
                                               prob=true_positive_prob,
                                               overdispersion=true_positive_overdispersion)
  n_true_positive = do.call(what=true_positive_par_list$fun_name,
                            args=true_positive_par_list$fun_args)

  false_positive_par_list = generate_counts_par(n_students=n_students,
                                                size=n_words-n_positive,
                                                prob=1-true_negative_prob,
                                                overdispersion=true_negative_overdispersion)
  n_false_positive = do.call(what=false_positive_par_list$fun_name,
                             args=false_positive_par_list$fun_args)

  observation = n_true_positive+n_false_positive

  result = data.frame("passage"=passage_name,
                      "X"=n_positive,
                      "Y"=observation,
                      "N"=n_words,
                      "n_positive"=n_positive,
                      "n_negative"=n_words-n_positive)

  return(result)
}


