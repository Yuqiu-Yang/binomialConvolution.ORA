library(dplyr)
library(extraDistr)

#'
#' @param n_students Number of students. This is the number of observations you want to generate
#' @param size
#' @param prob
#' @param overdispersion
#' @return
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
    icc = overdispersion/(size-1)
    par_alpha = prob * (1-icc) / icc
    par_list$fun_args$alpha = par_alpha
    par_list$fun_args$beta = par_alpha * (1/prob - 1)
  }else{
    par_list$fun_name = "rbinom"
    par_list$fun_args$prob = prob
  }
  return(par_list)
}


#'
#' @param positive_prob The probability of correct prounciation
#' @param true_positive_prob
#' @param true_negative_prob
#' @param positive_overdispersion
#' @param true_positive_overdispersion
#' @param true_negative_overdispersion
#' @param passage_name
#' @param n_students
#' @param n_words
#' @return description
#' @export
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


simulate_passages <- function(positive_prob,
                              true_positive_prob,
                              true_negative_prob,
                              positive_overdispersion=0,
                              true_positive_overdispersion=0,
                              true_negative_overdispersion=0,
                              passage_name=as.character(1:2),
                              n_students=rep(40,2),
                              n_words=rep(50,2))
{
  P = length(passage_name)
  passage_data = c()
  for(p in 1 : P)
  {
    result = generate_counts(positive_prob=positive_prob,
                             true_positive_prob=true_positive_prob,
                             true_negative_prob=true_negative_prob,
                             positive_overdispersion=positive_overdispersion,
                             true_positive_overdispersion=true_positive_overdispersion,
                             true_negative_overdispersion=true_negative_overdispersion,
                             passage_name=passage_name[p],
                             n_students=n_students[p],
                             n_words=n_words[p])
    passage_data = rbind(passage_data, result)
  }
  passage_data$passage = factor(passage_data$passage,
                                levels=passage_name)

  return(passage_data)
}




resample_passage <- function(passage_data,
                              true_positive_prob=NA,
                              true_negative_prob=NA,
                              passage_name="1",
                              sample_prob=NA)
{
  passage_data = passage_data[which(passage_data$passage == passage_name), ]
  n_students = nrow(passage_data)
  n_words = passage_data$N[1]

  if(is.na(true_positive_prob) | is.na(true_negative_prob))
  {
    method = "m-out-of-n"
    m=floor(n_students * sample_prob)
  }else{
    method = "semi-parametric"
    m=n_students
  }

  ind=sample(n_students, size=m, replace=TRUE)

  n_positive=passage_data$X[ind]
  if(method=="semi-parametric")
  {
    true_positive_par_list = generate_counts_par(n_students=n_students,
                                                 size=n_positive,
                                                 prob=true_positive_prob,
                                                 overdispersion=0)
    n_true_positive = do.call(what=true_positive_par_list$fun_name,
                              args=true_positive_par_list$fun_args)

    false_positive_par_list = generate_counts_par(n_students=n_students,
                                                  size=n_words-n_positive,
                                                  prob=1-true_negative_prob,
                                                  overdispersion=0)
    n_false_positive = do.call(what=false_positive_par_list$fun_name,
                               args=false_positive_par_list$fun_args)

    observation = n_true_positive+n_false_positive
  }else{
    observation = passage_data$Y[ind]
  }

  result = data.frame("passage"=passage_name,
                      "X"=n_positive,
                      "Y"=observation,
                      "N"=n_words,
                      "n_positive"=n_positive,
                      "n_negative"=n_words-n_positive)

  return(result)
}


resample_passages <- function(passage_data,
                               true_positive_prob=NA,
                               true_negative_prob=NA,
                               sample_prob=NA)
{
  passage_name = unique(passage_data$passage)
  result = c()
  for(p_name in passage_name)
  {
    p_data = resample_passage(passage_data=passage_data,
                               true_positive_prob=true_positive_prob,
                               true_negative_prob=true_negative_prob,
                               passage_name=p_name,
                               sample_prob=sample_prob)
    result = rbind(result, p_data)
  }
  return(result)
}



bootstrap_passages <- function(passage_data,
                               true_positive_prob=NA,
                               true_negative_prob=NA,
                               n_bootstrap=50,
                               sample_prob=NA)
{
  result_list = list()
  for(b in 1 : n_bootstrap)
  {
    result_list[[b]] = resample_passages(passage_data=passage_data,
                                         true_positive_prob=true_positive_prob,
                                         true_negative_prob=true_negative_prob,
                                         sample_prob=sample_prob)
  }
  return(result_list)
}



get_passage_moments <- function(passage_data)
{
  # For each passage we compute the passage length and sample size
  Np = aggregate(passage_data$N, by=list(passage=passage_data$passage), FUN=function(x)x[1])$x
  np = aggregate(passage_data$passage, by=list(passag=passage_data$passage), FUN=function(x)length(x))$x
  # based on this we can compute the weights
  wp = Np*np/(sum(Np*np))

  # Then we compute for each passage the sample mean variance covariance
  mXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=mean)$x
  vXp = aggregate(passage_data$X, by=list(passage=passage_data$passage), FUN=var)$x
  mYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=mean)$x
  vYp = aggregate(passage_data$Y, by=list(passage=passage_data$passage), FUN=var)$x
  covp = passage_data %>% group_by(passage) %>% summarise(sig=cov(X,Y))
  cXYp = covp$sig

  mX = sum(mXp * wp)
  vX = sum(vXp * wp)
  cXY = sum(cXYp * wp)
  N = sum(Np * wp)
  mY = sum(mYp * wp)
  vY = sum(vYp * wp)

  return(list("Np"=Np,
              "np"=np,
              "wp"=wp,
              "mXp"=mXp,
              "vXp"=vXp,
              "mYp"=mYp,
              "vYp"=vYp,
              "cXYp"=cXYp,
              "mX"=mX,
              "vX"=vX,
              "cXY"=cXY,
              "N"=N,
              "mY"=mY,
              "vY"=vY))
}





