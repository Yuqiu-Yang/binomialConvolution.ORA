

estimate_linear_regression <- function(passage_data,
                                       significance_level=0.05)
{
  Y = passage_data$Y
  X1 = passage_data$n_positive
  X2 = passage_data$n_negative

  # Fit a regression model without an intercept
  mod = lm(Y ~ 0 + X1 + X2)
  # Extract pi.hat and se values
  prob_est = c(mod$coefficients[1], 1 - mod$coefficients[2])

  Z = matrix(cbind(X1,X2),ncol=2)
  true_positive_prob = pmin(prob_est[1],1)
  true_negative_prob = pmin(prob_est[2],1)
  alpha.hat = matrix(c(true_positive_prob*(1-true_positive_prob),true_negative_prob*(1-true_negative_prob)),nrow=2)
  SE = sqrt(diag(solve(t(Z)%*%Z) %*% t(Z) %*% diag(as.numeric(Z %*% alpha.hat)) %*% Z %*% solve(t(Z)%*%Z)))

  prob_est = c(true_positive_prob, true_negative_prob)
  ul = prob_est + qnorm(1-significance_level/2)*SE
  ll = prob_est - qnorm(1-significance_level/2)*SE
  ul = pmax(pmin(ul, 1),0)
  ll = pmax(pmin(ll, 1),0)
  return(list("pi.hat"=unname(prob_est),
              "pi.hat.ul"=unname(ul),
              "pi.hat.ll"=unname(ll)))
}

