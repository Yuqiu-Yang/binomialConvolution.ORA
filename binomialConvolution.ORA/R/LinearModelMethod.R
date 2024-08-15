setwd("/Users/cpotgieter/Library/CloudStorage/Box-Box/Fall 2023/Stat Model Estimation")
load("passage_list.RData")

# Initialize empty storage frame
data <- data.frame()

# Loop over Passage_01 to Passage_10
for (i in 1:10) {

  passage_name <- paste0("Passage_", sprintf("%02d", i))
  
  data_extract <- passage_list[[passage_name]]
  
  temp_data <- data.frame(Passage = rep(i, length(data_extract$X)),
                          X = data_extract$X,
                          Y.Human = data_extract$Y.Human,
                          Y.AI = data_extract$Y.AI,
                          N = data_extract$N)

  data <- rbind(data, temp_data)
}

#View(data)

# Individual passage estimates (Human)

results <- data.frame(Passage = integer(),
                      pi.tp = numeric(),
                      pi.tn = numeric(),
                      se.tp = numeric(),
                      se.tn = numeric())

# Loop over each passage 
for (i in 1:10) {
  # Subset the data for the current passage
  data_red <- data[data$Passage == i,]
  Y <- data_red$Y.Human
  X1 <- data_red$X
  X2 <- data_red$N - data_red$X
  
  # Fit a regression model without an intercept
  mod <- lm(Y ~ 0 + X1 + X2)
  
  # Extract pi.hat and se values
  pi.hat <- c(mod$coefficients[1], 1 - mod$coefficients[2])
  pi.se <- as.numeric(summary(mod)$coefficients[, 2])
  
  # Combine the passage number, pi.hat, and se values into a single row
  result_row <- data.frame(Passage = i,
                           pi.tp = pi.hat[1],
                           pi.tn = pi.hat[2],
                           se.tp = pi.se[1],
                           se.tn = pi.se[2])
  
  # Append the result_row to the results data frame
  results <- rbind(results, result_row)
}

# Combined LS estimators -- Human
Y <- data$Y.Human
X1 <- data$X
X2 <- data$N - data$X
mod_combined <- lm(Y~0+X1+X2)
pi.hat <- c(mod_combined$coefficients[1],1-mod_combined$coefficients[2])
pi.se <- as.numeric(summary(mod_combined)$coefficients[,2])

result_row <- data.frame(Passage = "All",
                             pi.tp = pi.hat[1],
                             pi.tn = pi.hat[2],
                             se.tp = pi.se[1],
                             se.tn = pi.se[2])
results <- rbind(results, result_row)
row.names(results) <- NULL
results.human <- results

# Individual passage estimates (AI)

results <- data.frame(Passage = integer(),
                      pi.tp = numeric(),
                      pi.tn = numeric(),
                      se.tp = numeric(),
                      se.tn = numeric())

# Loop over each passage 
for (i in 1:10) {
  # Subset the data for the current passage
  data_red <- data[data$Passage == i,]
  Y <- data_red$Y.AI
  X1 <- data_red$X
  X2 <- data_red$N - data_red$X
  
  # Fit a regression model without an intercept
  mod <- lm(Y ~ 0 + X1 + X2)
  
  # Extract pi.hat and se values
  pi.hat <- c(mod$coefficients[1], 1 - mod$coefficients[2])
  pi.se <- as.numeric(summary(mod)$coefficients[, 2])
  
  # Combine the passage number, pi.hat, and se values into a single row
  result_row <- data.frame(Passage = i,
                           pi.tp = pi.hat[1],
                           pi.tn = pi.hat[2],
                           se.tp = pi.se[1],
                           se.tn = pi.se[2])
  
  # Append the result_row to the results data frame
  results <- rbind(results, result_row)
}

#row.names(results) <- NULL
#results

# Combined LS estimators -- Human
Y <- data$Y.AI
X1 <- data$X
X2 <- data$N - data$X
mod_combined <- lm(Y~0+X1+X2)
pi.hat <- c(mod_combined$coefficients[1],1-mod_combined$coefficients[2])
pi.se <- as.numeric(summary(mod_combined)$coefficients[,2])

result_row <- data.frame(Passage = "All",
                         pi.tp = pi.hat[1],
                         pi.tn = pi.hat[2],
                         se.tp = pi.se[1],
                         se.tn = pi.se[2])
results <- rbind(results, result_row)
row.names(results) <- NULL
results.ai <- results

results.human
results.ai
