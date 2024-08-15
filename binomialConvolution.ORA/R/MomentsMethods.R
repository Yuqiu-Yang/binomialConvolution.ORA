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

MOM_crude <- function(X,Y,N,n) {
  pi.tp.hat <- mean(Y)/N + cov(X,Y)/var(X)*(1-mean(X)/N)
  pi.tn.hat <- (1-mean(Y)/N) + cov(X,Y)/var(X)*mean(X)/N
  pi.hat <- c(pi.tp.hat, pi.tn.hat)
  
  mX <- mean(X)
  mY <- mean(Y)
  vX <- var(X)
  cXY <- cov(X,Y)
  data.omega <- cbind(X,Y,(X-mX)^2,(X-mX)*(Y-mY))
  Omega <- cov(data.omega)
  
  A <- rbind(
    c(-1/N*cXY/vX, 1/N, -cXY/vX^2*(1-mX/N),1/vX*(1-mX/N)),
    c(1/N*cXY/vX, -1/N, -cXY/vX^2*mX/N,1/vX*mX/N)
  )
  
  CV <- A %*% Omega %*% t(A) / n
  
  SE <- sqrt(diag(CV))
  
  return(list(pi.hat=pi.hat,SE=SE))
}

GMM_dist <- function(par,data,Emp.Mom,invOmega) {
  
  # data must be structured as a list with array items X and Y
  # and scalar items N and n
  
  mu.X <- par[1]
  var.X <- par[2]
  pi.tp <- par[3]
  pi.tn <- par[4]
  
  mu.Y <- mu.X * pi.tp + (N - mu.X) * (1 - pi.tn)
  var.Y <- mu.X * pi.tp * (1 - pi.tp) + (N - mu.X) * pi.tn * (1 - pi.tn) +
    var.X * (pi.tp + pi.tn - 1)^2
  cov.XY <- var.X * (pi.tp + pi.tn - 1)
  
  X <- data$X
  Y <- data$Y
  N <- data$N
  n <- data$n
  
  if (missing(Emp.Mom)) {
    Emp.Mom <- c(mean(X),var(X),mean(Y),var(Y),cov(X,Y))
    #Emp.Mom <- c(mean(X),var(X),mean(Y),cov(X,Y))
  }
  
  if (missing(invOmega)) {
    mX <- mean(X)
    mY <- mean(Y)
    n <- length(X)
    data.omega <- cbind(X,(X-mX)^2,Y,(Y-mY)^2,(X-mX)*(Y-mY))
    Omega <- cov(data.omega)
    invOmega <- solve(Omega)
  }
  
  Th.Mom <- c(mu.X, var.X, mu.Y, var.Y, cov.XY)
  
  D <- n * as.numeric( t(Emp.Mom - Th.Mom) %*% invOmega %*% (Emp.Mom - Th.Mom) )
  
  return(D)
  
}

## Simulated data

n <- 40
N <- 50
X <- rbinom(n, N, 0.95)
TP <- rbinom(n, X, 0.99)
TN <- rbinom(n, N-X, 1-0.7)
Y <- TP + TN

mom.ests <- MOM_crude(X,Y,N,n)

data.sub <- list(X = X,
                 Y = Y,
                 N = N,
                 n = n)

gmm.eval <- optim(c(mean(X),var(X),mom.ests$pi.hat), 
                  fn = GMM_dist, data = data.sub,
                  method = "L-BFGS-B", 
                  lower = c(0,0,0,0), 
                  upper = c(data.sub$N,Inf,1,1),
                  hessian = TRUE)


gmm.eval$par[c(3,4)]
sqrt(diag(solve(gmm.eval$hessian)))[c(3,4)]


index <- which(data$Passage==10)
data.sub <- list(X = data$X[index],
                 Y = data$Y.AI[index],
                 N = data$N[index[1]],
                 n = length(index))


mom.ests <- MOM_crude(data.sub$X, data.sub$Y, data.sub$N, data.sub$n)

gmm.eval <- optim(c(mean(data.sub$X),var(data.sub$X),mom.ests$pi.hat), 
                  fn = GMM_dist, data = data.sub,
                  method = "L-BFGS-B", 
                  lower = c(0,0,0,0), 
                  upper = c(data.sub$N,Inf,1,1),
                  hessian = TRUE)

gmm.eval$value
gmm.eval$par[c(3,4)]
sqrt(diag(solve(gmm.eval$hessian)))[c(3,4)]
