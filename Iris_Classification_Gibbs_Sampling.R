# install.packages("MCMCpack")
# install.packages("mvtnorm")
library(MCMCpack)
library(mvtnorm)

# Preprocessing
Y <- iris[, 1:4]; Y <- as.matrix(Y)
Z <- matrix(rep(0, 450), 150, 3)
cj <- factor(iris[, 5], labels = 1:3)
colnames(Z) <- c("setosa", "versicolor", "virginica")

# Preparing label matrix Z
# This part is optional to view real class labels
for(i in 1:150) {
  j <- cj[i]
  Z[i, j] <- 1
}

# Initialization of parameters
w <- matrix(rep(1/3, 450), 150, 3)
mu <- list(m_1 = as.matrix(colMeans(Y) + rnorm(4, 0, 1), 4, 1), 
           m_2 = as.matrix(colMeans(Y) + rnorm(4, 0, 1), 4, 1),
           m_3 = as.matrix(colMeans(Y) + rnorm(4, 0, 1), 4, 1))
mu_0 <- list(m_1 = as.matrix(colMeans(Y), 4, 1), 
             m_2 = as.matrix(colMeans(Y), 4, 1),
             m_3 = as.matrix(colMeans(Y), 4, 1))
sigma <- list(s_1 = cov(Y)/3, s_2 = cov(Y)/3, s_3 = cov(Y)/3)
sigma_0 <- list(s_1 = cov(Y)/3, s_2 = cov(Y)/3, s_3 = cov(Y)/3)

maxiter <- 2000 # Maximum number of iterations

w.result <- matrix(rep(1/3, 450), 150, 3)
mu.result <- list(m_1 = as.matrix(colMeans(Y) + rnorm(4, 0, 1), 4, 1), 
                  m_2 = as.matrix(colMeans(Y) + rnorm(4, 0, 1), 4, 1),
                  m_3 = as.matrix(colMeans(Y) + rnorm(4, 0, 1), 4, 1))
sigma.result <- list(s_1 = cov(Y)/3, s_2 = cov(Y)/3, s_3 = cov(Y)/3)

# Gibbs sampling main algorithm
for(t in 1:maxiter) {
  # Assigning cluster membership probabilities 
  # By using Mahalanobis distance between a data point and cluster means
  d <- matrix(rep(0, 450), 150, 3)
  for(i in 1:3) {
    d[, i] <-  sqrt(mahalanobis(x = Y, center = mu[[i]], cov = sigma[[i]]))
  }
  for(i in 1:150) {
    d[i, ] <- prod(d[i, ]) / d[i, ]
  }
  d <- d / rowSums(d)
  w <- w * d; w <- w / rowSums(w)
  I <- matrix(rep(0, 450), 150, 3)
  # Updating cluster probabilities by using Dirichlet-Categorical model
  for(i in 1:150) {
    theta <- rdirichlet(1, w[i, ]*150)
    I[i, ] <- as.numeric(theta == max(theta))
  }
  Cj <- colSums(I)
  # Updating the parameters of Gaussians
  sigma_hat <- list(s_1 = 0, s_2 = 0, s_3 = 0)
  mu_hat <- list(s_1 = 0, s_2 = 0, s_3 = 0)
  for(i in 1:3) {
    y <- Y - matrix(rep(1, 150), 150, 1) %*% t(mu[[i]])
    lambda <- 0
    for(j in 1:150) {
      lambda <- lambda + y[j, ] %*% t(y[j, ])
    }
    sigma[[i]] <- riwish(4 + Cj[i]*4, lambda)
    sigma_hat[[i]] <- solve((solve(sigma_0[[i]]) + Cj[i]*solve(sigma[[i]])))
    mu_hat[[i]] <- sigma_hat[[i]] %*% (t(t(mu_0[[i]]) %*% solve(sigma_0[[i]])) + solve(sigma[[i]]) %*% t(t(I[, i]) %*% Y))
    mu[[i]] <- t(rmvnorm(1, mu_hat[[i]], sigma_hat[[i]]))
  }
  # Accumulating the parameters
  # To be divided by number of iterations later on
  # In order to take their average
  if(t > 1) {
    w.result <- w.result + w
    for(i in 1:3) {
      mu.result[[i]] <- mu.result[[i]] + mu[[i]]
      sigma.result[[i]] <- sigma.result[[i]] + sigma[[i]]
    }
  }
}

# Dividing the sum of parameters along the path by number of iterations
# In order to obtain their estimate (average)
w.result <- w.result / maxiter 
for(i in 1:3) {
  mu.result[[i]] <- mu.result[[i]] / maxiter
  sigma.result[[i]] <- sigma.result[[i]] / maxiter 
}

# Test of results
# With respect to cluster probabilities in the last iteration
colSums(w[1:50,])
colSums(w[51:100,])
colSums(w[101:150,])
