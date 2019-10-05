# install.packages("MCMCpack")
library(MCMCpack)

# Preprocessing
Y <- iris[, 1:4]; Y <- as.matrix(Y)
Z <- matrix(rep(0, 450), 150, 3)
Cj <- factor(iris[, 5], labels = 1:3)
colnames(Z) <- c("setosa", "versicolor", "virginica")

# Preparing label matrix Z
# This part is optional to view real class labels
for(i in 1:150) {
  j <- Cj[i]
  Z[i, j] <- 1
}

# Initialization of parameters
w <- matrix(rep(1/3, 450), 150, 3)
mu <- list(m_1 = rnorm(4, colMeans(Y), 1), 
           m_2 = rnorm(4, colMeans(Y), 1),
           m_3 = rnorm(4, colMeans(Y), 1))
sigma <- list(s_1 = cov(Y)/3, s_2 = cov(Y)/3, s_3 = cov(Y)/3)

maxiter = 2000

# Expectation-Maximization main algorithm
for(t in 1:maxiter) {
  # Assigning cluster membership probabilities 
  # By using Mahalanobis distance between a data point and cluster means
  d <- matrix(rep(450), 150, 3)
  for(i in 1:3) {
    d[, i] <-  sqrt(mahalanobis(x = Y, center = mu[[i]], cov = sigma[[i]]))
  }
  for(i in 1:150) {
    d[i, ] <- prod(d[i, ]) / d[i, ]
  }
  d <- d / rowSums(d)
  w <- w * d
  # Updating cluster probabilities by using Dirichlet-Categorical model
  for(i in 1:150) {
    w[i, ] <- rdirichlet(1, w[i, ]*d[i, ]*150)
  }
  w <- w / rowSums(w)
  # Updating the parameters of Gaussians
  for(i in 1:3) {
    mu[[i]] <- w[, i] %*%  Y / sum(w[, i])
    X <- Y - matrix(mu[[i]], 150, 4, byrow = T)
    z <- 0
    for(j in 1:150) {
      z <- z + (w[j, i] * X[j, ] %*% t(X[j, ]))
    }
    sigma[[i]] <- z / sum(w[, i])
  }
}

# Test of results
colSums(w[1:50,])
colSums(w[51:100,])
colSums(w[101:150,])
