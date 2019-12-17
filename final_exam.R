library(bayesplot)
library(rstan)
library(beanplot)
stroke <- within(read.csv("/Users/liudingning/Desktop/MAcourses/MA 578/stroke.csv"), subject <- factor(subject))
source("/Users/liudingning/Desktop/MAcourses/MA 578/bslm.R")
# setup design matrices



tem <- stroke %>% group_by(week) %>% mutate(weekly_score = sum(score))
df <- data.frame(score=tem$weekly_score[1:8], week=c(1,2,3,4,5,6,7,8))
y <- df$score
X <- model.matrix(~ week, data = df)
n <- nrow(X)

# function

#jeffery prior
log_Ptheta <- function(theta) {
  log(((n-1)/((1-theta)^2)+(1/(theta^2))))/2
}

# sample from inverse (scaled) chi-square with parameters `nu` and `tau2`;
# nu_tau2 = nu * tau2 for convenience
rinvsquare <- function (ns, nu, nu_tau2) {1 / rgamma(ns, nu / 2, nu_tau2 / 2)}


#RSS
RSS <- function(y, x, beta) {
  crossprod(y-x%*%beta) 
}

#likelihood
lhood_ytheta <- function(y, x, beta, theta, sigma2){
  t <- (y-x%*%beta)/sqrt(sigma2)
  t_m <- sum(t)/n
  s_t <- crossprod(t- t_m)/n
  -log(theta)/2-(n-1)*log(1-theta)/2-((n-1)*s_t/(1-theta) + t_m^2/theta)/2
}


#rou
rou <- function(theta){
  (n*theta-1)/(n-1)
}

#r_theta
r_theta <- function(rou) {
  (1-rou)*diag(n) + rou*matrix(1, nrow = n, ncol = n)
}

#get mu
mvn <- function(mu, sigma) {
  C <- chol(sigma)
  z <- rnorm(length(mu), backsolve(C, mu, transpose=TRUE), rep(1, length(mu)))
  return(crossprod(C,z))
}

mcmc_array <- function (ns, nchains = 1, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

plot_trace <- function (x, ...)
  plot(x, type = "l", col = "gray", xlab = "iteration", ...)
plot_hist <- function (x, ...)
  hist(x, col = "gray", border = "white", main = "", ...)


# init
bl <- lm(y ~ X - 1)
beta <- coef(bl)
sigma2 <- deviance(bl) / n
theta <- 0.5
m <- 100
theta_0 <- seq(0.01, 0.98, length=m) 

ns <- 2000 # number of mcmc samples
params <- c("beta_intercept", "beta_week","sigma2", "rou","theta") 
sims <- mcmc_array(ns, params = params)   #simulation 1 chain




for (i in 1:ns) {
  rou_0 <- rou(theta)
  r_theta_0 <- r_theta(rou_0)
  C <- chol(r_theta_0)
  y_tilde <- backsolve(C, y, transpose=TRUE)
  X_tilde <- backsolve(C, X, transpose=TRUE)
  
  

  # [ beta | theta, sigma2, y_tilde ]
  #beta <- rnorm(1, solve(t(X_tilde)%*%X_tilde)%*%t(X_tilde)%*%y_tilde, sigma2*solve(t(X_tilde)%*%X_tilde))
  beta <- mvn(backsolve(crossprod(X_tilde), crossprod(X_tilde,y_tilde)),  sigma2*solve(crossprod(X_tilde)))
  
  
  # [ sigam2 | beta, theta, y_tilde ]
  RSS_tilde <- RSS(y_tilde, X_tilde, beta)
  sigma2 <- rinvsquare(1, n, RSS_tilde)
  
  # [ theta | beta, sigma2, y_tilde ]
  #iterate
  lj <- rep(0, m)
  for (j in 1:m) {
    c_0 <- chol(r_theta(rou(theta_0[j])))
    lj[j] <- lhood_ytheta(y, X, beta, theta_0[j],sigma2) + log_Ptheta(theta_0[j])
  }
  lj_t <- lj - max(lj)
  theta <- sample(theta_0, 1, prob = exp(lj_t - log(sum(exp(lj_t)))))
  rou_0 <- rou(theta)
  
  sims[i, 1, ] <- c(beta, sigma2, rou_0, theta)
}



r <- 1000:2000 # discard burn-in

monitor(sims)
par(mfrow = c(2, 2))
plot_trace(sims[r, , 1], ylab = 'beta 0')
plot_trace(sims[r, , 2], ylab = 'beta 1')
plot_trace(sqrt(sims[r, , 3]), ylab = expression(sigma))
plot_trace(sims[r, , 4], ylab = 'rho')

plot_hist(sims[r, , 1], ylab = 'beta 0')
plot_hist(sims[r, , 2], ylab = 'beta 1')
plot_hist(sqrt(sims[r, , 3]), ylab = expression(sigma))
plot_hist(sims[r, , 4], ylab = 'rho')

acf(sims[r, , 1], main = 'beta 0')
acf(sims[r, , 2], main = 'beta 1')
acf(sqrt(sims[r, , 3]), main = expression(sigma))
acf(sims[r, , 4], main = 'rho')

# comparison
library(lme4)
re <- lmer(score ~ (1 + week | subject), data = stroke)
summary(re)



#posterior predictive checks
p <- ncol(X)
y_rep <- matrix(nrow = ns, ncol = n)
y_rep_0 <- matrix(nrow = ns, ncol = n)
for (is in 1:ns){
  y_rep_0[is,] <- rnorm(n, X %*% sims[is, 1, 1:p], sqrt(sims[is, 1, p + 1]))
  r_theta_0 <- r_theta(sims[is, 1, p + 2])
  c <- chol(r_theta_0)
  y_rep[is,] <- crossprod(c, y_rep_0[is,])
}
par(mfrow=c(1,1))
boxplot(y_rep, outline = F); points(y, pch = 19, col = "red")


#outlier analysis
e_i <- matrix(nrow = ns, ncol = n)
for (is in 1:ns) {
  e_i[is,] <- (y - y_rep[is,])/sqrt(sims[is,1,p+1])  }

plot(e_i[r], col='black', xlab='Index', ylab="Error")
abline(h=2, lty=2,col = "red"); abline(h=-2, lty=2,col = "red")
hist(e_i[r], main='Histogram of Error', xlab = 'Error')   


alpha <- 0.05
k_alpha <- -qnorm(alpha/2)
mean(abs(e_i[r])<k_alpha) #prob of not-outlier observations 


#(d)

temp1<-data.frame(nscore=c(320, 335, 370, 420, 495, 535, 595, 650, 295, 320, 370, 410, 455, 495, 510, 515, 245, 310, 330, 365, 390, 410, 445, 460), week=rep(1:8,3), group= c(rep('A',8),rep('B',8),rep('C',8)))

# setup design matrices
y_d <- temp1$nscore
X_0 <- model.matrix(~ week + group, data = temp1)
X_1 <- model.matrix(~ week + group + week:group, data = temp1)





n <- nrow(X_0)
bl_0 <- lm(y_d ~ X_0 - 1)
beta <- coef(bl_0)
sigma2 <- deviance(bl) / n
theta <- 0.5
m <- 100
theta_0 <- seq(0.01, 0.98, length=m) 

ns <- 2000 # number of mcmc samples
params_0 <- c("beta_intercept", "beta_week","beta_B","beta_C","sigma2", "rou","theta") 
sims_0 <- mcmc_array(ns, params = params_0)   #simulation 1 chain




for (i in 1:ns) {
  rou_0 <- rou(theta)
  r_theta_0 <- r_theta(rou_0)
  C <- chol(r_theta_0)
  y_tilde <- backsolve(C, y_d, transpose=TRUE)
  X_tilde <- backsolve(C, X_0, transpose=TRUE)
  
  
  
  # [ beta | theta, sigma2, y_tilde ]
  #beta <- rnorm(1, solve(t(X_tilde)%*%X_tilde)%*%t(X_tilde)%*%y_tilde, sigma2*solve(t(X_tilde)%*%X_tilde))
  beta <- mvn(backsolve(crossprod(X_tilde), crossprod(X_tilde,y_tilde)),  sigma2*solve(crossprod(X_tilde)))
  
  
  # [ sigam2 | beta, theta, y_tilde ]
  RSS_tilde <- RSS(y_tilde, X_tilde, beta)
  sigma2 <- rinvsquare(1, n, RSS_tilde)
  
  # [ theta | beta, sigma2, y_tilde ]
  #iterate
  lj <- rep(0, m)
  for (j in 1:m) {
    c_0 <- chol(r_theta(rou(theta_0[j])))
    lj[j] <- lhood_ytheta(y_d, X_0, beta, theta_0[j],sigma2) + log_Ptheta(theta_0[j])
  }
  lj_t <- lj - max(lj)
  theta <- sample(theta_0, 1, prob = exp(lj_t - log(sum(exp(lj_t)))))
  rou_0 <- rou(theta)
  
  sims_0[i, 1, ] <- c(beta, sigma2, rou_0, theta)
}



r <- 1000:2000 # discard burn-in

monitor(sims_0)
par(mfrow = c(2, 3))
plot_trace(sims_0[r, , 1], ylab = 'beta 0')
plot_trace(sims_0[r, , 2], ylab = 'beta 1')
plot_trace(sims_0[r, , 3], ylab = 'beta B')
plot_trace(sims_0[r, , 4], ylab = 'beta C')
plot_trace(sqrt(sims_0[r, , 5]), ylab = expression(sigma))
plot_trace(sims_0[r, , 6], ylab = 'rho')

plot_hist(sims_0[r, , 1], ylab = 'beta 0')
plot_hist(sims_0[r, , 2], ylab = 'beta 1')
plot_hist(sims_0[r, , 3], ylab = 'beta B')
plot_hist(sims_0[r, , 4], ylab = 'beta C')
plot_hist(sqrt(sims_0[r, , 5]), ylab = expression(sigma))
plot_hist(sims_0[r, , 6], ylab = 'rho')



n <- nrow(X_1)
bl_1 <- lm(y_d ~ X_1 - 1)
beta <- coef(bl_1)
sigma2 <- deviance(bl) / n
theta <- 0.5
m <- 100
theta_0 <- seq(0.01, 0.98, length=m) 

ns <- 2000 # number of mcmc samples
params_1 <- c("beta_intercept", "beta_week","beta_B","beta_C","beta_week_B","beta_week_C","sigma2", "rou","theta") 
sims_1 <- mcmc_array(ns, params = params_1)   #simulation 1 chain




for (i in 1:ns) {
  rou_0 <- rou(theta)
  r_theta_0 <- r_theta(rou_0)
  C <- chol(r_theta_0)
  y_tilde <- backsolve(C, y_d, transpose=TRUE)
  X_tilde <- backsolve(C, X_1, transpose=TRUE)
  
  
  
  # [ beta | theta, sigma2, y_tilde ]
  #beta <- rnorm(1, solve(t(X_tilde)%*%X_tilde)%*%t(X_tilde)%*%y_tilde, sigma2*solve(t(X_tilde)%*%X_tilde))
  beta <- mvn(backsolve(crossprod(X_tilde), crossprod(X_tilde,y_tilde)),  sigma2*solve(crossprod(X_tilde)))
  
  
  # [ sigam2 | beta, theta, y_tilde ]
  RSS_tilde <- RSS(y_tilde, X_tilde, beta)
  sigma2 <- rinvsquare(1, n, RSS_tilde)
  
  # [ theta | beta, sigma2, y_tilde ]
  #iterate
  lj <- rep(0, m)
  for (j in 1:m) {
    c_0 <- chol(r_theta(rou(theta_0[j])))
    lj[j] <- lhood_ytheta(y_d, X_1, beta, theta_0[j],sigma2) + log_Ptheta(theta_0[j])
  }
  lj_t <- lj - max(lj)
  theta <- sample(theta_0, 1, prob = exp(lj_t - log(sum(exp(lj_t)))))
  rou_0 <- rou(theta)
  
  sims_1[i, 1, ] <- c(beta, sigma2, rou_0, theta)
}



monitor(sims_0)
monitor(sims_1)

