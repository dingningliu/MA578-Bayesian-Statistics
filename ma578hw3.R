#3.12
#read data
fatal <- data.frame(year = 1:10, accidents = c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22) )


#(e)
# check frequentist estimates
linear_model <- lm(accidents ~ year, data = fatal)
summary(linear_model)
cat('Alpha: ', coef(linear_model)[1])
cat('Beta: ', coef(linear_model)[2])
cat('Standard Error of Alpha: ', summary(linear_model)$coefficients[1,2])
cat('Standard Error of Beta: ', summary(linear_model)$coefficients[2,2])

#(f)
# P(alpha, beta) propto 1
#Sketch its contours 
m <- 100 # number of grid subdivisions
alpha <- seq(20, 40, length=m) 
beta <- seq(-2, 0.5, length=m) 


#loglikelihood
a <- coef(linear_model)[1]
b <- coef(linear_model)[2]
lhood <- function(a, b) {     #log-likelihood
  sum(-(a+b*fatal$year) + fatal$accidents * log(a+b*fatal$year))
}

#log-prior
lprior = log(1)

#joint posterior
lj <- matrix(nrow=m, ncol=m)
for (i in 1:m) {
  for (j in 1:m) {
    lj[i, j] = lprior + lhood(alpha[i], beta[j])
  }
}
pab <- exp(lj - log(sum(exp(lj)))) 
# plot
#image(alpha, beta, pab)
#points(coef(linear_model)[1], coef(linear_model)[2], pch=3)
#contour(alpha, beta, pab, add=TRUE)
contour(alpha, beta, pab, xlim = c(20,40), ylim = c(-2.5,0.5), xlab='alpha', ylab='beta')


#(g)  
#sample alpha beta
pa <- rowSums(pab) # marginal on alpha
ns <- 1000 # #samples
alpha.s <- beta.s <- numeric(ns)
for (s in 1:ns) {
  ia <- sample.int(m, 1, prob=pa) # sample alpha
  ib <- sample.int(m, 1, prob=pab[ia,]) # sample beta | alpha
  alpha.s[s] <- alpha[ia]; beta.s[s] <- beta[ib]
}
hist(alpha.s+beta.s*11, main = expression(paste('Histogram of ', alpha+11*beta)), xlab = expression(alpha+11*beta))


#(h)
y_1986 <- rpois(1000, alpha.s+beta.s*11)
quantile(y_1986,c(0.025,0.975))
