## Generate some data (ANCOVA-like)
## - y: an integer response
## - fact: a categorical factor with three levels
## - x: a continuous covariate
rm(list=ls(all=TRUE))
library(coda)
library(rjags)
library(runjags)
library(MASS)
library(ggplot2)

set.seed(666)
Nlevels <- 5
Nsamp <- 200
Ntot <- Nlevels*Nsamp
fact <- gl(Nlevels,Nsamp,Ntot)
x <- rnorm(Ntot,0,5)
Xmat<-model.matrix(~fact+x)
int <- 2
betas <- seq(from=0.3,to=0.8,length=Nlevels-1)
slope <- 0.05
l_mu <- Xmat %*% c(int,betas,slope)
mu <- exp(l_mu)
mu[mu<0|is.na(mu)]<-0
disp <- 5 #size/shape
prob <- disp/(disp+mu) #prob
y <- rnbinom(Ntot,disp,prob)
data <- data.frame(y=y,x=x,fact=fact)

##visualize the data (can skip this step)
ggplot(data,aes(y=y, x=x, colour=fact))+geom_point()+geom_smooth(method="lm")
X <- model.matrix(~fact+x)
data.list <- list(Y=y, X=X, n=length(y),ngroups=dim(Xmat)[2])

model <- "model {
  for (i in 1:n){
    Y[i] ~ dnegbin(p[i], r)
    p[i] <- r / (r + mu[i])
    log(mu[i]) <- eta[i] ## could be bounded on the logarithmic scale, if eta is failing to converge, e.g., max(-20, min(20, eta[i]))
    eta[i]<-inprod(beta,X[i,]) #matrix multiplication on the design matrix X
  }
  for (i in 1:ngroups) {
    beta[i] ~ dnorm(0,1e-6) # coefficient terms, including slope and intercept and factors
  }
  r ~ dgamma(0.001, 0.001) #dispersion parameter
}
"

model2 <- "model {
  for (i in 1:n){
Y[i] ~ dpois(lambda[i])
lambda[i] <- exp(eta[i]) * r ## could be bounded on the logarithmic scale, if eta is failing to converge, e.g., max(-20, min(20, eta[i]))
eta[i]<-inprod(beta,X[i,]) #matrix multiplication on the design matrix X
}
for (i in 1:ngroups) {
beta[i] ~ dnorm(0,1e-6) # coefficient terms, including slope and intercept and factors
}
r ~ dgamma(0.001, 0.001) #dispersion parameter
}
"

model.jags <- run.jags(model=model,inits=NULL, data=data.list, monitor=c("beta","r"),method="rjags",modules=c("glm","bugs"),
                   n.chain=2,sample=1000,adapt=1000, burnin=1000, thin=2)

model.jags2 <- run.jags(model=model2,inits=NULL, data=data.list, monitor=c("beta","lambda"),method="rjags",modules=c("glm","bugs"),
                       n.chain=2,sample=1000,adapt=1000, burnin=1000, thin=2)

plot(model.jags)
plot(model.jags2)
#autocorr.diag(as.mcmc(model.jags))
summary(model.jags)
summary(model.jags2)

r <- mean(as.matrix(model.jags$mcmc)[,"r"])
model.jags.mcmc <- as.matrix(model.jags$mcmc)[,1:dim(Xmat)[2]]

r2 <- mean(as.matrix(model.jags2$mcmc)[,"r"])
model.jags.mcmc2 <- as.matrix(model.jags2$mcmc)[,1:dim(Xmat)[2]]

apply(exp(model.jags.mcmc %*% t(model.matrix(~gl(Nlevels,1,Nlevels)+rep(0,Nlevels)))), 2, mean)
apply(exp(model.jags.mcmc %*% t(model.matrix(~gl(Nlevels,1,Nlevels)+rep(0,Nlevels)))), 2, quantile,probs=c(0.025,0.975))
apply(exp(model.jags.mcmc2 %*% t(model.matrix(~gl(Nlevels,1,Nlevels)+rep(0,Nlevels)))), 2, mean)

data.glmNB<-glm.nb(y~fact+x)

summary(data.glmNB)
c(exp(c(int,betas,slope)),disp)
c(exp(coef(data.glmNB)),data.glmNB$theta)
apply(cbind(exp(model.jags.mcmc),as.matrix(model.jags$mcmc)[,"r"]),2,mean)
apply(cbind(exp(model.jags.mcmc2),as.matrix(model.jags2$mcmc)[,"r"]),2,mean)

apply(cbind(exp(model.jags.mcmc),as.matrix(model.jags$mcmc)[,"r"]),2,quantile,probs=c(0.025,0.5,0.975))
