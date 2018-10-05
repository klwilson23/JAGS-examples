rm(list=ls(all=TRUE))

# Set the working directory, example below:
setwd('~/Google Drive/Biphasic Growth Methods Manuscript/Drafts/Coded Examples/JAG outcome')


library(runjags) #load, install, or require the 'runjags' library
library(stats4) #load the 'stats4' library
library(coda)


#############################
## simulation of regression #
#############################
N <- 50
Ngroups <- 10
m.true <- 2.5
m.cv <- 0.25
m_i.true <- rnorm(Ngroups,mean=m.true,sd=m.true*m.cv)
b.true <- 6
b.cv <- 1.5
b_i.true <- rnorm(Ngroups,mean=b.true,sd=b.true*b.cv)
x <- seq(from=1,to=100,length=N)
sd.true <- 10

data <- NULL
for(i in 1:Ngroups)
{
  y.true <- m_i.true[i]*x + b_i.true[i]
  y.ran <- rnorm(N,mean=y.true,sd=sd.true)
  sub.data <- cbind(x,y.true,y.ran,rep(i,N))
  data <- rbind(data,sub.data)
}

data <- as.data.frame(data)
colnames(data) <- c("x","y.true","y.ran","group")
plot(data$x,data$y.true,type='c',ylim=c(range(data$y.ran)),col="black",lwd=2,lty=2)
points(data$x,data$y.ran,col=data$group)

############################
## JAGS model script #######
############################

sim3 <- "model {
#run through all the individual observations
for(i in 1:N) {
y[i] ~ dnorm(mu[i],tau)
mu[i] <- m_i[group[i]]*x[i]+b_i[group[i]]
} # end the calculation of the likelihood

for(j in 1:Ngroups) {
m_i[j] ~ dnorm(m,tau.m)
b_i[j] ~ dnorm(b,tau.b)
}

m ~ dunif(0,5)
b ~ dnorm(0,1e-4)
m.cv ~ dunif(0,3)
b.cv ~ dunif(0,3)
sd ~ dunif(0,1000)

tau <- 1/sd^2
tau.m <- 1/(m*m.cv)^2
tau.b <- 1/(b*b.cv)^2
}"

##########################################
# initialize the parameters of the model #
##########################################
fit.lm <- lm(y.ran~x,data=data)
m <- as.numeric(coef(fit.lm)[2])
b <- as.numeric(coef(fit.lm)[1])
sd <- mean(abs(data$y.ran[1:10]))*sd(data$y.ran)/mean(data$y.ran)

inits1 <- list(
  m=m,
  b=b,
  m_i=rep(m,Ngroups),
  b_i=rep(b,Ngroups),
  m.cv=0.20,
  b.cv=0.30,
  sd=sd,
  .RNG.name="base::Wichmann-Hill", 
  .RNG.seed=100
)

inits3 <- inits2 <- inits1

rngList <- function(x,y){ #this function creates randomly 'jittered' starting values for each MCMC chain
  lis <- lapply(x, lapply, length) #get the lower order dimensions of the list x
  names(lis) <- lapply(x, length) #get the names of those dimensions of the list x
  l_el <- length(names(lis)) #get the maximum number elements of the highest order dimensions of the list x
  for(i in 1:(l_el-2)) #loop through those dimensions which need to be 'jittered'
  {
    x[[i]] <- x[[i]]*(1+runif(1,-0.15,0.15)) #jitter values of the list by +/- 15%
  }
  x[[l_el]] <- round(runif(1,1,1e7),0) #have a random RNG seed for the MCMC chain
  return(x)
}

inits2 <- rngList(inits2,inits1) # jitter chain 2, based on values of chain 1
inits3 <- rngList(inits3,inits1) # jitter chain 3, based on values of chain 1

inits <- list(inits1,inits2,inits3) # compile all initial values into one list


mon_names <- names(inits1)[-c(length(inits1),length(inits1)-1)]


##############################################
# compile the data for the model into a list #
##############################################

data.compiled <- list(N=length(data$x),
             Ngroups=length(unique(data$group)),
             x=data$x,
             y=data$y.ran,
             group=data$group)

##########################################
# run the JAGS model #####################
##########################################

Nsamp <- 200 # how many posterior samples does each chain need to get, after thinning and burin-in and adaptation?
thin_rt <- 50 # place some sort of thinning rate?
burnins <- 0.50*round(Nsamp*thin_rt,0) # how long is the burnin, this bases it on the number of total posterior draws?
adaptin <- round(1*burnins,0)

start <- proc.time();
fit <- run.jags(model=sim3, monitor=mon_names, 
                data=data.compiled, n.chains=3, method="rjags", inits=inits,
                plots=F,silent.jag=F, modules=c("bugs","glm","dic"),
                sample=Nsamp,adapt=adaptin,burnin=burnins,thin=thin_rt,summarise=TRUE)
end <- (proc.time() - start)

print((end[3]/60)/60)
print((start[3]/(3*(Nsamp*thin_rt+burnins))))
summary(fit)
################################################
# store the outcome of the model's MCMC chains #
################################################
coeff <- as.mcmc.list(fit, vars=mon_names)
mcmc_res <- rbind(coeff[[1]],coeff[[2]],coeff[[3]])

####################################################
# graph the marginal posteriors for each parameter #
####################################################
summary(mcmc_res)
summary(fit)
colnames(mcmc_res)
hist(mcmc_res[,"m"],col=rgb(0.1,0.1,0.1,0.2),main='',xlab="Marginal posterior distribution of m")
hist(coeff[[1]][,"m"],add=TRUE,col=rgb(1,0.1,0.1,0.2))
hist(coeff[[2]][,"m"],add=TRUE,col=rgb(0.3,0.1,0.1,0.2))
hist(coeff[[3]][,"m"],add=TRUE,col=rgb(0.6,0.1,0.1,0.2))
abline(v=quantile(mcmc_res[,"m"],probs=c(0.025,0.5,0.975)),lty=c(2,1,2),col="black")
abline(v=m.true,lty=1,lwd=3,col="purple")

hist(mcmc_res[,"b"],main='',xlab="Marginal posterior distribution of b")
abline(v=quantile(mcmc_res[,"b"],probs=c(0.025,0.5,0.975)),lty=c(2,1,2),col="black")
abline(v=b.true,lty=1,lwd=3,col="purple")

hist(mcmc_res[,"sd"],main='',xlab="Marginal posterior distribution of sigma")
abline(v=quantile(mcmc_res[,"sd"],probs=c(0.025,0.5,0.975)),lty=c(2,1,2),col="black")
abline(v=sd.true,lty=1,lwd=3,col="purple")

######################################################
# run diagnostics on the MCMC and posterior outcomes #
######################################################
# http://svitsrv25.epfl.ch/R-doc/library/coda/html/00Index.html

layout(matrix(1,nrow=1,ncol=1))
plot(coeff)
geweke.diag(coeff)
geweke.plot(coeff)

traceplot(coeff)
gelman.diag(coeff)
gelman.plot(coeff,ask=FALSE)
heidel.diag(coeff)

#layout(matrix(1:9,nrow=3,ncol=3,byrow=T))
cumuplot(coeff,auto.layout=T,ask=FALSE)
effectiveSize(coeff)

autocorr.plot(coeff,auto.layout=T,ask=FALSE)
autocorr.diag(coeff)
crosscorr(coeff)
layout(matrix(1,nrow=1,ncol=1))
crosscorr.plot(coeff)

# codamenu()

y.hat <- sapply(x,FUN=function(x){as.vector(mcmc_res[,"m"])*x + as.vector(mcmc_res[,"b"])})
pred.inter <- apply(y.hat,2,FUN=quantile,probs=c(0.025,0.5,0.975))
xs <- c(x,rev(x))

plot(x,pred.inter[2,],ylim=range(pred.inter),col=0,xlab="Independent Variable",ylab="Dependent Variable")
polygon(xs,c(pred.inter[1,],rev(pred.inter[3,])),col="grey50")
lines(x,pred.inter[2,],col="black",lwd=2)
lines(x,y.true,col="blue",lwd=2,lty=2)
abline(fit.lm,lty=3,col="orange",lwd=3)
points(x,y.ran,pch=21,bg="white")
