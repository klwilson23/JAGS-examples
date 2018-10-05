rm(list=ls(all=TRUE))

# Set the working directory, example below:
setwd('C:/Users/Kyle/Google Drive/Biphasic Growth Methods Manuscript/Drafts/Coded Examples/JAG outcome')


library(runjags) #load, install, or require the 'runjags' library
library(stats4) #load the 'stats4' library
library(coda)


#############################
## simulation of regression #
#############################
N <- 100
m.true <- 2.5
b.true <- 6
x <- seq(from=1,to=100,length=N)
y.true <- m.true*x + b.true
sd.true <- 30
y.ran <- rnorm(N,mean=y.true,sd=sd.true)
plot(x,y.true,type='l',ylim=c(range(y.ran)))
points(x,y.ran)

############################
## JAGS model script #######
############################

sim1 <- "model {
#run through all the individual observations
for(i in 1:N) {
y[i] ~ dnorm(mu[i],tau)
mu[i] <- m*x[i]+b
} # end the calculation of the likelihood

#m ~ dnorm(0,1e-4)
m ~ dnorm(10,2)T(0,)
b ~ dnorm(0,1e-4)

sd ~ dunif(0,1000)

tau <- 1/sd^2

}"

##########################################
# initialize the parameters of the model #
##########################################
fit.lm <- lm(y.ran~x)
m <- diff(range(y.ran))/diff(range(x))
b <- min(y.ran)
sd <- abs(y.ran[1])*sd(y.ran)/mean(y.ran)
inits1 <- list(
  m=m,
  b=b,
  sd=sd,
  .RNG.name="base::Wichmann-Hill", 
  .RNG.seed=735
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

data <- list(N=length(x),
             x=x,
             y=y.ran)

##########################################
# run the JAGS model #####################
##########################################

Nsamp <- 300 # how many posterior samples does each chain need to get, after thinning and burin-in and adaptation?
thin_rt <- 50 # place some sort of thinning rate?
burnins <- 0.20*round(Nsamp*thin_rt,0) # how long is the burnin, this bases it on the number of total posterior draws?
adaptin <- round(1*burnins,0)

start <- proc.time();
fit <- run.jags(model=sim1, monitor=mon_names, 
                    data=data, n.chains=3, method="rjags", inits=inits,
                    plots=F,silent.jag=F, #modules=c("bugs","glm","dic")
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
gelman.plot(coeff)
heidel.diag(coeff)

layout(matrix(1:9,nrow=3,ncol=3,byrow=T))
cumuplot(coeff,auto.layout=FALSE)
effectiveSize(coeff)

autocorr.plot(coeff,auto.layout=FALSE)
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
