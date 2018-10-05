
modelName <- "model {
#run through all the individual observations
for(i in 1:N) {
y[i] ~ dnorm(mu[i],tau)
mu[i] <- 
} # end the calculation of the likelihood

###########################
### priors by population ##
###########################

#for(j in 1:Ngroup) {
#
#} # end the priors for each group



############################################
# hyper priors on parameters of model ######
############################################

m ~ 
b ~ 

#########################################################
# below are the half-t priors on variance parameters ####
#########################################################

sd ~ 

#############################################
# conversion to JAGS precision parameters ###
#############################################

tau <- 

###############################################
# re-parameterize  parameters of the model ####
###############################################
}"
