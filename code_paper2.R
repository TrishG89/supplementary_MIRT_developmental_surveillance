#libraries and functions
source("DBDA2E-utilities.R") # Must be in R's current working directory - functions provided in "Doing Bayesian Data Analysis" by John Kruschke.
fileNameRoot = "Jags-ItemResponseTheory-"

logistic = function( x , g=1 , t=0 ) {
  return( 1/(1+exp(-(g*(x-t)))) )
}
library(tidyverse)
library(coda)
library(e1071)
library(lattice)
library(factoextra)
library(RColorBrewer)
library(cluster)
library(NbClust)

#load library for Rstan
#loading package
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

#Data set up --------------------------------------------------------------------------------------------------------------------------------------------------
#load data
milestones<-read_csv("data_for_IRM.csv")
milestones$areacode<- factor(milestones$areacode, levels = c("A_1.1", "A_1.2", "A_1.3",
                                                             "A_2.1", "A_2.2", "A_2.3",
                                                             "A_3.1", "A_3.2", "A_3.3",
                                                             "A_4.1", "A_4.2", "A_4.3",
                                                             "A_5.1", "A_5.2", "A_5.3",
                                                             "A_6.1", "A_6.2", "A_6.3",
                                                             "A_7.1", "A_7.2", "A_7.3",  
                                                             "A_8.1", "A_8.2", "A_8.3",
                                                             "A_9.1", "A_9.2", "A_9.3",
                                                             "A_10.1","A_10.2","A_10.3",
                                                             "A_11.1","A_11.2","A_11.3",
                                                             "A_12.1","A_12.2","A_12.3",
                                                             "A_13.1","A_13.2","A_14.1",
                                                             "A_14.2","A_15.1","A_15.2",
                                                             "A_16.1","A_16.2","A_17.1",
                                                             "A_17.2","A_18.1","A_18.2",
                                                             "A_19.1","A_20.1","A_21.1",
                                                             "A_22.1","A_23.1","A_24.1",
                                                             "A_25.1","A_28.1","A_31.1","A_34.1",
                                                             "H_1.1", "H_1.2", "H_1.3",
                                                             "H_2.1", "H_2.2", "H_2.3",
                                                             "H_3.1", "H_3.2", "H_3.3",
                                                             "H_4.1", "H_4.2", "H_4.3",
                                                             "H_5.1", "H_5.2", "H_5.3",
                                                             "H_6.1", "H_6.2", "H_6.3",
                                                             "H_7.1", "H_7.2", "H_7.3",  
                                                             "H_8.1", "H_8.2", "H_8.3",
                                                             "H_9.1", "H_9.2", "H_9.3",
                                                             "H_10.1","H_10.2","H_10.3",
                                                             "H_11.1","H_11.2","H_11.3",
                                                             "H_12.1","H_12.2","H_12.3",
                                                             "H_13.1","H_13.2","H_14.1",
                                                             "H_14.2","H_15.1","H_15.2",
                                                             "H_16.1","H_16.2","H_17.1",
                                                             "H_17.2","H_18.1","H_18.2",
                                                             "H_19.1","H_20.1","H_21.1",
                                                             "H_22.1","H_23.1","H_24.1",
                                                             "H_25.1","H_28.1","H_31.1","H_34.1",
                                                             "M_1.1", "M_1.2", "M_1.3",
                                                             "M_2.1", "M_2.2", "M_2.3",
                                                             "M_3.1", "M_3.2", "M_3.3",
                                                             "M_4.1", "M_4.2", "M_4.3",
                                                             "M_5.1", "M_5.2", "M_5.3",
                                                             "M_6.1", "M_6.2", "M_6.3",
                                                             "M_7.1", "M_7.2", "M_7.3",  
                                                             "M_8.1", "M_8.2", "M_8.3",
                                                             "M_9.1", "M_9.2", "M_9.3",
                                                             "M_10.1","M_10.2","M_10.3",
                                                             "M_11.1","M_11.2","M_11.3",
                                                             "M_12.1","M_12.2","M_12.3",
                                                             "M_13.1","M_13.2","M_14.1",
                                                             "M_14.2","M_15.1","M_15.2",
                                                             "M_16.1","M_16.2","M_17.1",
                                                             "M_17.2","M_18.1","M_18.2",
                                                             "M_19.1","M_20.1","M_21.1",
                                                             "M_22.1","M_23.1","M_24.1",
                                                             "M_25.1","M_28.1","M_31.1","M_34.1",
                                                             "S_1.1", "S_1.2", "S_1.3",
                                                             "S_2.1", "S_2.2", "S_2.3",
                                                             "S_3.1", "S_3.2", "S_3.3",
                                                             "S_4.1", "S_4.2", "S_4.3",
                                                             "S_5.1", "S_5.2", "S_5.3",
                                                             "S_6.1", "S_6.2", "S_6.3",
                                                             "S_7.1", "S_7.2", "S_7.3",  
                                                             "S_8.1", "S_8.2", "S_8.3",
                                                             "S_9.1", "S_9.2", "S_9.3",
                                                             "S_10.1","S_10.2","S_10.3",
                                                             "S_11.1","S_11.2","S_11.3",
                                                             "S_12.1","S_12.2","S_12.3",
                                                             "S_13.1","S_13.2","S_14.1",
                                                             "S_14.2","S_15.1","S_15.2",
                                                             "S_16.1","S_16.2","S_17.1",
                                                             "S_17.2","S_18.1","S_18.2",
                                                             "S_19.1","S_20.1","S_21.1",
                                                             "S_22.1","S_23.1","S_24.1",
                                                             "S_25.1","S_28.1","S_31.1","S_34.1",
                                                             "T_1.1", "T_1.2", "T_1.3",
                                                             "T_2.1", "T_2.2", "T_2.3",
                                                             "T_3.1", "T_3.2", "T_3.3",
                                                             "T_4.1", "T_4.2", "T_4.3",
                                                             "T_5.1", "T_5.2", "T_5.3",
                                                             "T_6.1", "T_6.2", "T_6.3",
                                                             "T_7.1", "T_7.2", "T_7.3",  
                                                             "T_8.1", "T_8.2", "T_8.3",
                                                             "T_9.1", "T_9.2", "T_9.3",
                                                             "T_10.1","T_10.2","T_10.3",
                                                             "T_11.1","T_11.2","T_11.3",
                                                             "T_12.1","T_12.2","T_12.3",
                                                             "T_13.1","T_13.2","T_14.1",
                                                             "T_14.2","T_15.1","T_15.2",
                                                             "T_16.1","T_16.2","T_17.1",
                                                             "T_17.2","T_18.1","T_18.2",
                                                             "T_19.1","T_20.1","T_21.1",
                                                             "T_22.1","T_23.1","T_24.1",
                                                             "T_25.1","T_28.1","T_31.1","T_34.1",
                                                             "V_1.1", "V_1.2", "V_1.3",
                                                             "V_2.1", "V_2.2", "V_2.3",
                                                             "V_3.1", "V_3.2", "V_3.3",
                                                             "V_4.1", "V_4.2", "V_4.3",
                                                             "V_5.1", "V_5.2", "V_5.3",
                                                             "V_6.1", "V_6.2", "V_6.3",
                                                             "V_7.1", "V_7.2", "V_7.3",  
                                                             "V_8.1", "V_8.2", "V_8.3",
                                                             "V_9.1", "V_9.2", "V_9.3",
                                                             "V_10.1","V_10.2","V_10.3",
                                                             "V_11.1","V_11.2","V_11.3",
                                                             "V_12.1","V_12.2","V_12.3",
                                                             "V_13.1","V_13.2","V_14.1",
                                                             "V_14.2","V_15.1","V_15.2",
                                                             "V_16.1","V_16.2","V_17.1",
                                                             "V_17.2","V_18.1","V_18.2",
                                                             "V_19.1","V_20.1","V_21.1",
                                                             "V_22.1","V_23.1","V_24.1",
                                                             "V_25.1","V_28.1","V_31.1","V_34.1"))

milestones<-arrange(milestones, child_id, areacode)
#calculate the proportion of achieved milestones per child and proportion of subjects who achieved each milestone
subjPropCorr = aggregate(status ~ child_id , data=milestones , FUN="mean") 
itemPropCorr = aggregate(status ~ areacode , data=milestones , FUN="mean" ) 
#Arrange data into list for STAN ---------------------------------------------------------------------------------------------------------------------------------------------
milestones<-as.data.frame(milestones)
milestones_wide<-spread(milestones, areacode, status)
#store child-id's
child_id<-milestones_wide$child_id
milestones_wide<-select(milestones_wide, -child_id)
#response vector
y = milestones$status
y<-as.matrix(y)
#item ID for each observation
itemID = as.numeric(factor(milestones[,"areacode"]))
#subject ID for each observation
subjID = as.numeric(factor(milestones[,"child_id"]))
#total number of milestones
Nitem = length(unique(itemID)) 
#total number of observations
Ntotal = nrow(milestones)
#total number of children
Nsubj = length(unique(subjID))


##map each milestone to the functional domain that it measures

itemmap<-itemID

itemmap<-replace(itemmap, itemmap %in% c(1:58), 1)
itemmap<-replace(itemmap, itemmap %in% c(59:116), 2)
itemmap<-replace(itemmap, itemmap %in% c(117:174), 3)
itemmap<-replace(itemmap, itemmap %in% c(175:232), 4)
itemmap<-replace(itemmap, itemmap %in% c(233:290), 5)
itemmap<-replace(itemmap, itemmap %in% c(291:348), 6)
#number of functional domains
K=6

#change to index names
I = Nsubj
J=Nitem
N = Ntotal
ii = subjID
jj = itemID
K=K #number of latent groups
kk = itemmap #maps item to latent variable
#hyperprior for Inverse Wishart for covariance
W=diag(6)
#store the data in a list
stan_input<-list(I=I, J=J, N=N, ii=ii, jj=jj, y=y, K=K, kk=kk, W=W)

#save data list
saveRDS(stan_input, "stan_input.rds")

#fit stan model ------------------------------------------------------------------------------------------------------------------------------------------------------------
twopl.fit<-stan(file = "2pl_multivariate_hierarchical_model.stan", data =stan_input, 
                pars = c("alpha", "beta", "mu_alpha", "sigma_alpha", "mu_theta", "Sigma_theta", "theta", "cor_sigma"),
                iter=50000, chains=4) 
#note this model can take many hours to run - the model for the paper took approximately 6 hours.

#check results
print(twopl.fit)
summary(twopl.fit)

#save MCMC draws
codaSamples = As.mcmc.list(twopl.fit)
save(codaSamples, file = "2pl_multivariate_hierarchical_reparam_cormat.Rdata")

#extract and store the posterior mean ability and standardise the posterior means ----------------------------------------------------------------------------------------------
parameterNames = varnames(codaSamples)
mcmcMat = as.matrix(codaSamples)
dim(mcmcMat)
#calculate the posterior mean for all parameters
mcmcMean = apply( mcmcMat , 2 , mean )
#extract the posterior means for the ability parameters
subjInfo = cbind( 
  round( mcmcMean[ grep( "theta" , names(mcmcMean) ) ] , 1 ) )
#split into the functional domains - check parameter names to identify the position on where to subset
subjInfo<-as.data.frame(subjInfo)
subset_auditory_subject<-subjInfo[43:157, ]
subset_hands_subject<-subjInfo[158:272, ]
subset_movement_subject<-subjInfo[273:387, ]
subset_speech_subject<-subjInfo[388:502, ]
subset_tactile_subject<-subjInfo[503:617, ]
subset_vision_subject<-subjInfo[618:732, ]

all_posterior_means<-cbind(subset_auditory_subject, subset_hands_subject, subset_movement_subject,
                           subset_speech_subject, subset_tactile_subject, subset_vision_subject)

all_posterior_means<-as.data.frame(all_posterior_means)
#cluster the posterior mean ability parameters using hierarchical clustering
#scale the data
scale_z<-scale(all_posterior_means)
#create dissimilarity matrices using Euclidean distance (root sum-of-squares of difference)
diss_euc_z<-daisy(scale_z, metric = "euclidean")
#hierarchical clustering using Ward linkage
hc_ward_euc_z <- hclust(diss_euc_z, method = "ward.D2" )

#plot dendrogram
plot(hc_ward_euc_z, cex = 0.6)
rect.hclust(hc_ward_euc_z, k = 3, border = 2:4)

#evaluate clustering
#30 different clustering validation tools
res.nbclust <- NbClust(data = scale_z, diss = NULL, distance = "euclidean", 
                       min.nc = 2, max.nc = 20, method = "ward.D2", index = "all", alphaBeale = 0.1)
