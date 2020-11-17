###########################
##  Packages & Settings
###########################

rm(list=ls())  # remove all variables 

library(rstan)
library(ggplot2)
library(gridExtra)
library(formattable)
library(bayesplot)
library(cowplot)
library(rstanarm)

# source HDIofMCMC.R to calculate HDI
setwd("/home/heesunpark/project/Meta-d_STAN/codes/")
source("HDIofMCMC.R") 


PATH_DATA = "/home/heesunpark/project/Meta-d_STAN/example_data/"
PATH_FIGURE = "/home/heesunpark/project/Meta-d_STAN/figure/"

DATA_NAME = "contrast_perception_counts.txt" #"orientation_perception_counts.txt" / "orientation_memory_counts.txt"
FIGURE_VERSION = "Indiv_con_percept_"

######################################################
##  Load the Data and Process for use in STAN model
######################################################

data = read.table(paste0(PATH_DATA, DATA_NAME), header=T, sep="") 

#data setting
data <- as.matrix(data) 

allSubjs =  c(1:dim(data)[1]) # all subject IDs
N = length(allSubjs)      # number of subjects
Trials = sum(data[1,]) # number of trials per subject (=300)
nratings = length(data[1,])/4

counts <- array(0, c(N, length(data[1,])))
H <- array(0, c(N))
FA <- array(0, c(N))
CR <- array(0, c(N))
M <- array(0, c(N))

for (i in 1:N){
  counts[i,] <- data[i,]
}

for (i in 1:N){
  H[i] <- sum(counts[i,((nratings*3)+1):(nratings*4)])
  FA[i] <- sum(counts[i,(nratings+1):(nratings*2)])
  CR[i] <- sum(counts[i, 1:(nratings)])
  M[i] <- sum(counts[i, ((nratings*2)+1):(nratings*3)])
}


# Adjust to ensure non-zero counts for type 1 d' point estimate
d1 <- data.frame()
c1 <- data.frame()

for (i in 1:N){
  adj_f <- 1/((nratings)*2)
  nR_S1_adj = counts[i,1:(nratings*2)] + adj_f
  nR_S2_adj = counts[i,((nratings*2)+1):(nratings*4)] + adj_f
  
  ratingHR <- matrix()
  ratingFAR <- matrix()
  for (c in 2:(nratings*2)) {
    ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
    ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
    
  }
  t1_index <- nratings
  d1_tmp <- qnorm(ratingHR[(t1_index)]) - qnorm(ratingFAR[(t1_index)])
  d1 <- rbind(d1, d1_tmp)
  c1_tmp <- -0.5 * (qnorm(ratingHR[(t1_index)]) + qnorm(ratingFAR[(t1_index)]))
  c1 <- rbind(c1, c1_tmp)
}

d1 <- as.vector(as.matrix(d1))
c1 <- as.vector(as.matrix(c1))


dataList <- list(
  N = N,
  H = H,   
  FA = FA,
  CR = CR,
  M = M,
  nratings = nratings,
  counts = counts,
  d1 = d1,
  c1 = c1
)


#####################################
##  Run STAN model
#####################################

output = stan("metad_indiv.stan", data = dataList, pars = c("meta_d", "cS1", "cS2", "log_lik", "y_pred", "Mratio"),
              iter = 8000, warmup=4000, chains=4, cores=20)

# save(output, file="/Users/heesunpark/project/Meta-d_STAN/Indiv.RData")


###############################################
##  Check trace plots
###############################################

#traceplot
trace_meta <- traceplot(output, pars="meta_d")
ggsave(paste0(PATH_FIGURE,FIGURE_VERSION,'_trace_meta.jpg'), trace_meta,width=40, height=24, units="cm",dpi=300)

trace_cS1 <- traceplot(output, pars="cS1")
ggsave(paste0(PATH_FIGURE,FIGURE_VERSION,'_trace_cS1.jpg'), trace_cS1, width=40, height=24, units="cm", dpi=300)

trace_cS2 <- traceplot(output, pars="cS2")
ggsave(paste0(PATH_FIGURE,FIGURE_VERSION,'_trace_cS2.jpg'), trace_cS2, width=40, height=24, units="cm", dpi=300)


##################################################
##  Posterior predictive check
##################################################


#extract data
extract_output <- rstan::extract(output)

for (j in 1:(nratings*4)){
  
  PP_dens_tmp <- ppc_dens_overlay(counts[,j], extract_output$counts_pred[1:50,,j])+ggtitle(paste0("counts[",j,"]")) #for all subj, count[,i]
  assign(paste0("PP_dens", j), PP_dens_tmp)
  PP_hist_tmp <- ppc_hist(counts[,j],extract_output$counts_pred[1:50,,j])+ggtitle(paste0("counts[",j,"]")) #for all subj, count[,i]
  assign(paste0("PP_hist",j), PP_hist_tmp)
}

#save distribution
a1 <- grid.arrange(get(paste("PP_dens",1,sep="")), get(paste("PP_dens",2,sep="")), get(paste("PP_dens",3,sep="")), get(paste("PP_dens",4,sep="")), ncol=4)
a2 <- grid.arrange(get(paste("PP_dens",5,sep="")), get(paste("PP_dens",6,sep="")), get(paste("PP_dens",7,sep="")), get(paste("PP_dens",8,sep="")), ncol=4)
a3 <- grid.arrange(get(paste("PP_dens",9,sep="")), get(paste("PP_dens",10,sep="")), get(paste("PP_dens",11,sep="")), get(paste("PP_dens",12,sep="")), ncol=4)
a4 <- grid.arrange(get(paste("PP_dens",13,sep="")), get(paste("PP_dens",14,sep="")), get(paste("PP_dens",15,sep="")), get(paste("PP_dens",16,sep="")), ncol=4)

ppc_distribution <- grid.arrange(a1,a2,a3,a4, nrow=4)

ggsave(paste0(PATH_FIGURE, FIGURE_VERSION,'_PPC_distribution.jpg'), ppc_distribution, width=40, height=24, units="cm", dpi=300)

#save hist
for (j in 1:(nratings*4)){
  ggsave(paste0(PATH_FIGURE, FIGURE_VERSION,'_PPC_hist_counts[',j,'].jpg'), get(paste("PP_hist",j,sep="")), width=40, height=24, units="cm", dpi=300)
}


##################################################
## Plot distributions
##################################################

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")
pd_meta_d <- mcmc_areas(output, regex_pars = 'meta_d', prob = 0.95) + plot_title
ggsave(paste0(PATH_FIGURE, FIGURE_VERSION,'_posterior_distribution_meta_d.jpg'), pd_meta_d, width=40, height=24, units="cm", dpi=300)

pd_cS1 <- mcmc_intervals(output, pars = c(paste0('cS1[', c(1:K), ',1]'), paste0('cS1[', c(1:K), ',2]'), paste0('cS1[', c(1:K), ',3]')), prob = 0.95) + plot_title
ggsave(paste0(PATH_FIGURE, FIGURE_VERSION,'_posterior_distribution_cS1.jpg'), pd_cS1, width=40, height=24, units="cm", dpi=300)

pd_cS2 <- mcmc_intervals(output, pars = c(paste0('cS2[', c(1:K), ',1]'), paste0('cS2[', c(1:K), ',2]'), paste0('cS2[', c(1:K), ',3]')), prob = 0.95) + plot_title
ggsave(paste0(PATH_FIGURE, FIGURE_VERSION,'_posterior_distribution_cS2.jpg'), pd_cS2, width=40, height=24, units="cm", dpi=300)

pd_Mratio <- mcmc_areas(output, regex_pars = 'Mratio', prob = 0.95) + plot_title
ggsave(paste0(PATH_FIGURE, FIGURE_VERSION,'_posterior_distribution_Mratio.jpg'), pd_Mratio,width=40, height=24, units="cm", dpi=300)



##################################################
## Check rhats
##################################################

rhats <- rhat(output) # bayesplot package

rhats_meta_d <- NULL
for (i in 1:K){rhats_meta_d <- rbind(rhats_meta_d, rhats[paste0("meta_d[",i,"]")])}
metarhat <- mcmc_rhat(as.vector(rhats_meta_d)) + ggtitle("rhat of meta-d")

rhats_cS1_cS2 <- NULL
for (i in 1:K){
  rhats_cS1_cS2 <- rbind(rhats_cS1_cS2, rhats[paste0("cS1[",i,",1]")], rhats[paste0("cS1[",i,",2]")], rhats[paste0("cS1[",i,",3]")], rhats[paste0("cS2[",i,",1]")], rhats[paste0("cS2[",i,",2]")], rhats[paste0("cS2[",i,",3]")])
}
cS1S2rhat <- mcmc_rhat(as.vector(rhats_cS1_cS2)) + ggtitle("rhat of cS1 and cS2")
rhatplots <- grid.arrange(metarhat, cS1S2rhat, ncol=3)
ggsave(paste0(PATH_FIGURE,FIGURE_VERSION,"_rhat.jpg"), rhatplots, width=33, height=10, units="cm",dpi=300)