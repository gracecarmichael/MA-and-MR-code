library(MendelianRandomization)
library(tidyverse)
library(devtools)
#install_github("phenoscanner/phenoscanner")
library(phenoscanner)
library(MRInstruments) 
library(TwoSampleMR)
library(readr)
library(LDlinkR)
library(doParallel)
library(kableExtra)

### Load in the data ###
# Using GWAS Atlas data 

# These are the 77 SNPS we need from the obesity study
BMI_SNPS = read_csv("topSNPs.csv")
BMI_SNPS = as.data.frame(BMI_SNPS)
# These are all the SNPs associated with obesity
nd = read.table("allSNPs.uniq", header = T)
# extract full information about 77 SNPs from study
dat = nd[which(nd$SNP %in% BMI_SNPS$rsID), ]
# extract information on proxy SNP
rs9581854 = nd[which(nd$SNP == 'rs12016871'), ]
rs9581854$SNP = 'rs9581854'
# merge data and proxy
dat = rbind(dat,rs9581854)
dat = dat[order(dat[,1],decreasing = TRUE),]
BMI_SNPS = BMI_SNPS[order(BMI_SNPS[,3],decreasing = TRUE),]
# extract chromosome from dataset with 77 SNPs and add to dataset with full information
dat = cbind(dat,BMI_SNPS$CHR,rep('bmi',77))
# respecify column names
colnames(dat) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "N",'chr', 'phenotype')

# format the data for the TwoSampleMR package
bmi_exp_dat <- format_data(dat, phenotype_col = 'phenotype')


### Check for linkage disequilibrium on each chromosome ###

chr_vec = bmi_exp_dat$chr.exposure

# Register the cluster with 2 cores
cl <- makeCluster(6)
registerDoParallel(cl)

# The two methods for calculating the LD matrix
calc_ld = ld_matrix
calc_ld2 = LDmatrix

LD_values = foreach(i=1:19, .combine='c', .multicombine=TRUE) %dopar% {
  chr <- which(chr_vec == i)
  chr <- bmi_exp_dat$SNP[chr]
  
  # LDlinkR version
  result <- tryCatch(
    {
      calc_ld2(chr, 
               pop = "CEU", 
               r2d = "r2", 
               token = "b22ecb07e3d4", 
               file = FALSE, 
               genome_build = "grch37")
    },
    error=function(cond) {
      # Choose a return value in case of error
      return(NULL)
    },
    finally={
      message("Finished processing")
    }
  )    
  
  return(list(result))
}

LD_value2 = foreach(i=1:19, .combine='c', .multicombine=TRUE) %dopar% {
  chr <- which(chr_vec == i)
  chr <- bmi_exp_dat$SNP[chr]
  
  # LDlinkR version
  result <- tryCatch(
    {
      calc_ld(chr)
    },
    error=function(cond) {
      # Choose a return value in case of error
      return(NULL)
    },
    finally={
      message("Finished processing")
    }
  )    
  
  return(list(result))
}

stopCluster(cl)

# Just change i to view the matrices

#LD_values[[i]]

for (i in 1:19) {
  print(LD_value2[[i]])
}



### Get the outcome data ###
asthma_out_dat <- extract_outcome_data(
  snps = bmi_exp_dat$SNP,
  outcomes = 'ebi-a-GCST90014325', # Study on asthma
  proxies = TRUE, # Use a proxy if we cannot find the SNP from the exposure data in the outcome data
  rsq = 0.7 # This is the r2 value from the study
)


### Harmonise the datasets ###
bmi_asthma_dat = harmonise_data(bmi_exp_dat, asthma_out_dat, action = 2)
bmi_asthma_dat$exposure = "Body mass index"
bmi_asthma_dat$outcome = "Asthma"

#### FREQUENTIST MR ####

### Conducting the MR using the Mendelian Randomization package ###

# The values that the TwoSampleMR package actually uses for MR - remove palindromic SNP
bmi_asthma_dat_keep = bmi_asthma_dat[which(bmi_asthma_dat$mr_keep == T),]

# specify inputs for model
bx = bmi_asthma_dat_keep$beta.exposure
bxse = bmi_asthma_dat_keep$se.exposure
by = bmi_asthma_dat_keep$beta.outcome
byse = bmi_asthma_dat_keep$se.outcome

# create input data for model
MRInputObject = mr_input(bx = bx, 
                         bxse = bxse,
                         by = by,
                         byse = byse,
                         exposure = "Obesity: Body mass index",
                         outcome = "Asthma",
                         snps = bmi_asthma_dat_keep$SNP)

# produce outputs for all methods (IVE and MR-Egger included)
MRAllObject_all <- mr_allmethods(MRInputObject, method = "all")
MRAllObject_all



### Plots summarising results ###

# Scatter plot for each method
mr_plot(MRInputObject, orientate = F, interactive = F, line = "ivw")
mr_plot(MRInputObject, orientate = F, interactive = F, line = "egger")

# Forest plot A - SNPs 
mr_forest(MRInputObject)

# Forest plot B - Methods
mr_forest(MRInputObject,
          snp_estimates=FALSE,
          methods = c("ivw", "egger"))

# Funnel plot
mr_funnel(MRInputObject)

# Leave-one-out plot
mr_loo(MRInputObject)


### Conducting the MR using TwoSampleMR package ###

# Perform the MR
res <- mr(bmi_asthma_dat)
# Get results on the OR scale
res_or = generate_odds_ratios(res)
res_or = res_or[c(1,3), ]
res_or

## Additional analyses (model checks)

# Heterogeneity
mr_heterogeneity(bmi_asthma_dat)

# Horizontal pleiotropy
mr_pleiotropy_test(bmi_asthma_dat)


# Leave-one-out analysis 
res_loo <- mr_leaveoneout(bmi_asthma_dat)


### Plots summarising results ###
# Scatter plot
p1 <- mr_scatter_plot(res_or, bmi_asthma_dat)
p1[[1]]

# Forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]

#LOO plot
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

# Funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]




#### BAYESIAN MR ####

library('mrbayes')
library(gridExtra)

# Set the random seed for reproducible results
set.seed(2022)

# Convert an object of class MRInput from the 
# MendelianRandomization package to the mrbayes mr_format class
bayes_data = mrinput_mr_format(MRInputObject)


# Bayesian MR Egger model with a 
# choice of prior distributions fitted using RStan.
bayes_egger = mr_egger_stan(
  bayes_data,
  prior = 2,
  n.chains = 3,
  n.burn = 1000,
  n.iter = 5000,
  seed = 12345,
  rho = 0.5
)


# Bayesian inverse variance weighted model with a 
# choice of prior distributions fitted using RStan.
bayes_ivw = mr_ivw_stan(
  bayes_data,
  prior = 2,
  n.chains = 3,
  n.burn = 1000,
  n.iter = 5000,
  seed = 12345
)

# Save and load the results
#save(bayes_ivw, file = 'bayes_ivw.Rdata')
#save(bayes_egger, file = 'bayes_egger.Rdata')

load('bayes_ivw.Rdata')
load('bayes_egger.Rdata')

# Display the Bayesian analysis results 
bayes_ivw
bayes_egger

# Save the results so that they can be converted
# to the odds-ratio scale
sum_ivw = summary(bayes_ivw)
exp(sum_ivw$summary[1,c(4,8)])

sum_egger = summary(bayes_egger)
round(sum_egger$summary[1:3, c(4, 8)], 2)

# Get the posterior samples from the models
post.samps.ivw = as.data.frame(rstan::extract(bayes_ivw))
names(post.samps.ivw) = c("estimate")

post.samps.egger = as.data.frame(rstan::extract(bayes_egger))[,1:3]

# Convert the results to the odds-ratio scale
exp.samps.ivw = exp(post.samps.ivw)
exp.samps.egger = exp(post.samps.egger)


head(rstan::As.mcmc.list(bayes_ivw))
head(as.data.frame(rstan::extract(bayes_ivw)))




### Plots summarising results ###

# Traceplots to confirm that the chains have converged
rstan::traceplot(bayes_ivw)
rstan::traceplot(bayes_egger)


# Plot of the IVW Estimate
ivw = ggplot(aes(x = estimate), data = post.samps.ivw) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samps.ivw$estimate)) +
  labs(x = "Overall IVW Estimate",
       y = element_blank()) +
  theme_minimal()

# Plot of the IVW Estimate (OR)
ivw.or = ggplot(aes(x = estimate), data = exp.samps.ivw) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(exp.samps.ivw$estimate)) +
  geom_vline(xintercept = 1, color = "gray", size = 1, linetype = 2) +
  labs(x = "Overall IVW Estimate (Odds-Ratio)",
       y = element_blank()) +
  theme_minimal()

# Plot of the MR Egger Intercept
egger.int = ggplot(aes(x = intercept), data = post.samps.egger) +
  geom_vline(xintercept = 0, color = "gray", size = 1, linetype = 2) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samps.egger$intercept)) +
  labs(x = "Egger Intercept",
       y = element_blank()) +
  theme_minimal()

# Plot of the MR Egger Estimate
egger.est = ggplot(aes(x = estimate), data = post.samps.egger) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samps.egger$estimate)) +
  labs(x = "Overall Egger Estimate",
       y = element_blank()) +
  theme_minimal()


# Plot of the MR Egger Estimate (OR)
egger.est.or = ggplot(aes(x = estimate), data = exp.samps.egger) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(exp.samps.egger$estimate)) +
  geom_vline(xintercept = 1, color = "gray", size = 1, linetype = 2) +
  labs(x = "Overall Egger Estimate (Odds-Ratio)",
       y = element_blank()) +
  theme_minimal()

# Plot of the MR Egger Sigma
egger.sig = ggplot(aes(x = sigma), data = post.samps.egger) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samps.egger$sigma)) +
  labs(x = "Egger Sigma",
       y = element_blank()) +
  theme_minimal()



# Display the plots
ivw.or
egger.int
egger.est.or
egger.sig



