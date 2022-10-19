setwd("/Users/jared/Library/Group Containers/UBF8T346G9.OneDriveStandaloneSuite/OneDrive.noindex/OneDrive/Documents/University stuff/Honours Year/Thesis/Submissions/R code")

library(ggplot2)
library(metafor)
library("rstan") # observe startup messages
library("StanHeaders")
library(brms)
library(clubSandwich)
library(tidybayes)
library(tidyverse)
library(extraDistr)
library(tidybayes)
library(dplyr)
library(ggridges)
library(glue)
library(stringr)
library(forcats)
library(kableExtra)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Read in the data
dat = read.table('MA_data_v2.txt', sep = ',', header = T)
# Convert standard errors to variance
dat$var = (dat$StdErr)^2


### FREQUENTIST APPROACH ###

# Data without the Barros study.
dat2 = dat[-c(7,8),]

## Fit models

# Three layer model
mod1 = rma.mv(yi = EffectSize,
              V = var,
              slab = Paper,
              data = dat,
              random = ~1 | Paper/id,
              method = 'REML',
              test = 't')
summary(mod1)

# Two layer model
mod1_red = rma.mv(yi = EffectSize,
                  V = var,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't',
                  sigma2 = c(0,NA))

summary(mod1_red)
# Check if 3 layer model is significantly better than 2 layer
anova(mod1,mod1_red) 
# two layer model is favored over 3 but we have reason to 
# still stick with 3 level due to multiple estimates from each study. We know the effect 
# sizes from studies are not independent so nested model makes more sense.

# Set the correlation for the covariance matrix
rho = 0.5
# Impute the covariance matrix
V <- with(dat, 
          impute_covariance_matrix(vi = var,
                                   cluster = Paper,
                                   r = rho))
# Re-fit the three-level model using the covariance matrix from above
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't',
                  dfs="contain")
summary(che.mod)

# Get more accurate confidence intervals and p-values 
CI <- conf_int(che.mod, 
               vcov = "CR2")
p.val <- coef_test(che.mod, 
                   vcov = "CR2")

# Without the Barros Paper
V2 <- with(dat2, 
          impute_covariance_matrix(vi = var,
                                   cluster = Paper,
                                   r = rho))
che.mod2 <- rma.mv(yi = EffectSize,
                  V = V2,
                  slab = Paper,
                  data = dat2,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't',
                  dfs="contain")
summary(che.mod2)

CI2 <- conf_int(che.mod2, 
               vcov = "CR2")
p.val2 <- coef_test(che.mod2, 
                   vcov = "CR2")


## Produce summary plots of results

# Forest plot on log(OR) scale
forest.rma(che.mod)
# Forest plot on OR scale
forest.rma(che.mod, transf = exp, xlab = "Observed Outcome (Odds-Ratio)", showweights = TRUE, header = TRUE, addpred = F)
abline(v = 1, col = "red", lty = 2)

# Radial plot with Barros
radial.rma(che.mod, zlab = "Z-score", xlab = "Inverse Standard Error", center = T, transf = exp)
# Radial plot without Barros
radial.rma(che.mod2, zlab = "Z-score", xlab = "Inverse Standard Error", center = T, transf = exp)

# Funnel plot - Log-Odds scale
funnel(che.mod, studlab = TRUE)
# Funnel plot - OR scale
funnel.rma(che.mod, atransf = exp, layout = "RevMan5", xlab = "Observed Outcome (Odds-Ratio)", label = FALSE) 


## Model checking

# Assess Influence of observations and studies
cd = cooks.distance(che.mod)
plot(cd, type = 'o', pch = 19, xlab = "Observed Outcome", ylab = "Cook's Distance")

cd_study = cooks.distance(che.mod, cluster = Paper)
plot(cd_study, type = 'o', pch = 19, xlab = "Study", ylab = "Cook's Distance")

# Re-fit model without points and studies that seem influential from plots
# Data without Lu et al.
data_no_lu = dat[-c(1,2),]

# Data without outcome 7
data_no_7 = dat[-7,]

# Data without outcome 8
data_no_8 = dat[-8,]

# Without the Lu Paper
V3 <- with(data_no_lu, 
           impute_covariance_matrix(vi = var,
                                    cluster = Paper,
                                    r = rho))
che.mod3 <- rma.mv(yi = EffectSize,
                   V = V3,
                   slab = Paper,
                   data = data_no_lu,
                   random = ~1 | Paper/id,
                   method = 'REML',
                   test = 't',
                   dfs="contain")
summary(che.mod3)

CI3 <- conf_int(che.mod3, 
                vcov = "CR2")
p.val3 <- coef_test(che.mod3, 
                    vcov = "CR2")

# Without outcome 7
V4 <- with(data_no_7, 
           impute_covariance_matrix(vi = var,
                                    cluster = Paper,
                                    r = rho))
che.mod4 <- rma.mv(yi = EffectSize,
                   V = V4,
                   slab = Paper,
                   data = data_no_7,
                   random = ~1 | Paper/id,
                   method = 'REML',
                   test = 't',
                   dfs="contain")
summary(che.mod4)

CI4 <- conf_int(che.mod4, 
                vcov = "CR2")
p.val4 <- coef_test(che.mod4, 
                    vcov = "CR2")


# Without outcome 8
V5 <- with(data_no_8, 
           impute_covariance_matrix(vi = var,
                                    cluster = Paper,
                                    r = rho))
che.mod5 <- rma.mv(yi = EffectSize,
                   V = V5,
                   slab = Paper,
                   data = data_no_8,
                   random = ~1 | Paper/id,
                   method = 'REML',
                   test = 't',
                   dfs="contain")
summary(che.mod5)

CI5 <- conf_int(che.mod5, 
                vcov = "CR2")
p.val5 <- coef_test(che.mod2, 
                    vcov = "CR2")


# Assess Leverage of observations

hat = hatvalues(che.mod)
plot(hat, type = 'o')

# Residuals
std_res = rstandard(che.mod, cluster = Paper)
stdnt_res = rstudent(che.mod, progbar=FALSE, cluster = Paper,
                     reestimate=TRUE, parallel="no", ncpus=1)


## Variance contributions - can run Function for variance component calculation (below)
# or install Xcode, then devtools, then dmetar through devtools
i2 <- var.comp(mod1)
#summary(i2) - can use this if dmetar is installed
i2$results
i2$plot

# Function for variance component calculation

mlm.variance.distribution = var.comp = function(x){
  
  m = x
  
  # Check class
  if (!(class(m)[1] %in% c("rma.mv", "rma"))){
    stop("x must be of class 'rma.mv'.")
  }
  
  # Check for three level model
  if (m$sigma2s != 2){
    stop("The model you provided does not seem to be a three-level model. This function can only be used for three-level models.")
  }
  
  # Check for right specification (nested model)
  if (sum(grepl("/", as.character(m$random[[1]]))) < 1){
    stop("Model must contain nested random effects. Did you use the '~ 1 | cluster/effect-within-cluster' notation in 'random'? See ?metafor::rma.mv for more details.")
  }
  
  # Get variance diagonal and calculate total variance
  n = m$k.eff
  vector.inv.var = 1/(diag(m$V))
  sum.inv.var = sum(vector.inv.var)
  sum.sq.inv.var = (sum.inv.var)^2
  vector.inv.var.sq = 1/(diag(m$V)^2)
  sum.inv.var.sq = sum(vector.inv.var.sq)
  num = (n-1)*sum.inv.var
  den = sum.sq.inv.var - sum.inv.var.sq
  est.samp.var = num/den
  
  # Calculate variance proportions
  level1=((est.samp.var)/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level2=((m$sigma2[2])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level3=((m$sigma2[1])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  
  # Prepare df for return
  Level=c("Level 1", "Level 2", "Level 3")
  Variance=c(level1, level2, level3)
  df.res=data.frame(Variance)
  colnames(df.res) = c("% of total variance")
  rownames(df.res) = Level
  I2 = c("---", round(Variance[2:3], 2))
  df.res = as.data.frame(cbind(df.res, I2))
  
  totalI2 = Variance[2] + Variance[3]
  
  
  # Generate plot
  df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                   "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                   "Type" = rep(1,2))
  
  df2 = data.frame("Level" = rownames(df.res),
                   "Variance" = df.res[,1],
                   "Type" = rep(2,3))
  
  df = as.data.frame(rbind(df1, df2))
  
  
  g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
    coord_cartesian(ylim = c(0,1), clip = "off") +
    geom_bar(stat="identity", position="fill", width = 1, color="black") +
    scale_y_continuous(labels = scales::percent)+
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.line.y = element_blank(),
          axis.title.y=element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_line(lineend = "round"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.background = element_rect(linetype="solid",
                                           colour ="black"),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"cm"),
          axis.ticks.length=unit(.25, "cm"),
          plot.margin = unit(c(1,3,1,1), "lines")) +
    scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                                 "deepskyblue1", "deepskyblue2")) +
    
    # Add Annotation
    
    # Total Variance
    annotate("text", x = 1.5, y = 1.05,
             label = paste("Total Variance:",
                           round(m$sigma2[1]+m$sigma2[2]+est.samp.var, 3))) +
    
    # Sampling Error
    annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
             label = paste("Sampling Error Variance: \n", round(est.samp.var, 3)), size = 3) +
    
    # Total I2
    annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
             label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
    annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
             label = paste("Variance not attributable \n to sampling error: \n", round(m$sigma2[1]+m$sigma2[2],3)), size = 3) +
    
    # Level 1
    annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                       round(df$Variance[3],2), "%", sep=""), size = 3) +
    
    # Level 2
    annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
             label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
    
    # Level 3
    annotate("text", x = 2, y = (df[5,2]/2)/100,
             label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)
  
  returnlist = list(results = df.res,
                    totalI2 = totalI2,
                    plot = g)
  class(returnlist) = c("mlm.variance.distribution", "list")
  
  invisible(returnlist)
  
  returnlist
  
}


## Sensitivity analysis for different values of rho
rho = 0.2
V <- with(dat, 
          impute_covariance_matrix(vi = dat$var,
                                   cluster = Paper,
                                   r = rho))
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't',)
summary(che.mod)

conf_int(che.mod, 
         vcov = "CR2")
coef_test(che.mod, 
          vcov = "CR2")
rho = 0.4
V <- with(dat, 
          impute_covariance_matrix(vi = dat$var,
                                   cluster = Paper,
                                   r = rho))
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't',
                  dfs="contain")
summary(che.mod)

conf_int(che.mod, 
         vcov = "CR2")
coef_test(che.mod, 
          vcov = "CR2")
rho = 0.5
V <- with(dat, 
          impute_covariance_matrix(vi = dat$var,
                                   cluster = Paper,
                                   r = rho))
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't')
summary(che.mod)

conf_int(che.mod, 
         vcov = "CR2")
coef_test(che.mod, 
          vcov = "CR2")
rho = 0.6
V <- with(dat, 
          impute_covariance_matrix(vi = dat$var,
                                   cluster = Paper,
                                   r = rho))
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't')
summary(che.mod)

conf_int(che.mod, 
         vcov = "CR2")
coef_test(che.mod, 
          vcov = "CR2")
rho = 0.8
V <- with(dat, 
          impute_covariance_matrix(vi = dat$var,
                                   cluster = Paper,
                                   r = rho))
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't')
summary(che.mod)

conf_int(che.mod, 
         vcov = "CR2")
coef_test(che.mod, 
          vcov = "CR2")


rho = 1
V <- with(dat, 
          impute_covariance_matrix(vi = dat$var,
                                   cluster = Paper,
                                   r = rho))
V = make.positive.definite(V) 
che.mod <- rma.mv(yi = EffectSize,
                  V = V,
                  slab = Paper,
                  data = dat,
                  random = ~1 | Paper/id,
                  method = 'REML',
                  test = 't')
summary(che.mod)

conf_int(che.mod, 
         vcov = "CR2")
coef_test(che.mod, 
          vcov = "CR2")



### BAYESIAN APPROACH ###

## Plot of possible prior choices

# HC(0, 1), HC(0, 0.5), HC(0, 0.1)
qhcauchy(0.95, 1) # Imposes the least restrictions
qhcauchy(0.95, 0.5)
qhcauchy(0.95, 0.1) # Imposes the most restrictions


x = rhcauchy(1e5, 2)
curve(dhcauchy(x, 1), 0, 15, col = "red", xlim = c(0, 15), xlab = "", ylab = "Density", pty = 2)
curve(dhcauchy(x, 0.5), 0, 15, col = "black", add = TRUE, pty = 2)
curve(dhcauchy(x, 0.1), 0, 15, col = "blue", add = TRUE, pty = 2)
legend("topright", legend=c("HC(0, 1)", "HC(0, 0.5)", "HC(0, 0.1)"),
       col=c("red", "black", "blue"), lty=1, cex=0.8)


## Set the priors for the effect size and the between-study and between-outcome heterogeneity
priors2 <- c(prior(normal(0,1), class = Intercept),
            prior(cauchy(0,0.5), class = sd))

## Fit the model with HC(0,0.5)
iters = 10000
burn = 5000

m.brm <- brm(EffectSize|se(StdErr) ~ 1 + (1|Paper) + (1 | id),
             data = dat,
             prior = priors2,
             iter = iters, warmup = burn, cores = 4, chains = 4,
             seed = 14)

m.brm


## Sensitivity analysis for different scale parameters in half cauchy prior

iters = 10000
burn = 5000

priors.sens1 = c(prior(normal(0,1), class = Intercept),
                 prior(cauchy(0,0.1), class = sd))

m.sens1 <- brm(EffectSize|se(StdErr) ~ 1 + (1|Paper) + (1 | id),
              data = dat,
              prior = priors.sens1,
              iter = iters, warmup = burn, cores = 6, chains = 3,
              seed = 14)

priors.sens2 = c(prior(normal(0,1), class = Intercept),
                 prior(cauchy(0,0.3), class = sd))

m.sens2 <- brm(EffectSize|se(StdErr) ~ 1 + (1|Paper) + (1 | id),
               data = dat,
               prior = priors.sens2,
               iter = iters, warmup = burn, cores = 6, chains = 3,
               seed = 14)


priors.sens3 = c(prior(normal(0,1), class = Intercept),
                 prior(cauchy(0,0.5), class = sd))

m.sens3 <- brm(EffectSize|se(StdErr) ~ 1 + (1|Paper) + (1 | id),
               data = dat,
               prior = priors.sens3,
               iter = iters, warmup = burn, cores = 6, chains = 3,
               seed = 14)


priors.sens4 = c(prior(normal(0,1), class = Intercept),
                 prior(cauchy(0,0.7), class = sd))

m.sens4 <- brm(EffectSize|se(StdErr) ~ 1 + (1|Paper) + (1 | id),
               data = dat,
               prior = priors.sens4,
               iter = iters, warmup = burn, cores = 6, chains = 3,
               seed = 14)


priors.sens5 = c(prior(normal(0,1), class = Intercept),
                 prior(cauchy(0,1), class = sd))

m.sens5 <- brm(EffectSize|se(StdErr) ~ 1 + (1|Paper) + (1 | id),
               data = dat,
               prior = priors.sens5,
               iter = iters, warmup = burn, cores = 6, chains = 3,
               seed = 14)

m.sens1
m.sens2
m.sens3
m.sens4
m.sens5

## Validation - checking convergence

pp_check(m.brm)

# Summary plots of results

post.samples <- posterior_samples(m.brm, c("^b", "^sd"))
names(post.samples)
names(post.samples) <- c("pooled_effect", "outcome_heterogeneity", "study_heterogeneity")
# Transform posterior samples to OR scale
post.samples$pooled_effect = exp(post.samples$pooled_effect)


# Posterior plot for overall pooled effect size
ggplot(aes(x = pooled_effect), data = post.samples) +
  geom_vline(xintercept = 1, color = "gray", size = 1, linetype = 2) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samples$pooled_effect)) +
  labs(x = "Overall Pooled Effect",
       y = element_blank()) +
  xlim(0, 3) +
  theme_minimal()

# Posterior plot for heterogeneity between outcomes
ggplot(aes(x = outcome_heterogeneity), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post.samples$outcome_heterogeneity)) +        # add point at mean
  labs(x = "Heterogeneity between outcomes",
       y = element_blank()) +
  xlim(0, 1) +
  theme_minimal()

# Posterior plot for heterogeneity between studies
ggplot(aes(x = study_heterogeneity), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post.samples$study_heterogeneity)) +        # add point at mean
  labs(x = "Heterogeneity between studies",
       y = element_blank()) +
  xlim(-0.25, 1.5) +
  theme_minimal()


# Forest Plots
papers = m.brm %>%
  spread_draws(b_Intercept, r_Paper[Paper,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = exp(b_Intercept) + r_Paper)

papers$rank = papers$mu
pooledp = papers
pooledp$Paper = "Pooled Effect"
pooledp$rank = -100

ids = m.brm %>%
  spread_draws(b_Intercept, r_id[id,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = exp(b_Intercept) + r_id)

ids$rank = ids$mu
pooledi = ids
pooledi$id = "Pooled Effect"
ids$id = as.character(ids$id)
pooledi$rank = -100

together_id = bind_rows(ids, pooledi)
together_paper = bind_rows(papers, pooledp)

# Forest plot OR scale for outcomes
together_id %>%
  ungroup() %>%
  
  # plot
  ggplot(aes(x = mu,  y = reorder(id, rank))) +
  geom_vline(xintercept = 1, color = "gray", size = 1, linetype = 2) +
  geom_vline(xintercept = exp(fixef(m.brm)[1, 1]), color = "white", size = 1) +
  geom_vline(xintercept = exp(fixef(m.brm)[1, 3:4]), color = "white", linetype = 2) +
  stat_halfeye(.width = .95, size = 2/3) +
  labs(x = expression(italic("Odds Ratio")),
       y = NULL) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0))


# Forest plot OR scale for papers
together_paper %>%
  ungroup() %>%
  
  # plot
  ggplot(aes(x = mu,  y = reorder(Paper, rank))) +
  geom_vline(xintercept = 1, color = "gray", size = 1, linetype = 2) +
  geom_vline(xintercept = exp(fixef(m.brm)[1, 1]), color = "white", size = 1) +
  geom_vline(xintercept = exp(fixef(m.brm)[1, 3:4]), color = "white", linetype = 2) +
  stat_halfeye(.width = .95, size = 2/3) +
  labs(x = expression(italic("Odds Ratio")),
       y = NULL) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0))
