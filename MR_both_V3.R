library(MendelianRandomization)
library(tidyverse)
library(devtools)
#install_github("phenoscanner/phenoscanner")
library(phenoscanner)
library(MRInstruments) 
library(TwoSampleMR)
library(readr)
library(LDlinkR)

### Load in the data
### Using GWAS Atlas data - 
### cross reference with the MR Instruments data to fill in the gaps

# These are the 77 SNPS we need from the Obesity study
BMI_SNPS <- read_csv("topSNPtable_atlasID_143.csv")
BMI_SNPS = as.data.frame(BMI_SNPS)
nd = read.table("Noelle_data.uniq", header = T)
dat = nd[which(nd$SNP %in% BMI_SNPS$rsID), ]
rs9581854 = nd[which(nd$SNP == 'rs12016871'), ]
rs9581854$SNP = 'rs9581854'

dat = rbind(dat,rs9581854)
dat = dat[order(dat[,1],decreasing = TRUE),]
BMI_SNPS = BMI_SNPS[order(BMI_SNPS[,3],decreasing = TRUE),]

bmi_exp_dat = cbind(dat,BMI_SNPS$CHR,rep('bmi',77))

colnames(bmi_exp_dat) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "N",'chr', 'phenotype')

# Now we have extra info about 74 of the 77 SNPs (we will drop the remaining three due to missing data)
# This is exactly the same data (although the p-vals are rounded slightly in the bmi_gwas data, which is fine)
bmi_exp_dat <- format_data(dat, phenotype_col = 'phenotype')
bmi_exp_dat
### Check for linkage disequilibrium on each chromosome
chr1 <- which(bmi_exp_dat$chr.exposure == '1')
chr1 <- bmi_exp_dat$SNP[chr1]
LDmatrix(chr1, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr2 <- which(bmi_exp_dat$chr.exposure == '2')
chr2 <- bmi_exp_dat$SNP[chr2]
LDmatrix(chr2, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr3 <- which(bmi_exp_dat$chr.exposure == '3')
chr3 <- exposure_dat$SNP[chr3]
LDmatrix(chr3, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr4 <- which(bmi_exp_dat$chr.exposure == '4')
chr4 <- bmi_exp_dat$SNP[chr4]
LDmatrix(chr4, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr5 <- which(bmi_exp_dat$chr.exposure == '5')
chr5 <- bmi_exp_dat$SNP[chr5]
LDmatrix(chr5, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr6 <- which(bmi_exp_dat$chr.exposure == '6')
chr6 <- bmi_exp_dat$SNP[chr6]
LDmatrix(chr6, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr7 <- which(bmi_exp_dat$chr.exposure == '7')
chr7 <- bmi_exp_dat$SNP[chr7]
LDmatrix(chr4, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr8 <- which(bmi_exp_dat$chr.exposure == '8')
chr8 <- bmi_exp_dat$SNP[chr8]
LDmatrix(chr8, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr9 <- which(bmi_exp_dat$chr.exposure == '9')
chr9 <- bmi_exp_dat$SNP[chr9]
LDmatrix(chr9, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr10 <- which(bmi_exp_dat$chr.exposure == '10')
chr10 <- bmi_exp_dat$SNP[chr10]
LDmatrix(chr10, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr11 <- which(bmi_exp_dat$chr.exposure == '11')
chr11 <- bmi_exp_dat$SNP[chr11]
LDmatrix(chr11, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr12 <- which(bmi_exp_dat$chr.exposure == '12')
chr12 <- bmi_exp_dat$SNP[chr12]
LDmatrix(chr12, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr13 <- which(bmi_exp_dat$chr.exposure == '13')
chr13 <- bmi_exp_dat$SNP[chr13]
LDmatrix(chr13, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr14 <- which(bmi_exp_dat$chr.exposure == '14')
chr14 <- bmi_exp_dat$SNP[chr14]
LDmatrix(chr14, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr15 <- which(bmi_exp_dat$chr.exposure == '15')
chr15 <- bmi_exp_dat$SNP[chr15]
LDmatrix(chr15, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr16 <- which(bmi_exp_dat$chr.exposure == '16')
chr16 <- bmi_exp_dat$SNP[chr16]
LDmatrix(chr16, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr17 <- which(bmi_exp_dat$chr.exposure == '17')
chr17 <- bmi_exp_dat$SNP[chr17]
LDmatrix(chr17, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr18 <- which(bmi_exp_dat$chr.exposure == '18')
chr18 <- bmi_exp_dat$SNP[chr18]
LDmatrix(chr18, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

chr19 <- which(bmi_exp_dat$chr.exposure == '19')
chr19 <- bmi_exp_dat$SNP[chr19]
LDmatrix(chr19, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")

LDmatrix(bmi_exp_dat$SNP, 
         pop = "CEU", 
         r2d = "r2", 
         token = "b22ecb07e3d4", 
         file = FALSE, 
         genome_build = "grch37")
# These show no linkage disequilibrium in out 77 SNPs


# Get the outcome data 
asthma_out_dat <- extract_outcome_data(
  snps = bmi_exp_dat$SNP,
  outcomes = 'ebi-a-GCST90014325', # Study on asthma
  proxies = TRUE, # Use a proxy if we cannot find the SNP from the exposure data in the outcome data
  rsq = 0.7 # This is the r2 value from the study
)

# Harmonise the datasets - Try to infer the forward strand alleles using allele frequency information (set action = 2)
bmi_asthma_dat = harmonise_data(bmi_exp_dat, asthma_out_dat, action = 2)
bmi_asthma_dat$exposure = "Body mass index"
bmi_asthma_dat$outcome = "Asthma"

### Conducting the MR using the Mendelian Randomization package

# The values that the TwoSampleMR package actually uses for MR
bmi_asthma_dat_keep = bmi_asthma_dat[which(bmi_asthma_dat$mr_keep == T),]

MRInputObject = mr_input(bx = bmi_asthma_dat_keep$beta.exposure, 
                         bxse = bmi_asthma_dat_keep$se.exposure,
                         by = bmi_asthma_dat_keep$beta.outcome,
                         byse = bmi_asthma_dat_keep$se.outcome,
                         exposure = "Obesity: Body mass index",
                         outcome = "Asthma",
                         snps = bmi_asthma_dat_keep$SNP)

# Conduct the MR
MRAllObject_all <- mr_allmethods(MRInputObject, method = "all")
MRAllObject_all

## Plots

# Scatter plot
mr_plot(MRInputObject)

mr_plot(MRAllObject_all, orientate = F, interactive = F)

# Explanation: Scatter plot created by mr_plot command applied to a MRInput object. Estimated genetic
# associations with the outcome (vertical axis) are plotted against predicted associations with the outcome from the
# inverse-variance weighted method (horizontal axis). Error bars are 95% confidence intervals, and the
# diagonal line has gradient 1.


# Forest plot A - SNPs 
mr_forest(MRInputObject)

# Forest plot B - Methods
mr_forest(MRInputObject,
          snp_estimates=FALSE,
          methods = c("ivw", "median", "wmedian", "egger", "maxlik", "mbe", "conmix"))

# Forest plot explanations: Forest plots created by mr_forest command. A: comparison of variant-specific estimates
# plus inverse-variance weighted (IVW) estimate (default options). B: comparison of estimates from different
# methods with variant-specific estimates switched off. Points represent estimates and horizontal error bars are 95%
# confidence intervals (CI). These are all log-odds

# Funnel plot
mr_funnel(MRInputObject)

# Explanation: Funnel plot created by mr_funnel command. Points represent variant-specific estimates and horizontal
# error bars are 95% confidence intervals (CI).

# Leave-one-out plot
mr_loo(MRInputObject)

# Explanation: Leave-one-out plot created by mr_loo command. Points represent estimates from the inverse-variance
# weighted (IVW) method, omitting the variant indicated. Horizontal error bars are 95% confidence intervals (CI).
# The mr_loo function allows the user to investigate sensitivity of the IVW estimate to individual data points. This
# is implemented by calculating the IVW estimate omitting each variant from the analysis in turn (loo stands for
# "leave one out"). The IVW estimate based on all the variants is also plotted for reference.



### Conducting the MR using TwoSampleMR package

# Perform the MR
res <- mr(bmi_asthma_dat)
res

## Sensitivity analyses

# Heterogeneity
mr_heterogeneity(bmi_asthma_dat)

# Horizontal pleiotropy
mr_pleiotropy_test(bmi_asthma_dat)

# Single SNP analysis
# The method used to perform the single SNP MR is the Wald ratio by default
res_single <- mr_singlesnp(bmi_asthma_dat)

# Leave-one-out analysis - It is possible to perform a leave-one-out analysis, where the MR is performed again but leaving out each SNP in turn, 
# to identify if a single SNP is driving the association. It uses the inverse variance weighted method.
res_loo <- mr_leaveoneout(bmi_asthma_dat)

## Plots

p1 <- mr_scatter_plot(res, bmi_asthma_dat)
p1[[1]]

# The scatterplot suggests a positive causal relationship of the SNP effects on BMI against the SNP effects on asthma. Each point displayed on the 
# graph represents a single genetic variant. The horizontal and vertical lines extending from each point represent the 95% confidence interval for the genetic associations. 
# The horizontal axis of the graph displays the estimated genetic associations with the exposure (BMI), and the vertical axis displays the estimated genetic associations 
# with the outcome (asthma).

p2 <- mr_forest_plot(res_single)
p2[[1]]

# The forest plot shows the causal estimate using each SNP alone as well as the overall causal estimate using all the SNPs with MR-Egger and IVW. 
# The error bars represent the 95% confidence intervals.

p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

# Explanation: The estimated causal effect is shown for each excluded SNP and the overall estimate using all the SNPs is shown in red. 
# The error bars represent the 95% confidence intervals.

p4 <- mr_funnel_plot(res_single)
p4[[1]]


