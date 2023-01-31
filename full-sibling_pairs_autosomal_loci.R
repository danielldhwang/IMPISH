# This is the R scripts used to imputed parental genotypes at autosomal loci using genotypes of offspring full sibling pairs
# Please refer to Hwang et al. (2020) PLOS GENETICS

library(lmerTest)
library(MASS)
rm(list=ls())

Nrep = 1000	#Number of simulation replicates
N = 2000	#Sample size
Vm = 0.001	#Variance explained by maternal effect
Vf = 0.001	#Variance explained by fetal effect
p = 0.1		#Increaser allele frequency
Opp = FALSE	#Boolean variable to denote whether maternal and fetal effects in opposite directions (TRUE = opposite directions)
rho = 0.2	#Covariance between sibling residual variances
q <- 1-p  	#Decreaser allele frequency

#Create vectors to save results from regression analyses
#For analysis using simulated parental and fetal genotypes
beta_genotyped_mat <- vector(length = Nrep)		#Fitting regression with maternal genotypes known (maternal coefficient)
beta_genotyped_fet <- vector(length = Nrep)		#Fitting regression with maternal genotypes known (fetal coefficient)
se_genotyped_mat <- vector(length = Nrep)
se_genotyped_fet <- vector(length = Nrep)
pval_genotyped_mat <- vector(length = Nrep)
pval_genotyped_fet <- vector(length = Nrep)

#Analysis using imputed parental genotypes and simulated fetal genotypes
beta_imputed_mat <- vector(length = Nrep)			#Fitting regression with genotypes imputed (maternal coefficient)
beta_imputed_fet <- vector(length = Nrep)			#Fitting regression with genotypes imputed (fetal coefficient)
se_imputed_mat <- vector(length = Nrep)
se_imputed_fet <- vector(length = Nrep)
pval_imputed_mat <- vector(length = Nrep)
pval_imputed_fet <- vector(length = Nrep)

#Analysis using only sib genotypes
beta_sib <- vector(length = Nrep)				#Fitting incorrect model to just sibs
se_sib <- vector(length = Nrep)
pval_sib <- vector(length = Nrep)

#ANOVA for full model against null model
pval_imputed_omni <- vector(length = Nrep)
pval_genotyped_omni <- vector(length = Nrep)

a <- sqrt(1/(2*p*q))   	#Create genetic variable of variance one. Assume no dominance. 
BetaM <- sqrt(Vm) 	#Path coefficient for maternal effect
BetaF <- sqrt(Vf) 	#Path coefficient for fetal effect
if(Opp == TRUE) {  	#If Opp is true then maternal and fetal effects have opposite directions of effect
  BetaF = -BetaF
}

Ve <- (1 - BetaM^2 - BetaF^2 - 2*0.5*BetaM*BetaF) #Residual variance in trait

for(j in 1:Nrep) {
  #Sample mothers' genotypes
  Zm <- sample(x = c("AA","Aa","aa"), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  #Sample fathers' genotypes
  Zp <- sample(x = c("AA","Aa","aa"), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  
  Z_sib1 <- vector(length = N)
  Z_sib2 <- vector(length = N)
  Z_imp <- vector(length = N)		#Imputed vector of genotypes at mother's locus
  
  #Simulate sib 1 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sib1[i] <- switch(paste(Zm[i],Zp[i]),
                        'AA AA'='AA',
                        'AA Aa'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                        'AA aa'='Aa',
                        'Aa AA'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                        'Aa Aa'=ifelse(r[i] <= 0.25, 'aa', ifelse(r[i] > 0.75, 'AA', 'Aa')),
                        'Aa aa'=ifelse(r[i] <= 0.5, 'aa', 'Aa'),
                        'aa AA'='Aa',
                        'aa Aa'=ifelse(r[i] <= 0.5, 'aa', 'Aa'),
                        'aa aa'='aa')
  }
  
  #Simulate sib 2 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sib2[i] <- switch(paste(Zm[i], Zp[i]),
                        'AA AA'='AA',
                        'AA Aa'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                        'AA aa'='Aa',
                        'Aa AA'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                        'Aa Aa'=ifelse(r[i] <= 0.25, 'aa', ifelse(r[i] > 0.75, 'AA', 'Aa')),
                        'Aa aa'=ifelse(r[i] <= 0.5, 'aa', 'Aa'),
                        'aa AA'='Aa',
                        'aa Aa'=ifelse(r[i] <= 0.5, 'aa', 'Aa'),
                        'aa aa'='aa')
  }
  
  #Change coding to 0, 1 and 2 in order to calculate allele frequencies
  Z_sib1_012 <- ifelse(Z_sib1=='AA', 0, ifelse(Z_sib1=='Aa', 1, 2))
  Z_sib2_012 <- ifelse(Z_sib2=='AA', 0, ifelse(Z_sib2=='Aa', 1, 2))
  
  #Calculate observed allele frequency based on sibling genotypes
  q_obs <- sum(Z_sib1_012+Z_sib2_012)/(2*N*2)
  p_obs <- 1 - q_obs
  
  #Impute mother's genotypes
  #Calculate conditional probability of mother's genotypes given sibs genotypes
  p_AA_AA_AA = p_obs*(2*p_obs+2)/(p_obs+1)^2		#P(Mum = AA | Sib = AA and Sib = AA) 
  p_Aa_AA_AA = 2/(p_obs+1) - 1					#P(Mum = Aa | Sib = AA and Sib = AA)
  p_aa_AA_AA = 0							#P(Mum = aa | Sib = AA and Sib = AA)
  p_AA_AA_Aa = p_obs/(p_obs+1)					#P(Mum = AA | Sib = AA and Sib = Aa)
  p_Aa_AA_Aa = 1/(p_obs+1)					#P(Mum = Aa | Sib = AA and Sib = Aa)
  p_aa_AA_Aa = 0							#P(Mum = aa | Sib = AA and Sib = Aa)
  p_AA_AA_aa = 0							#P(Mum = AA | Sib = AA and Sib = aa)
  p_Aa_AA_aa = 1							#P(Mum = Aa | Sib = AA and Sib = aa)
  p_aa_AA_aa = 0							#P(Mum = aa | Sib = AA and Sib = aa)
  p_AA_Aa_Aa = (0.5*(p_obs-2)*p_obs)/(p_obs^2-p_obs-1)	#P(Mum = AA | Sib = Aa and Sib = Aa)
  p_Aa_Aa_Aa = -0.5 / ((p_obs-1)*p_obs-1)				#P(Mum = Aa | Sib = Aa and Sib = Aa)
  p_aa_Aa_Aa = (0.5*(p_obs-1)*(p_obs+1))/(p_obs^2-p_obs-1)	#P(Mum = aa | Sib = Aa and Sib = Aa)
  p_AA_aa_Aa = 0							#P(Mum = AA | Sib = aa and Sib = Aa)
  p_Aa_aa_Aa = 1/(q_obs+1)						#P(Mum = Aa | Sib = aa and Sib = Aa)
  p_aa_aa_Aa = q_obs/(q_obs+1)					#P(Mum = aa | Sib = aa and Sib = Aa)
  p_AA_aa_aa = 0							#P(Mum = AA | Sib = aa and Sib = aa) 
  p_Aa_aa_aa = 2/(q_obs+1) - 1					#P(Mum = Aa | Sib = aa and Sib = aa)
  p_aa_aa_aa = q_obs*(2*q_obs+2)/(q_obs+1)^2			#P(Mum = aa | Sib = aa and Sib = aa)
  
  for (i in 1:N) {
    Z_imp[i] <- switch(paste(Z_sib1[i], Z_sib2[i]),
                       'AA AA'= -a*p_AA_AA_AA + 0*p_Aa_AA_AA + a*p_aa_AA_AA,
                       'AA Aa'= -a*p_AA_AA_Aa + 0*p_Aa_AA_Aa + a*p_aa_AA_Aa,
                       'AA aa'= -a*p_AA_AA_aa + 0*p_Aa_AA_aa + a*p_aa_AA_aa,
                       'Aa AA'= -a*p_AA_AA_Aa + 0*p_Aa_AA_Aa + a*p_aa_AA_Aa,
                       'Aa Aa'= -a*p_AA_Aa_Aa + 0*p_Aa_Aa_Aa + a*p_aa_Aa_Aa,
                       'Aa aa'= -a*p_AA_aa_Aa + 0*p_Aa_aa_Aa + a*p_aa_aa_Aa,
                       'aa AA'= -a*p_AA_AA_aa + 0*p_Aa_AA_aa + a*p_aa_AA_aa,
                       'aa Aa'= -a*p_AA_aa_Aa + 0*p_Aa_aa_Aa + a*p_aa_aa_Aa,
                       'aa aa'= -a*p_AA_aa_aa + 0*p_Aa_aa_aa + a*p_aa_aa_aa)
  }
  
  #Change coding AA/Aa/aa to the genetic value -a/0/a
  Z_sib1 <- ifelse(Z_sib1=='AA', -a, ifelse(Z_sib1=='Aa', 0, a))
  Z_sib2 <- ifelse(Z_sib2=='AA', -a, ifelse(Z_sib2=='Aa', 0, a))
  Zm <- ifelse(Zm=='AA', -a, ifelse(Zm=='Aa', 0, a))
  Zp <- ifelse(Zp=='AA', -a, ifelse(Zp=='Aa', 0, a))
  
  #Create correlated error variables for sib 1 and sib 2
  Sigma <- matrix(c(Ve , rho, rho, Ve),2,2)
  e <- mvrnorm(n = N, mu = c(0, 0), Sigma)
  
  #Simulate offspring outcome
  Y_sib1 <- BetaM*Zm + BetaF*Z_sib1 + e[,1]
  Y_sib2 <- BetaM*Zm + BetaF*Z_sib2 + e[,2]
  
  #Create a test dataset of simulated and imputed genotypes and phenotype (offspring outcome)
  test <- data.frame(rbind(cbind(1:N, Y_sib1, Zm, Z_imp, Z_sib1), 
                           cbind(1:N, Y_sib2, Zm, Z_imp, Z_sib2)))
  colnames(test) <- c("fam", "Y","Zm", "Z_imp", "Z")
  
  #Run analyses
  lmer_genotyped <- lmer(Y ~ Zm + Z + (1|fam), data = test, REML = FALSE);  results_genotyped <- summary(lmer_genotyped)
  lmer_imputed <- lmer(Y ~ Z_imp + Z + (1|fam), data = test, REML = FALSE); results_imputed <- summary(lmer_imputed)
  lmer_sib <- lmer(Y ~ Z + (1|fam), data = test, REML = FALSE); results_sib <- summary(lmer_sib)
  lmer_null <- lmer(Y ~ (1|fam), data = test, REML = FALSE)
  
  beta_genotyped_mat[j] <- results_genotyped$coefficient[2,1]
  se_genotyped_mat[j] <- results_genotyped$coefficient[2,2]
  pval_genotyped_mat[j] <- results_genotyped$coefficient[2,5]
  
  beta_genotyped_fet[j] <- results_genotyped$coefficient[3,1]
  se_genotyped_fet[j] <- results_genotyped$coefficient[3,2]
  pval_genotyped_fet[j] <- results_genotyped$coefficient[3,5]
  
  pval_genotyped_omni[j] <- anova(lmer_genotyped,lmer_null)[2,8] 
  
  beta_imputed_mat[j] <- results_imputed$coefficient[2,1]
  se_imputed_mat[j] <- results_imputed$coefficient[2,2]
  pval_imputed_mat[j] <- results_imputed$coefficient[2,5]
  
  beta_imputed_fet[j] <- results_imputed$coefficient[3,1]
  se_imputed_fet[j] <- results_imputed$coefficient[3,2]
  pval_imputed_fet[j] <- results_imputed$coefficient[3,5]
  
  pval_imputed_omni[j] <- anova(lmer_imputed,lmer_null)[2,8] 
  
  beta_sib[j] <- results_sib$coefficient[2,1]
  se_sib[j] <- results_sib$coefficient[2,2]
  pval_sib[j] <- results_sib$coefficient[2,5]
}