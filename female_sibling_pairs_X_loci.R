# This is the R scripts used to imputed parental genotypes at X chromosomal loci using genotypes of offspring female sibling pairs
# Please refer to Hwang et al. (2020) PLOS GENETICS

library(lmerTest)
library(MASS)
rm(list=ls())

Nrep = 1000	#Number of simulation replicates
N = 2000	#Sample size
Vm = 0.001	#Variance explained by maternal effect
Vp = 0.001	#Variance explained by paternal effect
Vf = 0.001	#Variance explained by fetal effect
p = 0.1		#Increaser allele frequency
Opp = FALSE	#Boolean variable to denote whether maternal and fetal effects in opposite directions (TRUE = opposite directions)
rho = 0.2	#Covariance between sibling residual variances
q <- 1-p 	#Decreaser allele frequency

#Create vectors to save results from regression analyses
#For analysis using simulated genotypes
beta_genotyped_mat <- vector(length = Nrep) #Fitting regression with maternal and paternal genotypes known (maternal coefficient)
beta_genotyped_pat <- vector(length = Nrep)  #Fitting regression with maternal and paternal genotypes known (paternal coefficient)
beta_genotyped_fet <- vector(length = Nrep)   #Fitting regression with maternal and paternal genotypes known (fetal coefficient)
se_genotyped_mat <- vector(length = Nrep)
se_genotyped_pat <- vector(length = Nrep)
se_genotyped_fet <- vector(length = Nrep)
pval_genotyped_mat <- vector(length = Nrep)
pval_genotyped_pat <- vector(length = Nrep)
pval_genotyped_fet <- vector(length = Nrep)

#Analysis using imputed maternal and paternal genotypes
beta_imputed_mat <- vector(length = Nrep)	#Fitting regression with genotypes imputed (maternal coefficient)
beta_imputed_pat <- vector(length = Nrep)	#Fitting regression with genotypes imputed (paternal coefficient)
beta_imputed_fet <- vector(length = Nrep)	#Fitting regression with genotypes imputed (fetal coefficient)
se_imputed_mat <- vector(length = Nrep)
se_imputed_pat <- vector(length = Nrep)
se_imputed_fet <- vector(length = Nrep)
pval_imputed_mat <- vector(length = Nrep)
pval_imputed_pat <- vector(length = Nrep)
pval_imputed_fet <- vector(length = Nrep)

#Analysis using only sibling genotypes
beta_sib <- vector(length = Nrep)		#Fitting incorrect model including only siblings
se_sib <- vector(length = Nrep)
pval_sib <- vector(length = Nrep)

#ANOVA for full model against null model
pval_genotyped_omni <- vector(length = Nrep)
pval_imputed_omni <- vector(length = Nrep)

a <- sqrt(1/(2*p*q)) #Create genetic variable of variance one for maternal and fetal SNPs on X chromosome. Assume no dominance.
BetaM <- sqrt(Vm) #Path coefficient for maternal effect
BetaP <- sqrt(Vp/2) #Path coefficient for paternal effect
BetaF <- sqrt(Vf) #Path coefficient for fetal effect
if(Opp == TRUE) {  #If Opp is true then maternal and fetal effects have opposite directions of effect
  BetaF = -BetaF
}

Ve_x <- (1 - BetaM^2 - BetaF^2 - BetaP^2*2 - 2*0.5*BetaM*BetaF - 2*0.5*BetaP*BetaF*2) #Residual variance in trait for female sib

for(j in 1:Nrep) {
  #Sample mothers' genotypes
  Zm <- sample(x = c('AA','Aa','aa'), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  #Sample fathers' genotypes
  Zp <- sample(x = c('A','a'), size = N, replace = TRUE, prob = c(p, q))
  
  Z_sib1x <- vector(length = N)
  Z_sib2x <- vector(length = N)
  
  Zmxx_imp <- vector(length = N)		#Imputed vector of genotypes at mum's locus
  Zpxx_imp <- vector(length = N)		#Imputed vector of genotypes at dad's locus
  
  #Simulate female sib 1 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sib1x[i] <- switch(paste(Zm[i],Zp[i]),
                         'AA A'='AA',
                         'AA a'='Aa',
                         'Aa A'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                         'Aa a'=ifelse(r[i] <= 0.5, 'Aa', 'aa'),
                         'aa A'='Aa',
                         'aa a'='aa')
  }
  
  #Simulate female sib 2 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sib2x[i] <- switch(paste(Zm[i],Zp[i]),
                         'AA A'='AA',
                         'AA a'='Aa',
                         'Aa A'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                         'Aa a'=ifelse(r[i] <= 0.5, 'Aa', 'aa'),
                         'aa A'='Aa',
                         'aa a'='aa')
  }
  
  #Change coding to 0, 1 and 2 in order to calculate allele frequencies
  Z_sib1x_012 <- ifelse(Z_sib1x=='AA', 0,
                        ifelse(Z_sib1x=='Aa', 1, 2))
  Z_sib2x_012 <- ifelse(Z_sib2x=='AA', 0,
                        ifelse(Z_sib2x=='Aa', 1, 2))
  
  #Calculate observed allele frequency based on sibling genotypes
  q_obs <- sum(Z_sib1x_012, Z_sib2x_012)/(2*N*2)
  p_obs <- 1 - q_obs
  
  #Impute parental genotypes from female sibs
  #Calculate conditional probability of mother's genotypes given sibs genotypes
  p_AA_AA_AA = (2*p_obs)/(p_obs+1)		#P(Mum = AA | Sibx = AA and Sibx = AA) 
  p_Aa_AA_AA = (1-p_obs)/(p_obs+1)		#P(Mum = Aa | Sibx = AA and Sibx = AA)
  p_aa_AA_AA = 0					#P(Mum = aa | Sibx = AA and Sibx = AA)
  p_AA_AA_Aa = 0			    		#P(Mum = AA | Sibx = AA and Sibx = Aa)
  p_Aa_AA_Aa = 1    					#P(Mum = Aa | Sibx = AA and Sibx = Aa)
  p_aa_AA_Aa = 0					#P(Mum = aa | Sibx = AA and Sibx = Aa)
  p_AA_AA_aa = 0					#P(Mum = AA | Sibx = AA and Sibx = aa)
  p_Aa_AA_aa = 1					#P(Mum = Aa | Sibx = AA and Sibx = aa)
  p_aa_AA_aa = 0					#P(Mum = aa | Sibx = AA and Sibx = aa)
  p_AA_Aa_Aa = 2/3*p_obs				#P(Mum = AA | Sibx = Aa and Sibx = Aa)
  p_Aa_Aa_Aa = 1/3					#P(Mum = Aa | Sibx = Aa and Sibx = Aa)
  p_aa_Aa_Aa = 2/3*q_obs				#P(Mum = aa | Sibx = Aa and Sibx = Aa)
  p_AA_aa_Aa = 0					#P(Mum = AA | Sibx = aa and Sibx = Aa)
  p_Aa_aa_Aa = 1    					#P(Mum = Aa | Sibx = aa and Sibx = Aa)
  p_aa_aa_Aa = 0    					#P(Mum = aa | Sibx = aa and Sibx = Aa)
  p_AA_aa_aa = 0					#P(Mum = AA | Sibx = aa and Sibx = aa) 
  p_Aa_aa_aa = (1-q_obs)/(q_obs+1)			#P(Mum = Aa | Sibx = aa and Sibx = aa)
  p_aa_aa_aa = (2*q_obs)/(q_obs+1)	    		#P(Mum = aa | Sibx = aa and Sibx = aa)
  
  for (i in 1:N) {
    Zmxx_imp[i] <- switch(paste(Z_sib1x[i], Z_sib2x[i]),
                          'AA AA'= -a*p_AA_AA_AA + 0*p_Aa_AA_AA + a*p_aa_AA_AA,
                          'AA Aa'= -a*p_AA_AA_Aa + 0*p_Aa_AA_Aa + a*p_aa_AA_Aa,
                          'AA aa'= -a*p_AA_AA_aa + 0*p_Aa_AA_aa + a*p_aa_AA_aa,
                          'Aa AA'= -a*p_AA_AA_Aa + 0*p_Aa_AA_Aa + a*p_aa_AA_Aa,
                          'Aa Aa'= -a*p_AA_Aa_Aa + 0*p_Aa_Aa_Aa + a*p_aa_Aa_Aa,
                          'Aa aa'= -a*p_AA_aa_Aa + 0*p_Aa_aa_Aa + a*p_aa_aa_Aa,
                          'aa AA'= -a*p_AA_aa_AA + 0*p_Aa_aa_AA + a*p_aa_aa_AA,
                          'aa Aa'= -a*p_AA_aa_Aa + 0*p_Aa_aa_Aa + a*p_aa_aa_Aa,
                          'aa aa'= -a*p_AA_aa_aa + 0*p_Aa_aa_aa + a*p_aa_aa_aa)
  }
  
  #Calculate conditional probability of father's genotypes given sib's genotype
  p_dA_AA_AA <- 1 					#P(Dad = A | Sibx = AA and Sibx = AA) 
  p_da_AA_AA <- 0 					#P(Dad = a | Sibx = AA and Sibx = AA) 
  p_dA_AA_Aa <- 1 					#P(Dad = A | Sibx = AA and Sibx = Aa)
  p_da_AA_Aa <- 0 					#P(Dad = a | Sibx = AA and Sibx = Aa) 
  p_dA_Aa_AA <- 1 					#P(Dad = A | Sibx = Aa and Sibx = AA)
  p_da_Aa_AA <- 0 					#P(Dad = a | Sibx = Aa and Sibx = AA)
  p_dA_Aa_Aa <- 2/3-1/3*p_obs 			#P(Dad = A | Sibx = Aa and Sibx = Aa)
  p_da_Aa_Aa <- 2/3-1/3*q_obs 			#P(Dad = a | Sibx = Aa and Sibx = Aa)
  p_dA_Aa_aa <- 0 					#P(Dad = A | Sibx = Aa and Sibx = aa)
  p_da_Aa_aa <- 1 					#P(Dad = a | Sibx = Aa and Sibx = aa)
  p_dA_aa_Aa <- 0 					#P(Dad = A | Sibx = aa and Sibx = Aa)
  p_da_aa_Aa <- 1 					#P(Dad = a | Sibx = aa and Sibx = Aa)
  p_dA_aa_aa <- 0 					#P(Dad = A | Sibx = aa and Sibx = aa)
  p_da_aa_aa <- 1 					#P(Dad = a | Sibx = aa and Sibx = aa)
  
  for (i in 1:N) {
    Zpxx_imp[i] <- switch(paste(Z_sib1x[i], Z_sib2x[i]),
                          'AA AA'= -a*p_dA_AA_AA + a*p_da_AA_AA,
                          'Aa AA'= -a*p_dA_Aa_AA + a*p_da_Aa_AA,
                          'AA Aa'= -a*p_dA_AA_Aa + a*p_da_AA_Aa,
                          'Aa Aa'= -a*p_dA_Aa_Aa + a*p_da_Aa_Aa,
                          'aa Aa'= -a*p_dA_aa_Aa + a*p_da_aa_Aa,
                          'Aa aa'= -a*p_dA_Aa_aa + a*p_da_Aa_aa,
                          'aa aa'= -a*p_dA_aa_aa + a*p_da_aa_aa)
  }  
  
  #Convert AA/Aa/aa to genetic value 
  Zm <- ifelse(Zm=='AA', -a, ifelse(Zm=='Aa', 0, a))
  Zp <- ifelse(Zp=='A', -a, a)
  Z_sib1x <- ifelse(Z_sib1x=='AA', -a, ifelse(Z_sib1x=='Aa', 0, a))
  Z_sib2x <- ifelse(Z_sib2x=='AA', -a, ifelse(Z_sib2x=='Aa', 0, a))
  
  #Create correlated error variables for sib 1 and sib 2 in the model
  Sigma_xx <- matrix(c(Ve_x, rho, rho, Ve_x),2,2)
  e_xx <- mvrnorm(n = N, mu = c(0, 0), Sigma_xx)   
  
  #Simulate offspring outcome
  Y_sib1x_xx <- BetaM*Zm + BetaP*Zp + BetaF*Z_sib1x + e_xx[,1]
  Y_sib2x_xx <- BetaM*Zm + BetaP*Zp + BetaF*Z_sib2x + e_xx[,2]
  
  test <- data.frame(rbind(cbind(1:N, Y_sib1x_xx, Zm, Zp, Zmxx_imp, Zpxx_imp, Z_sib1x), 
                           cbind(1:N, Y_sib2x_xx, Zm, Zp, Zmxx_imp, Zpxx_imp, Z_sib2x)))
  colnames(test) <- c("fam", "Y_xx", "Zm", "Zp", "Zmxx_imp", "Zpxx_imp", "Z_xx")
  
  #Run analyses
  lmer_genotyped <- lmer(Y_xx ~ Zm + Zp + Z_xx + (1|fam), data = test, REML = FALSE);  results_genotyped <- summary(lmer_genotyped)
  lmer_imputed <- lmer(Y_xx ~ Zmxx_imp + Zpxx_imp + Z_xx + (1|fam), data = test, REML = FALSE); results_imputed <- summary(lmer_imputed)
  lmer_sib <- lmer(Y_xx ~ Z_xx + (1|fam), data = test, REML = FALSE); results_sib <- summary(lmer_sib)
  lmer_null <- lmer(Y_xx ~ (1|fam), data = test, REML = FALSE)
  
  beta_genotyped_mat[j] <- results_genotyped$coefficient[2,1]
  se_genotyped_mat[j] <- results_genotyped$coefficient[2,2]
  pval_genotyped_mat[j] <- results_genotyped$coefficient[2,5]
  
  beta_genotyped_pat[j] <- results_genotyped$coefficient[3,1]
  se_genotyped_pat[j] <- results_genotyped$coefficient[3,2]
  pval_genotyped_pat[j] <- results_genotyped$coefficient[3,5]
  
  beta_genotyped_fet[j] <- results_genotyped$coefficient[4,1]
  se_genotyped_fet[j] <- results_genotyped$coefficient[4,2]
  pval_genotyped_fet[j] <- results_genotyped$coefficient[4,5]
  
  pval_genotyped_omni[j] <- anova(lmer_genotyped,lmer_null)[2,8] 
  
  beta_imputed_mat[j] <- results_imputed$coefficient[2,1]
  se_imputed_mat[j] <- results_imputed$coefficient[2,2]
  pval_imputed_mat[j] <- results_imputed$coefficient[2,5]
  
  beta_imputed_pat[j] <- results_imputed$coefficient[3,1]
  se_imputed_pat[j] <- results_imputed$coefficient[3,2]
  pval_imputed_pat[j] <- results_imputed$coefficient[3,5]
  
  beta_imputed_fet[j] <- results_imputed$coefficient[4,1]
  se_imputed_fet[j] <- results_imputed$coefficient[4,2]
  pval_imputed_fet[j] <- results_imputed$coefficient[4,5]
  
  pval_imputed_omni[j] <- anova(lmer_imputed,lmer_null)[2,8]
  
  beta_sib[j] <- results_sib$coefficient[2,1]
  se_sib[j] <- results_sib$coefficient[2,2]
  pval_sib[j] <- results_sib$coefficient[2,5]
}