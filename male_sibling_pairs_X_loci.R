# This is the R scripts used to imputed parental genotypes at X chromosomal loci using genotypes of offspring male sibling pairs
# Please refer to Hwang et al. (2020) PLOS GENETICS

library(lmerTest)
library(MASS)
rm(list=ls())

Nrep = 1000	#Number of simulation replicates
N = 2000	#Sample size
Vm = 0.001	#Variance explained by maternal effect
Vp = 0.00	#Variance explained by paternal effect
Vf = 0.001	#Variance explained by fetal effect
p = 0.1		#Increaser allele frequency
Opp = FALSE	#Boolean variable to denote whether maternal and fetal effects in opposite directions (TRUE = opposite directions)
rho = 0.2	#Covariance between sibling residual variances
q <- 1-p 	#Decreaser allele frequency

#Create vectors to save results from regression analyses
#For analysis using simulated genotypes:
beta_genotyped_mat <- vector(length = Nrep)		#Fitting regression with maternal and paternal genotypes known (maternal coefficient)
beta_genotyped_fet <- vector(length = Nrep)	  #Fitting regression with maternal and paternal genotypes known (fetal coefficient)
se_genotyped_mat <- vector(length = Nrep)
se_genotyped_fet <- vector(length = Nrep)
pval_genotyped_mat <- vector(length = Nrep)
pval_genotyped_fet <- vector(length = Nrep)

#Analysis using imputed maternal and paternal genotypes: 
beta_imputed_mat <- vector(length = Nrep)	#Fitting regression with genotypes imputed (maternal coefficient)
beta_imputed_fet <- vector(length = Nrep)	#Fitting regression with genotypes imputed (fetal coefficient)
se_imputed_mat <- vector(length = Nrep)
se_imputed_fet <- vector(length = Nrep)
pval_imputed_mat <- vector(length = Nrep)
pval_imputed_fet <- vector(length = Nrep)

#Analysis using only sibling genotypes: 
beta_sib <- vector(length = Nrep)		#Fitting incorrect model including only siblings
se_sib <- vector(length = Nrep)
pval_sib <- vector(length = Nrep)

#ANOVA for full model against null model
pval_genotyped_omni <- vector(length = Nrep)
pval_imputed_omni <- vector(length = Nrep)


a <- sqrt(1/(2*p*q))  	#Create genetic variable of variance one for maternal SNPs on X chromosome. Assume no dominance.
BetaM <- sqrt(Vm) 	#Path coefficient for maternal effect
BetaP <- sqrt(Vp/2) 	#Path coefficient for paternal effect
BetaF <- sqrt(Vf/2) 	#Path coefficient for fetal effect
if(Opp == TRUE) {  	#If Opp is true then maternal and fetal effects have opposite directions of effect
  BetaF = -BetaF
}

Ve_y <- (1 - BetaM^2 - BetaF^2*2 - BetaP^2*2 - 2*1*BetaM*BetaF) #Residual variance in trait for male sib

for(j in 1:Nrep) {
  #Sample mothers' genotypes
  Zm <- sample(x = c('AA','Aa','aa'), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  #Sample fathers' genotypes
  Zp <- sample(x = c('A','a'), size = N, replace = TRUE, prob = c(p, q))
  
  Z_sib1y <- vector(length = N)
  Z_sib2y <- vector(length = N)
  Zmyy_imp <- vector(length = N)		#Imputed vector of genotypes at mum's locus
  
  #Simulate male sib1 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sib1y[i] <- switch(Zm[i],
                         'AA'='A',
                         'Aa'=ifelse(r[i] <= 0.5, 'A', 'a'),
                         'aa'='a')
  }
  
  #Simulate male sib2 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sib2y[i] <- switch(Zm[i],
                         'AA'='A',
                         'Aa'=ifelse(r[i] <= 0.5, 'A', 'a'),
                         'aa'='a')
  }
  
  #Change coding to 0, 1 and 2 in order to calculate allele frequencies
  Z_sib1y_01 <- ifelse(Z_sib1y=='A', 0, 2)
  Z_sib2y_01 <- ifelse(Z_sib2y=='A', 0, 2)
  
  #Calculate observed allele frequency based on sibling one's genotypes
  q_obs <- sum(c(Z_sib1y_01, Z_sib2y_01))/(2*N*2)
  p_obs <- 1 - q_obs
  
  #Impute parental genotypes from male-female sibs
  #Calculate conditional probability of mother's genotypes given sibs genotypes
  p_AA_A_A = (2*p_obs)/(p_obs+1)	  	#P(Mum = AA | Sib1y = A and Sib2y = A)
  p_Aa_A_A = (1-p_obs)/(p_obs+1)	  	#P(Mum = Aa | Sib1y = A and Sib2y = A)
  p_aa_A_A = 0                      			#P(Mum = aa | Sib1y = A and Sib2y = A)
  p_AA_A_a = 0                     			#P(Mum = AA | Sib1y = A and Sib2y = a) 
  p_Aa_A_a = 1            	        		#P(Mum = Aa | Sib1y = A and Sib2y = a) 
  p_aa_A_a = 0            	        		#P(Mum = aa | Sib1y = A and Sib2y = a) 
  p_AA_a_A = p_AA_A_a            	        
  p_Aa_a_A = p_Aa_A_a
  p_aa_a_A = p_aa_A_a
  p_AA_a_a = 0                      			#P(Mum = AA | Sib1y = a and Sib2y = a) 
  p_Aa_a_a = (1-q_obs)/(q_obs+1)    		#P(Mum = Aa | Sib1y = a and Sib2y = a) 
  p_aa_a_a = (2*q_obs)/(q_obs+1)    		#P(Mum = aa | Sib1y = a and Sib2y = a) 
  
  for (i in 1:N) {
    Zmyy_imp[i] <- switch(paste(Z_sib1y[i],Z_sib2y[i]),
                          'A A'= -a*p_AA_A_A + 0*p_Aa_A_A + a*p_aa_A_A,
                          'A a'= -a*p_AA_A_a + 0*p_Aa_A_a + a*p_aa_A_a,
                          'a A'= -a*p_AA_a_A + 0*p_Aa_a_A + a*p_aa_a_A,
                          'a a'= -a*p_AA_a_a + 0*p_Aa_a_a + a*p_aa_a_a)
  }
  
  #Convert AA/Aa/aa to genetic value
  Zm <- ifelse(Zm=='AA', -a, ifelse(Zm=='Aa', 0, a))
  Zp <- ifelse(Zp=='A', -a, a)
  Z_sib1y <- ifelse(Z_sib1y=='A', -a, a)
  Z_sib2y <- ifelse(Z_sib2y=='A', -a, a)
  
  #Create correlated error variables for sib 1 and sib 2 in the model
  Sigma_yy <- matrix(c(Ve_y, rho, rho, Ve_y),2,2)  #Male-male sibs
  e_yy <- mvrnorm(n = N, mu = c(0, 0), Sigma_yy)
  
  #Simulate offspring outcome
  Y_sib1y_yy <- BetaM*Zm + BetaP*Zp + BetaF*Z_sib1y + e_yy[,1]
  Y_sib2y_yy <- BetaM*Zm + BetaP*Zp + BetaF*Z_sib2y + e_yy[,2]
  
  test <- data.frame(rbind(cbind(1:N, Y_sib1y_yy, Zm, Zp, Zmyy_imp, Z_sib1y),
                           cbind(1:N, Y_sib2y_yy, Zm, Zp, Zmyy_imp, Z_sib2y)))
  colnames(test) <- c("fam", "Y", "Zm", "Zp", "Zmyy_imp", "Zy")
  
  #Run analyses
  lmer_genotyped <- lmer(Y ~ Zm + Zy + (1|fam), data = test, REML = FALSE);  results_genotyped <- summary(lmer_genotyped)
  lmer_imputed <- lmer(Y ~ Zmyy_imp + Zy + (1|fam), data = test, REML = FALSE); results_imputed <- summary(lmer_imputed)
  lmer_sib <- lmer(Y ~ Zy + (1|fam), data = test, REML = FALSE); results_sib <- summary(lmer_sib)
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