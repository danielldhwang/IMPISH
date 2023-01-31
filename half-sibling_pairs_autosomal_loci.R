# This is the R scripts used to imputed parental genotypes at autosomal loci using genotypes of offspring half-sibling pairs
# Please refer to Hwang et al. (2020) PLOS GENETICS

library(lmerTest)
library(MASS)
rm(list=(ls()))

Nrep = 1000	#Number of simulation replicates
N = 2000	#Sample size
Vm = 0.001	#Variance explained by maternal effect
Vp = 0.001	#Variance explained by paternal effect
Vf = 0.001	#Variance explained by fetal effect
p = 0.1		#Increaser allele frequency
Opp = FALSE	#Boolean variable to denote whether maternal and fetal effects in opposite directions (TRUE = opposite directions)
rho = 0.2	#Covariance between half sibling residual variances
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

#Analysis using only half sibling genotypes
beta_fet <- vector(length = Nrep)		#Fitting incorrect model including only half siblings
se_fet <- vector(length = Nrep)
pval_fet <- vector(length = Nrep)

#ANOVA for full model against null model
pval_genotyped_omni <- vector(length = Nrep)
pval_imputed_omni <- vector(length = Nrep)

a <- sqrt(1/2*p*q) 	#Create genetic variable of one. Assume no dominance.
BetaM <- sqrt(Vm) 	#Path coefficient for maternal effect
BetaP <- sqrt(Vp) 	#Path coefficient for paternal effect
BetaF <- sqrt(Vf) 	#Path coefficient for fetal effect
if(Opp == TRUE) {  	#If Opp is true then maternal and fetal effects have opposite directions of effect
  BetaF = -BetaF
}

Ve <- (1 - BetaM^2 - BetaP^2 - BetaF^2 - 2*0.5*BetaM*BetaF - 2*0.5*BetaP*BetaF) #Residual variance in trait for half siblings

for(j in 1:Nrep) {
  #Sample mothers' genotypes
  Zm <- sample(x = c('AA','Aa','aa'), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  #Sample fathers' genotypes
  Zp1 <- sample(x = c('AA','Aa','aa'), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  Zp2 <- sample(x = c('AA','Aa','aa'), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  
  Z_hsib1 <- vector(length = N)
  Z_hsib2 <- vector(length = N)
  Zm_imp <- vector(length = N)		#Imputed vector of genotypes at mum's locus
  Zp1_imp <- vector(length = N)		#Imputed vector of genotypes at dad1's locus
  Zp2_imp <- vector(length = N)		#Imputed vector of genotypes at dad2's locus
  
  #Simulate half sibling 1 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_hsib1[i] <- switch(paste(Zm[i],Zp1[i]),
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
  
  #Simulate half sibling 2 genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_hsib2[i] <- switch(paste(Zm[i],Zp2[i]),
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
  Z_hsib1_012 <- ifelse(Z_hsib1=='AA', 0, ifelse(Z_hsib1=='Aa', 1, 2))
  Z_hsib2_012 <- ifelse(Z_hsib2=='AA', 0, ifelse(Z_hsib2=='Aa', 1, 2))
  
  #Calculate observed allele frequency based on half sibling genotypes
  q_obs <- sum(Z_hsib1_012, Z_hsib2_012)/(2*N*2)
  p_obs <- 1 - q_obs
  
  #Impute mother's genotypes
  #Calculate conditional probability of mother's genotypes given half sibs genotypes
  p_AA_AA_AA = (2*p_obs)/(p_obs+1)				#P(Mum = AA | HSib = AA and HSib = AA) 
  p_Aa_AA_AA = (1-p_obs)/(p_obs+1)				#P(Mum = Aa | HSib = AA and HSib = AA)
  p_aa_AA_AA = 0							#P(Mum = aa | HSib = AA and HSib = AA)
  p_AA_AA_Aa = (2*p_obs)/(2*p_obs+1)				#P(Mum = AA | HSib = AA and HSib = Aa)
  p_Aa_AA_Aa = 1/(2*p_obs+1)					#P(Mum = Aa | HSib = AA and HSib = Aa)
  p_aa_AA_Aa = 0							#P(Mum = aa | HSib = AA and HSib = Aa)
  p_AA_AA_aa = 0							#P(Mum = AA | HSib = AA and HSib = aa)
  p_Aa_AA_aa = 1							#P(Mum = Aa | HSib = AA and HSib = aa)
  p_aa_AA_aa = 0							#P(Mum = aa | HSib = AA and HSib = aa)
  p_AA_Aa_Aa = ((2*p_obs)*(1-p_obs))/((4*p_obs)*(1-p_obs)+1)	#P(Mum = AA | HSib = Aa and HSib = Aa)
  p_Aa_Aa_Aa = (1-p_obs)/((p_obs-1)*(4*p_obs^2-4*p_obs-1))	#P(Mum = Aa | HSib = Aa and HSib = Aa)
  p_aa_Aa_Aa = ((2*q_obs)*(1-q_obs))/((4*q_obs)*(1-q_obs)+1)	#P(Mum = aa | HSib = Aa and HSib = Aa)
  p_AA_aa_Aa = 0							#P(Mum = AA | HSib = aa and HSib = Aa)
  p_Aa_aa_Aa = 1/(2*q_obs+1)					#P(Mum = Aa | HSib = aa and HSib = Aa)
  p_aa_aa_Aa = (2*q_obs)/(2*q_obs+1)				#P(Mum = aa | HSib = aa and HSib = Aa)
  p_AA_aa_aa = 0							#P(Mum = AA | HSib = aa and HSib = aa) 
  p_Aa_aa_aa = (1-q_obs)/(q_obs+1)					#P(Mum = Aa | HSib = aa and HSib = aa)
  p_aa_aa_aa = (2*q_obs)/(q_obs+1)					#P(Mum = aa | HSib = aa and HSib = aa)
  
  for (i in 1:N) {
    Zm_imp[i] <- switch(paste(Z_hsib1[i],Z_hsib2[i]),
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
  
  #Impute father's genotypes
  #Calculate conditional probability of father's genotypes given half sibling genotype 
  p_dAA_AA_AA <- p_obs 			#P(Dad = AA | HSib = AA and HSib = AA) 
  p_dAa_AA_AA <- q_obs  			#P(Dad = Aa | HSib = AA and HSib = AA) 
  p_daa_AA_AA <- 0  			#P(Dad = aa | HSib = AA and HSib = AA) 
  p_dAA_AA_Aa <- p_obs  			#P(Dad = AA | HSib = AA and HSib = Aa) 
  p_dAa_AA_Aa <- q_obs 			#P(Dad = Aa | HSib = AA and HSib = Aa) 
  p_daa_AA_Aa <- 0  				#P(Dad = aa | HSib = AA and HSib = Aa) 
  p_dAA_AA_aa <- p_obs 			#P(Dad = AA | HSib = AA and HSib = aa) 
  p_dAa_AA_aa <- q_obs 			#P(Dad = Aa | HSib = AA and HSib = aa) 
  p_daa_AA_aa <- 0 				#P(Dad = aa | HSib = AA and HSib = aa)   
  p_dAA_Aa_AA <- p_obs^2/(2*p_obs+1)  	#P(Dad = AA | HSib = Aa and HSib = AA) 
  p_dAa_Aa_AA <- 2*p_obs/(2*p_obs+1)  	#P(Dad = Aa | HSib = Aa and HSib = AA) 
  p_daa_Aa_AA <- (1-p_obs^2)/(2*p_obs+1) 	#P(Dad = aa | HSib = Aa and HSib = AA)  
  p_dAA_dAA_Aa_Aa <- (2*p_obs^3-p_obs^4)/(1+4*p_obs-4*p_obs^2)	#P(Dad1 = AA and Dad2 = AA | HSib1 = Aa and HSib2 = Aa) 
  p_dAA_dAa_Aa_Aa <- (2*p_obs^2-2*p_obs^3)/(1+4*p_obs-4*p_obs^2)	#P(Dad1 = AA and Dad2 = Aa | HSib1 = Aa and HSib2 = Aa)
  p_dAA_daa_Aa_Aa <- (p_obs^2*q_obs^2)/(1+4*p_obs-4*p_obs^2)	#P(Dad1 = AA and Dad2 = aa | HSib1 = Aa and HSib2 = Aa)
  p_dAa_dAA_Aa_Aa <- p_dAA_dAa_Aa_Aa 			 #P(Dad1 = Aa and Dad2 = AA | HSib1 = Aa and HSib2 = Aa)
  p_dAa_dAa_Aa_Aa <- (2*p_obs*q_obs)/(1+4*p_obs-4*p_obs^2) #P(Dad1 = Aa and Dad2 = Aa | HSib1 = Aa and HSib2 = Aa)
  p_dAa_daa_Aa_Aa <- (2*q_obs^2-2*q_obs^3)/(1+4*q_obs-4*q_obs^2)	#P(Dad1 = Aa and Dad2 = aa | HSib1 = Aa and HSib2 = Aa)
  p_daa_dAA_Aa_Aa <- p_dAA_daa_Aa_Aa  			#P(Dad1 = aa and Dad2 = AA | HSib1 = Aa and HSib2 = Aa)
  p_daa_dAa_Aa_Aa <- p_dAa_daa_Aa_Aa  			#P(Dad1 = aa and Dad2 = Aa | HSib1 = Aa and HSib2 = Aa)
  p_daa_daa_Aa_Aa <- (2*q_obs^3-q_obs^4)/(1+4*q_obs-4*q_obs^2)  #P(Dad1 = aa and Dad2 = aa | HSib1 = Aa and HSib2 = Aa)
  p_dAA_Aa_aa <- p_obs*(2-p_obs)/(3-2*p_obs) 			#P(Dad = AA | HSib = Aa and HSib = aa) 
  p_dAa_Aa_aa <- (2-2*p_obs)/(3-2*p_obs)				#P(Dad = Aa | HSib = Aa and HSib = aa) 
  p_daa_Aa_aa <- (1-p_obs)^2/(3-2*p_obs) 				#P(Dad = aa | HSib = Aa and HSib = aa) 
  p_dAA_aa_AA <- 0 							#P(Dad = AA | HSib = aa and HSib = AA) 
  p_dAa_aa_AA <- p_obs 						#P(Dad = Aa | HSib = aa and HSib = AA) 
  p_daa_aa_AA <- q_obs 						#P(Dad = aa | HSib = aa and HSib = AA) 
  p_dAA_aa_Aa <- 0 							#P(Dad = AA | HSib = aa and HSib = Aa) 
  p_dAa_aa_Aa <- p_obs 						#P(Dad = Aa | HSib = aa and HSib = Aa)
  p_daa_aa_Aa <- q_obs 						#P(Dad = aa | HSib = aa and HSib = Aa)
  
  p_dAA_aa_aa <- 0 							#P(Dad = AA | HSib = aa and HSib = aa)
  p_dAa_aa_aa <- p_obs 						#P(Dad = Aa | HSib = aa and HSib = aa)
  p_daa_aa_aa <- q_obs 						#P(Dad = aa | HSib = aa and HSib = aa)
  
  for (i in 1:N) {
    Zp1_imp[i] <- switch(paste(Z_hsib1[i],Z_hsib2[i]),
                         'AA AA'= -a*p_dAA_AA_AA + 0*p_dAa_AA_AA + a*p_daa_AA_AA,
                         'AA Aa'= -a*p_dAA_AA_Aa + 0*p_dAa_AA_Aa + a*p_daa_AA_Aa,
                         'AA aa'= -a*p_dAA_AA_aa + 0*p_dAa_AA_aa + a*p_daa_AA_aa,
                         'Aa AA'= -a*p_dAA_Aa_AA + 0*p_dAa_Aa_AA + a*p_daa_Aa_AA,
                         'Aa Aa'= -a*(p_dAA_dAA_Aa_Aa + p_dAA_dAa_Aa_Aa + p_dAA_daa_Aa_Aa) + 
                           0*(p_dAa_dAA_Aa_Aa + p_dAa_dAa_Aa_Aa + p_dAa_daa_Aa_Aa) +
                           a*(p_daa_dAA_Aa_Aa + p_daa_dAa_Aa_Aa + p_daa_daa_Aa_Aa),
                         'Aa aa'= -a*p_dAA_Aa_aa + 0*p_dAa_Aa_aa + a*p_daa_Aa_aa,
                         'aa AA'= -a*p_dAA_aa_AA + 0*p_dAa_aa_AA + a*p_daa_aa_AA,
                         'aa Aa'= -a*p_dAA_aa_Aa + 0*p_dAa_aa_Aa + a*p_daa_aa_Aa,
                         'aa aa'= -a*p_dAA_aa_aa + 0*p_dAa_aa_aa + a*p_daa_aa_aa)
    Zp2_imp[i] <- switch(paste(Z_hsib2[i],Z_hsib1[i]),
                         'AA AA'= -a*p_dAA_AA_AA + 0*p_dAa_AA_AA + a*p_daa_AA_AA,
                         'AA Aa'= -a*p_dAA_AA_Aa + 0*p_dAa_AA_Aa + a*p_daa_AA_Aa,
                         'AA aa'= -a*p_dAA_AA_aa + 0*p_dAa_AA_aa + a*p_daa_AA_aa,
                         'Aa AA'= -a*p_dAA_Aa_AA + 0*p_dAa_Aa_AA + a*p_daa_Aa_AA,
                         'Aa Aa'= -a*(p_dAA_dAA_Aa_Aa + p_dAA_dAa_Aa_Aa + p_dAA_daa_Aa_Aa) + 
                           0*(p_dAa_dAA_Aa_Aa + p_dAa_dAa_Aa_Aa + p_dAa_daa_Aa_Aa) +
                           a*(p_daa_dAA_Aa_Aa + p_daa_dAa_Aa_Aa + p_daa_daa_Aa_Aa),
                         'Aa aa'= -a*p_dAA_Aa_aa + 0*p_dAa_Aa_aa + a*p_daa_Aa_aa,
                         'aa AA'= -a*p_dAA_aa_AA + 0*p_dAa_aa_AA + a*p_daa_aa_AA,
                         'aa Aa'= -a*p_dAA_aa_Aa + 0*p_dAa_aa_Aa + a*p_daa_aa_Aa,
                         'aa aa'= -a*p_dAA_aa_aa + 0*p_dAa_aa_aa + a*p_daa_aa_aa)
  }  
  
  #Convert AA/Aa/aa to genetic value -a/0/a
  Zm <- ifelse(Zm=='AA', -a, ifelse(Zm=='Aa', 0, a))
  Zp1 <- ifelse(Zp1=='AA', -a, ifelse(Zp1=='Aa', 0, a))
  Zp2 <- ifelse(Zp2=='AA', -a, ifelse(Zp2=='Aa', 0, a))
  Z_hsib1 <- ifelse(Z_hsib1=='AA', -a, ifelse(Z_hsib1=='Aa', 0, a))
  Z_hsib2 <- ifelse(Z_hsib2=='AA', -a, ifelse(Z_hsib2=='Aa', 0, a))
  
  #Create correlated error variables for sib 1 and sib 2 in the model with both parents 
  Sigma <- matrix(c(Ve, rho, rho, Ve),2,2)
  e <- mvrnorm(n = N, mu = c(0, 0), Sigma)
  
  #Simulate offspring outcome
  Y_hsib1 <- BetaM*Zm + BetaP*Zp1 + BetaF*Z_hsib1 + e[,1]
  Y_hsib2 <- BetaM*Zm + BetaP*Zp2 + BetaF*Z_hsib2 + e[,2]
  
  test <- data.frame(rbind(cbind(1:N, Y_hsib1, Zm, Zp1, Zm_imp, Zp1_imp, Z_hsib1), 
                           cbind(1:N, Y_hsib2, Zm, Zp2, Zm_imp, Zp2_imp, Z_hsib2)))
  colnames(test) <- c("fam", "Y","Zm", "Zp", "Zm_imp", "Zp_imp", "Z")
  
  #Run analyses
  lmer_genotyped <- lmer(Y ~ Zm + Zp + Z + (1|fam), data = test, REML = FALSE);  results_genotyped <- summary(lmer_genotyped)
  lmer_imputed <- lmer(Y ~ Zm_imp + Zp_imp + Z + (1|fam), data = test, REML = FALSE); results_imputed <- summary(lmer_imputed)
  lmer_hsib <- lmer(Y ~ Z + (1|fam), data = test, REML = FALSE); results_hsib <- summary(lmer_hsib)
  lmer_null <- lmer(Y ~ (1|fam), data = test, REML = FALSE)
  
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
  
  beta_fet[j] <- results_hsib$coefficient[2,1]
  se_fet[j] <- results_hsib$coefficient[2,2]
  pval_fet[j] <- results_hsib$coefficient[2,5]
}