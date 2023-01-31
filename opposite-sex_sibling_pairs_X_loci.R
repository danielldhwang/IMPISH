# This is the R scripts used to imputed parental genotypes at X chromosomal loci using genotypes of offspring male-female sibling pairs (opposite sex)
# Please refer to Hwang et al. (2020) PLOS GENETICS

library(MASS)
library(lmerTest)
rm(list=ls())

Nrep = 1000 	#Number of simulation replicates
N = 2000 	#Sample size
Vm = 0.001 	#Variance explained by maternal effect
Vp = 0.001 	#Variance explained by paternal effect
Vf = 0.001 	#Variance explained by fetal effect
p = 0.5 	#Increaser allele frequency
Opp = FALSE #Boolean variable to denote whether maternal and fetal effects in opposite directions (TRUE = opposite directions)
rho = 0.2 	#Covariance between sibling residual variances
q <- 1-p 	#Decreaser allele frequency

#Create vectors to save results from regression analyses
#For analysis using simulated genotypes:
beta_genotyped_mat <- vector(length = Nrep) #Fitting regression with maternal and paternal genotypes known (maternal coefficient)
beta_genotyped_pat <- vector(length = Nrep)  #Fitting regression with maternal and paternal genotypes known (paternal coefficient)
beta_genotyped_fet <- vector(length = Nrep)   #Fitting regression with maternal and paternal genotypes known (fetal coefficient)
se_genotyped_mat <- vector(length = Nrep)
se_genotyped_pat <- vector(length = Nrep)
se_genotyped_fet <- vector(length = Nrep)
pval_genotyped_mat <- vector(length = Nrep)
pval_genotyped_pat <- vector(length = Nrep)
pval_genotyped_fet <- vector(length = Nrep)

#Analysis using imputed maternal and paternal genotypes: 
beta_imputed_mat <- vector(length = Nrep) #Fitting regression with genotypes imputed (maternal coefficient)
beta_imputed_pat <- vector(length = Nrep)  #Fitting regression with genotypes imputed (paternal coefficient)
beta_imputed_fet <- vector(length = Nrep)   #Fitting regression with genotypes imputed (fetal coefficient)
se_imputed_mat <- vector(length = Nrep)
se_imputed_pat <- vector(length = Nrep)
se_imputed_fet <- vector(length = Nrep)
pval_imputed_mat <- vector(length = Nrep)
pval_imputed_pat <- vector(length = Nrep)
pval_imputed_fet <- vector(length = Nrep)

#Analysis using only half-sib genotypes: 
beta_fet <- vector(length = Nrep) #Fitting incorrect model to just sibs
se_fet <- vector(length = Nrep)
pval_fet <- vector(length = Nrep)

#ANOVA for full model against null model
pval_imputed_omni <- vector(length = Nrep)
pval_genotyped_omni <- vector(length = Nrep)

a <- sqrt(1/(2*p*q)) 		#Create genetic variable of variance one for maternal SNPs on X chromosome. Assume no dominance.
BetaM <- sqrt(Vm) 		#Path coefficient for maternal effect
BetaP <- sqrt(Vp/2) 		#Path coefficient for paternal effect
BetaF <- sqrt(Vf/2) 		#Path coefficient for fetal effect. Assume **male** genotype
if(Opp == TRUE) {  		#If Opp is true then maternal and fetal effects have opposite directions of effect
  BetaF = -BetaF
}

Ve_x <- (1 - BetaM^2 - BetaF^2 - BetaP^2*2 - 2*0.5*BetaM*BetaF - 2*0.5*BetaP*BetaF*2) #Residual variance in trait for female sib
Ve_y <- (1 - BetaM^2 - BetaF^2*2 - BetaP^2*2 - 2*1*BetaM*BetaF) #Residual variance in trait for male sib

for(j in 1:Nrep) {
  #Sample mothers' genotypes
  Zm <- sample(x = c('AA','Aa','aa'), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
  #Sample fathers' genotypes
  Zp <- sample(x = c('A','a'), size = N, replace = TRUE, prob = c(p, q))
  
  Z_sibx <- vector(length = N)
  Z_siby <- vector(length = N)
  
  Zmxy_imp <- vector(length = N) #Imputed vector of genotypes at mum's locus
  Zpxy_imp <- vector(length = N) #Imputed vector of genotypes at dad's locus
  
  #Simulate female sib genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_sibx[i] <- switch(paste(Zm[i],Zp[i]),
                        'AA A'='AA',
                        'AA a'='Aa',
                        'Aa A'=ifelse(r[i] <= 0.5, 'AA', 'Aa'),
                        'Aa a'=ifelse(r[i] <= 0.5, 'Aa', 'aa'),
                        'aa A'='Aa',
                        'aa a'='aa')
  }
  
  #Simulate male sib genotype
  r <- runif(N)
  for (i in 1:N) { 
    Z_siby[i] <- switch(Zm[i],
                        'AA'='A',
                        'Aa'=ifelse(r[i] <= 0.5, 'A', 'a'),
                        'aa'='a')
  }
  
  #Change coding to 0, 1 and 2 in order to calculate allele frequencies
  Z_sibx_012 <- ifelse(Z_sibx=='AA', 0, ifelse(Z_sibx=='Aa', 1, 2))
  Z_siby_01 <- ifelse(Z_siby=='A', 0, 1)
  
  #Calculate observed allele frequency based on female sibling one's genotypes
  q_obs <- sum(Z_sibx_012)/(2*N)
  p_obs <- 1 - q_obs
  
  #Impute parental genotypes from male-female siblings
  #Calculate conditional probability of mother's genotypes given sibling genotypes
  p_AA_A_AA = (2*p_obs)/(p_obs+1)  	#P(Mum = AA | Siby = A and Sibx = AA)
  p_Aa_A_AA = (1-p_obs)/(p_obs+1)  	#P(Mum = Aa | Siby = A and Sibx = AA)
  p_aa_A_AA = 0 				#P(Mum = aa | Siby = A and Sibx = AA)
  p_AA_A_Aa = (2*p_obs)/(2*p_obs+1) 	#P(Mum = AA | Siby = A and Sibx = Aa) 
  p_Aa_A_Aa = 1/(2*p_obs+1)        		#P(Mum = Aa | Siby = A and Sibx = Aa) 
  p_aa_A_Aa = 0 				#P(Mum = aa | Siby = A and Sibx = Aa) 
  p_AA_A_aa = 0 				#P(Mum = AA | Siby = A and Sibx = aa) 
  p_Aa_A_aa = 1 				#P(Mum = Aa | Siby = A and Sibx = aa) 
  p_aa_A_aa = 0 				#P(Mum = aa | Siby = A and Sibx = aa) 
  
  p_AA_a_AA = 0 				#P(Mum = AA | Siby = a and Sibx = AA) 
  p_Aa_a_AA = 1 				#P(Mum = Aa | Siby = a and Sibx = AA) 
  p_aa_a_AA = 0 				#P(Mum = aa | Siby = a and Sibx = AA) 
  p_AA_a_Aa = 0 				#P(Mum = AA | Siby = a and Sibx = Aa) 
  p_Aa_a_Aa = 1/(2*q_obs+1)  		#P(Mum = Aa | Siby = a and Sibx = Aa) 
  p_aa_a_Aa = (2*q_obs)/(2*q_obs+1)  	#P(Mum = aa | Siby = a and Sibx = Aa) 
  p_AA_a_aa = 0 				#P(Mum = AA | Siby = a and Sibx = aa) 
  p_Aa_a_aa = (1-q_obs)/(q_obs+1) 		#P(Mum = Aa | Siby = a and Sibx = aa)
  p_aa_a_aa = (2*q_obs)/(q_obs+1) 		#P(Mum = aa | Siby = a and Sibx = aa)
  
  p_AA_AA_A <- p_AA_A_AA
  p_Aa_AA_A <- p_Aa_A_AA 
  p_aa_AA_A <- p_aa_A_AA
  p_AA_Aa_A <- p_AA_A_Aa 
  p_Aa_Aa_A <- p_Aa_A_Aa 
  p_aa_Aa_A <- p_aa_A_Aa
  p_AA_aa_A <- p_AA_A_aa
  p_Aa_aa_A <- p_Aa_A_aa
  p_aa_aa_A <- p_aa_A_aa
  
  p_AA_AA_a <- p_AA_a_AA
  p_Aa_AA_a <- p_Aa_a_AA
  p_aa_AA_a <- p_aa_a_AA
  p_AA_Aa_a <- p_AA_a_Aa
  p_Aa_Aa_a <- p_Aa_a_Aa
  p_aa_Aa_a <- p_aa_a_Aa
  p_AA_aa_a <- p_AA_a_aa
  p_Aa_aa_a <- p_Aa_a_aa
  p_aa_aa_a <- p_aa_a_aa
  
  for (i in 1:N) {
    Zmxy_imp[i] <- switch(paste(Z_siby[i], Z_sibx[i]),
                          'A AA'= -a*p_AA_A_AA + 0*p_Aa_A_AA + a*p_aa_A_AA,
                          'A Aa'= -a*p_AA_A_Aa + 0*p_Aa_A_Aa + a*p_aa_A_Aa,
                          'A aa'= -a*p_AA_A_aa + 0*p_Aa_A_aa + a*p_aa_A_aa,
                          'a AA'= -a*p_AA_a_AA + 0*p_Aa_a_AA + a*p_aa_a_AA,
                          'a Aa'= -a*p_AA_a_Aa + 0*p_Aa_a_Aa + a*p_aa_a_Aa,
                          'a aa'= -a*p_AA_a_aa + 0*p_Aa_a_aa + a*p_aa_a_aa)
  }
  
  #Calculate conditional probability of father's genotypes given sib's genotype
  p_A_A_AA <- 1 				#P(Dad = A | Siby = A and Sibx = AA)
  p_a_A_AA <- 0 				#P(Dad = a | Siby = A and Sibx = AA)
  p_A_A_Aa <- p_obs/(2*p_obs+1) 		#P(Dad = A | Siby = A and Sibx = Aa)
  p_a_A_Aa <- (p_obs+1)/(2*p_obs+1) 	#P(Dad = a | Siby = A and Sibx = Aa)
  p_A_A_aa <- 0 				#P(Dad = A | Siby = A and Sibx = aa)
  p_a_A_aa <- 1 				#P(Dad = a | Siby = A and Sibx = aa)
  p_A_a_AA <- 1 				#P(Dad = A | Siby = a and Sibx = AA)
  p_a_a_AA <- 0 				#P(Dad = a | Siby = a and Sibx = AA)
  p_A_a_Aa <- (q_obs+1)/(2*q_obs+1) 	#P(Dad = A | Siby = a and Sibx = Aa)
  p_a_a_Aa <- q_obs/(2*q_obs+1) 		#P(Dad = a | Siby = a and Sibx = Aa)
  p_A_a_aa <- 0 				#P(Dad = A | Siby = a and Sibx = aa)
  p_a_a_aa <- 1 				#P(Dad = a | Siby = a and Sibx = aa)
  
  for (i in 1:N) {
    Zpxy_imp[i] <- switch(paste(Z_siby[i], Z_sibx[i]),
                          'A AA'= -a*p_A_A_AA + a*p_a_A_AA,
                          'A Aa'= -a*p_A_A_Aa + a*p_a_A_Aa,
                          'A aa'= -a*p_A_A_aa + a*p_a_A_aa,
                          'a AA'= -a*p_A_a_AA + a*p_a_a_AA,
                          'a Aa'= -a*p_A_a_Aa + a*p_a_a_Aa,
                          'a aa'= -a*p_A_a_aa + a*p_a_a_aa)
  }    
  
  #Convert AA/Aa/aa to genetic value -a/0/a
  Zm <- ifelse(Zm=='AA', -a, ifelse(Zm=='Aa', 0, a))
  Zp <- ifelse(Zp=='A', -a, a)
  Z_sibx <- ifelse(Z_sibx=='AA', -a, ifelse(Z_sibx=='Aa', 0, a))
  Z_siby <- ifelse(Z_siby=='A', -a, a)
  
  #Create correlated error variables for sib 1 and sib 2 in the model
  Sigma_xy <- matrix(c(Ve_x, rho, rho, Ve_y),2,2)  #Female-male sibs
  e_xy <- mvrnorm(n = N, mu = c(0, 0), Sigma_xy)
  
  #Simulate offspring outcome
  Y_sibx_xy <- BetaM*Zm + BetaP*Zp + BetaF*Z_sibx + e_xy[,1]
  Y_siby_xy <- BetaM*Zm + BetaP*Zp + BetaF*Z_siby + e_xy[,2]
  
  testxy <- data.frame(rbind(cbind(1:N, Y_sibx_xy, Zm, Zp, Zmxy_imp, Zpxy_imp, Z_sibx), 
                             cbind(1:N, Y_siby_xy, Zm, Zp, Zmxy_imp, Zpxy_imp, Z_siby)))
  
  colnames(testxy) <- c("fam", "Y_xy", "Zm", "Zp", "Zmxy_imp", "Zpxy_imp", "Z_xy")
  
  #Run analyses
  lmer_genotyped <- lmer(Y_xy ~ Zm + Zp + Z_xy + (1|fam), data = testxy, REML = FALSE);  results_genotyped <- summary(lmer_genotyped)
  lmer_imputed <- lmer(Y_xy ~ Zmxy_imp + Zpxy_imp + Z_xy + (1|fam), data = testxy, REML = FALSE); results_imputed <- summary(lmer_imputed)
  lmer_sib <- lmer(Y_xy ~ Z_xy + (1|fam), data = testxy, REML = FALSE); results_sib <- summary(lmer_sib)
  lmer_null <- lmer(Y_xy ~ (1|fam), data = testxy, REML = FALSE)
  
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
  
  beta_fet[j] <- results_sib$coefficient[2,1]
  se_fet[j] <- results_sib$coefficient[2,2]
  pval_fet[j] <- results_sib$coefficient[2,5]
}
