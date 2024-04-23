##Diagnostics

##R^2 from linear regression of G on I/P as suggested by Eleanor Sanderson, this is a better method than a logistic approach.
##Not really necessary to have this on a per-SNP basis but if needed easy to implement using subset, as with main sim.

##How much variation is explained in incidence overall?

incidence_rsq <- summary(lm(I ~ G))$r.squared

##Progression

progression_rsq <- summary(lm(P ~ G))$r.squared

##Repetitions (nsims) For loop 

##Means for I and P

incidence_mean <- matrix(0, nSNPs, 1); colnames(incidence_mean) <- "Mean"

for (j in incidence.SNPs) {
incidence_mean[j, 1] (I[,j])

  }

for (j in progression.SNPs) {
incidence_mean[j, 1] (I[,j])

  }

progression_mean <- matrix(0, nSNPs, 1); colnames(progression_mean) <- "Mean"

for (j in progression.SNPs) {
progression_mean[j, 1] (P[,j])

  }

for (j in incidence.SNPs) {
progression_mean[j, 1] (P[,j])

  }

##Monte-Carlo SE

##Confidence interval coverage
