##Diagnostics

##R^2 for each SNP from linear regression of G on I/P as suggested by Eleanor Sanderson.

##Incidence

incidence_rsq <- matrix(0, nSNPs, 1); colnames(incidence_rsq) <- "Rsq"

for (j in incidence.SNPs) {
incidence_rsq_model <- lm(I ~ G[,j])
incidence_rsq[j, 1] <- summary(incidence_rsq_model)$r.squared
}

for (j in progression.SNPs) {
incidence_rsq_model <- lm(I ~ G[,j])
incidence_rsq[j, 1] <- summary(incidence_rsq_model)$r.squared
}

##Progression

progression_rsq <- matrix(0, nSNPs, 1); colnames(progression_rsq) <- "Rsq"

for (j in progression.SNPs) {
  progression_rsq_model <- lm(P ~ G[,j])
  progression_rsq[j, 1] <- summary(progression_rsq_model)$r.squared
}

for (j in incidence.SNPs) {
  progression_model <- glm(P ~ G[,j], family = binomial)
  progression_GWAS[j, 1] <- summary(progression_model)$coefficients["G[, j]", "Estimate"]
  progression_GWAS[j, 2] <- summary(progression_model)$coefficients["G[, j]", "Std. Error"]
  progression_GWAS[j, 3] <- summary(progression_model)$coefficients["G[, j]", "Pr(>|z|)"]
}

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
