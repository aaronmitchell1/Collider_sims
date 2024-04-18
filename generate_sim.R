## Auxiliary function to simulate logistic data.
expit <- function (x) exp(x) / (1 + exp(x))

## Set up the simulation - set sample size and common cause.
set.seed(54545)
n <- 1e5
U <- rnorm(n, 0, 1)

##Simulate one-sample genetic data from individual-level data - number of SNPs, simulate MAF, ratios for how many affect I and P.
nSNPs <- 100
maf <- runif(nSNPs, 0.05, 0.5)
incidence.SNPs <- 1:90
progression.SNPs <- 91:100
G <- matrix(0, n, nSNPs)
for (j in 1:nSNPs) G[, j] <- rbinom(n, 2, maf[j])
incidence.betas <- rep(0, nSNPs)
incidence.betas[incidence.SNPs] <- rnorm(length(incidence.SNPs), 0, 0.2)
incidence.probs <- expit(-1 + as.vector(G %*% incidence.betas) + 1 * U)
I <- rbinom(n, 1, incidence.probs)
progression.betas <- rep(0, nSNPs)
progression.betas[progression.SNPs] <- rnorm(length(progression.SNPs), 0, 0.2)
progression.probs <- expit(-1 + as.vector(G %*% progression.betas) + 1 * U)
P <- rbinom(n, 1, progression.probs)

##Simulate GWAS of incidence.
#Include effects of progression SNPs on incidence and vice versa to conduct subsequent analyses.

incidence_GWAS <- matrix(0, nSNPs, 3); colnames(incidence_GWAS) <- c("Estimate", "StdErr", "Pval")

for (j in incidence.SNPs) {
  incidence_model <- glm(I ~ G[,j], family = binomial)
  incidence_GWAS[j, 1] <- summary(incidence_model)$coefficients["G[, j]", "Estimate"]
  incidence_GWAS[j, 2] <- summary(incidence_model)$coefficients["G[, j]", "Std. Error"]
  incidence_GWAS[j, 3] <- summary(incidence_model)$coefficients["G[, j]", "Pr(>|z|)"]
}

for (j in progression.SNPs) {
  incidence_model <- glm(I ~ G[,j], family = binomial)
  incidence_GWAS[j, 1] <- summary(incidence_model)$coefficients["G[, j]", "Estimate"]
  incidence_GWAS[j, 2] <- summary(incidence_model)$coefficients["G[, j]", "Std. Error"]
  incidence_GWAS[j, 3] <- summary(incidence_model)$coefficients["G[, j]", "Pr(>|z|)"]
}

##Simulate GWAS of progression.

progression_GWAS <- matrix(0, nSNPs, 3); colnames(progression_GWAS) <- c("Estimate", "StdErr", "Pval")

for (j in progression.SNPs) {
  progression_model <- glm(P ~ G[,j], family = binomial)
  progression_GWAS[j, 1] <- summary(progression_model)$coefficients["G[, j]", "Estimate"]
  progression_GWAS[j, 2] <- summary(progression_model)$coefficients["G[, j]", "Std. Error"]
  progression_GWAS[j, 3] <- summary(progression_model)$coefficients["G[, j]", "Pr(>|z|)"]
}

for (j in incidence.SNPs) {
  progression_model <- glm(P ~ G[,j], family = binomial)
  progression_GWAS[j, 1] <- summary(progression_model)$coefficients["G[, j]", "Estimate"]
  progression_GWAS[j, 2] <- summary(progression_model)$coefficients["G[, j]", "Std. Error"]
  progression_GWAS[j, 3] <- summary(progression_model)$coefficients["G[, j]", "Pr(>|z|)"]
}
